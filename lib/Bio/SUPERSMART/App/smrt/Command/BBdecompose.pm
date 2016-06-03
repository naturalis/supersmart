package Bio::SUPERSMART::App::smrt::Command::BBdecompose;

use strict;
use warnings;

use List::MoreUtils 'uniq';
use File::Spec;

use Bio::Phylo::IO 'parse_tree';
use Bio::SUPERSMART::Config;
use Bio::SUPERSMART::Domain::MarkersAndTaxa;
use Bio::SUPERSMART::Service::MarkersAndTaxaSelector;
use Bio::SUPERSMART::Service::TreeService;

use Bio::SUPERSMART::Service::ParallelService;

use Bio::SUPERSMART::App::SubCommand;
use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);

# ABSTRACT: decomposes a backbone tree into individual clades

=head1 NAME

BBdecompose.pm - decomposes a backbone tree into individual clades

=head1 SYNOPSYS

smrt bbdecompose [-h ] [-v ] [-w <dir>] -b <file> -c <file> -a <file> -t <file> [-g ] [-o <file>]

=head1 DESCRIPTION

Given a rooted backbone phylogeny, a list of superclusters and a table of resolved taxa,
decomposes the backbone into its constituent, most recent, monophyletic clades, expands
these clades into all taxa from the provided table and assembles sets of alignments for
each clade that can be used for further tree inference within that clade.

Traverses the backbone tree to find the nearest monophyletic clade that groups
the exemplar leaves. In the default case, the clade is the genus that subtends the
two exemplars provided they are monophyletic in the backbone tree. If they are not, we
traverse upwards to find the nearest monophyletic set of genera.

Each clade is then expanded into its constituent set of species on the basis of the
taxon mapping file taxa. For those sets of species, the list of alignment values
is evaluated, and for each alignment whose average divergence does not exceed
CLADE_MAX_DISTANCE (given in the configuration file) but whose density of species
in the sets does exceed CLADE_MIN_DENSITY
the sequences for the focal species are written to a new file, in a directory that
groups them with the other relevant alignments for that clade.

=cut

sub options {
    my ($self, $opt, $args) = @_;
    my $tree_default    = "consensus.nex";
    my $taxa_default    = "species.tsv";
    my $aln_default     = "alignments";
    return (
        [
		     "backbone|b=s", 
		     "backbone tree as produced by e.g. 'smart bbinfer' and 'smrt consense', defaults to $tree_default", 
		 { arg => "file", default => $tree_default, galaxy_in => 1, galaxy_type => 'data'}
		],
        [
		     "alndir|a=s", 
		     "directory with alignments as produced by 'smrt aln'", 
		     { arg => "file", default => $aln_default, galaxy_in => 1, galaxy_type => 'data'}
		],
        [
		     "taxafile|t=s", 
		     "tsv (tab-seperated value) taxa file as produced by 'smrt taxize'", 
		     { arg => "file", default => $taxa_default, galaxy_in => 1, galaxy_type => 'data'}
		],
        [
		     "add_outgroups|g", 
		     "attempt to automatically add outgroup from sister genus for each clade, if sufficient marker overlap between clades", 
		     { galaxy_in => 1, galaxy_type => 'boolean' } 
		],
    );
}

sub validate {
    my ($self, $opt, $args) = @_;

    #  If alignment or taxa file is absent or empty, abort
    my @files = ( $opt->backbone, $opt->taxafile );
    foreach my $file ( @files ){
        $self->usage_error("need alignment and taxa files and backbone tree file") if not $file;
        $self->usage_error("file $file does not exist") unless (-e $file);
        $self->usage_error("file $file is empty") unless (-s $file);
    }
	$self->usage_error("need alndir") unless $opt->alndir;
}

sub run{
    my ($self, $opt, $args) = @_;

    # collect command-line arguments
	my $alndir = $self->process_inputdir( $opt->alndir );
    my $taxafile     = $opt->taxafile;
    my $backbone     = $opt->backbone;
    my $add_outgroup = $opt->add_outgroups;
    my $workdir      = $self->workdir;

    # instantiate helper objects
    my $mt     = Bio::SUPERSMART::Domain::MarkersAndTaxa->new;
    my $mts    = Bio::SUPERSMART::Service::MarkersAndTaxaSelector->new;
    my $ts     = Bio::SUPERSMART::Service::TreeService->new;
    my $config = Bio::SUPERSMART::Config->new;
    my $logger = $self->logger;

    # parse backbone tree
    $logger->info("Going to read backbone tree $backbone");
    (my $tree) = $ts->read_tree( '-file'   => $backbone );	
    $ts->remap_to_ti($tree);
	
    # parse taxon mapping
    $logger->info("Going to read taxa mapping $taxafile");
    my @taxa = $mt->parse_taxa_file( $taxafile );

    # now read the directory of alignments
    my @alignments = $mt->parse_aln_dir( $alndir );
  
	# decompose tree into clades and get the sets of species
    # extract_clades() returns a list of array references. the map operation results
    # therefore in a list of hash references with a single key, whose value is an array
    # reference, which should contain taxon IDs. XXX however, this is not the case, the
    # value is not an array reference with scalar taxon IDs, there is a deeper nesting
    # of arrays inside.
    my @clades = map { { 'ingroup' => $_ } } $ts->extract_clades($tree, @taxa);
	
    # get the exemplars
    for my $c ( @clades ) {
        my @ex = grep { defined $_ } map { $tree->get_by_name($_) } @{ $c->{'ingroup'} };
        $c->{'exemplars'} = [ keys %{{ map { $_->get_name => 1 } @ex }} ];
	}

    # get one outgroup species for each clade and append to species sets
	@clades = $self->_add_outgroups( \@clades, \@taxa, \@alignments ) if $add_outgroup;

	# collect alignments for clades and write them into respective clade directories
	my @species_for_clades = map {

		my $i = $_;
		$logger->info("Processing clade $i");
		
		# for each clade, collect a set of alignments
		my $clade = $clades[$i];
		my %ingroup  = map { $_ => 1 } @{ $clade->{'ingroup'} };
		my %outgroup = map { $_ => 1 } @{ $clade->{'outgroup'} };
		my %exemplars = map { $_ => 1 } @{ $clade->{'exemplars'} };
		
		my $mindens = $config->CLADE_MIN_DENSITY;
		my $maxdist = $config->CLADE_MAX_DISTANCE;

	    # check all alignments and assess whether they are suitable for the clade
	    my @clade_alignments = $mts->filter_clade_alignments( '-ingroup'  => $clade->{'ingroup'},
															  '-outgroup' => $clade->{'outgroup'},
															  '-alnfiles' => \@alignments,
															  '-clade'    => $i);

		# for the given set of taxa and alignments, get all subsets of taxa that share at least one marker
		my @all_taxa = map {@$_} map {values(%$_)} @clades;
		my %adj = $mts->generate_marker_adjacency_matrix(\@clade_alignments, [keys %outgroup, keys %ingroup]);
		my @subsets = @{$mts->get_connected_taxa_subsets(\%adj)};

		# sort subsets by size, decreasing
		@subsets = sort { scalar(@$b) <=> scalar(@$a)  } @subsets;

		# select the largest subset with exemplars present
		my @set;
		for my $s (@subsets) {
			if ( scalar grep( $exemplars{$_}, @$s)  >= 2) {
				@set = @$s;
				last;
			}
		}
		# skip if there are not enough species for clade tree
		my @names;
		if ( scalar (@set) < 3 ) {
			$logger->warn("Could not find sufficient data for species in clade $i. Skipping clade with taxa " . join(',', keys(%ingroup)));
		}
		else {
			_write_clade_alignments( $i, \@clade_alignments, \@set, $self->workdir );
			# write outgroup to file (skipped if already exists)
			_write_outgroup($i,[keys %outgroup],$workdir) if $add_outgroup;
			# write taxa file to clade directory
			@names = map { $mts->find_node($_)->taxon_name } @set;
			my @taxa_table = $mts->make_taxa_table( '-taxon_names' => \@names );
			my $taxafile = "clade$i/species.tsv";
			$mts->write_taxa_file( '-file' => $taxafile, '-table' => \@taxa_table );
		}
		{"clade$i"=> \@names}

	}  ( 0..$#clades);
	
	for ( @species_for_clades )  {
		my ($clade) = keys %$_;
		my @names = @{$_->{$clade}};
		$logger->info("Species in $clade: " . join(', ', @names));
	}

    $logger->info("DONE, results written into working directory $workdir");

    return 1;
}

# add outgroups to species sets representing clades
sub _add_outgroups {
	my ($self, $clades_ref, $taxa_ref, $alignments_ref) = @_;

	my @clades = @{$clades_ref};
	my @taxa = @{$taxa_ref};
	my @alignments = @{$alignments_ref};

	my $mt     = Bio::SUPERSMART::Domain::MarkersAndTaxa->new;
    my $mts    = Bio::SUPERSMART::Service::MarkersAndTaxaSelector->new;
    my $ts     = Bio::SUPERSMART::Service::TreeService->new;
	my $logger = $self->logger;

	# make a classification tree, source of candidate outgroups
	my $classtree = $ts->make_ncbi_tree( @taxa );
	
	# iterate over hashes, with key 'ingroup', value is an
	# array ref of taxon IDs
	my $counter = 0;
	for my $clade ( @clades ){
		my $ingroup = $clade->{'ingroup'};
		
		# TODO: Should extra_depth for more distant outgroup species be given as argument?
		my $extra_depth = 0;
		my @og = $mts->get_outgroup_taxa( $classtree, $ingroup, $extra_depth );
		
		# get the two species which occur in the most number of alignments
		my %aln_for_sp;
		for my $aln ( @alignments ) {
			my %fasta = $mt->parse_fasta_file($aln);
                for my $ogsp ( @og ) {
                    my @sp = grep { /taxon\|$ogsp[^\d]/ } keys %fasta;
                    @sp = map { $1 if $_ =~ m/taxon\|([0-9]+)/ } @sp;
                    $aln_for_sp{$_}++ for @sp;
                }
		}
		my @sorted_sp = sort { $aln_for_sp{$a} <=> $aln_for_sp{$b} } keys %aln_for_sp;
		my @outgroup = scalar @sorted_sp > 4 ? @sorted_sp[0..3] : @sorted_sp;
		
		$clade->{'outgroup'} = \@outgroup;
		$logger->info("Adding outgroup species " . join (', ', @outgroup) . " to clade #" . $counter);
		$counter++;
	}
	return (@clades);
}

# writes alignments for a given set of taxa in a subfolder within the clade
sub _write_clade_alignments {
	my ( $clade, $alns, $taxa, $workdir ) = @_;

	# make directory to write clade data to
	my $cladedir = File::Spec->catdir($workdir, "clade$clade");
	mkdir $cladedir if not -d $cladedir;

	# write all alignemnts to clade directory 
	my @alignments = @$alns;
	for my $i (0..$#alignments) {
		my %aln = %{$alignments[$i]};
		
		# make alignment directory 
		my $alndir = File::Spec->catdir($cladedir, 'alignments');
		mkdir $alndir if not -d $alndir;
		
		my $def = (keys %aln)[0];
		my $seed_gi = $1 if $def =~ /seed_gi\|([0-9]+)/;

		# file name e.g. 12345-clade0.fa
		my $filename = File::Spec->catfile($alndir, "$seed_gi-clade$clade.fa");
		open my $fh, '>', $filename or die $!;

		for my $defline( keys %aln ) {
			my $species = $1 if $defline =~ /taxon\|([0-9]+)/;

			# only write sequence if taxon is in the given list of taxa
			if ( grep($species, @$taxa) ) {
				print $fh '>', $defline, "\n", $aln{$defline}, "\n";
			}
		}
		close $fh;
	}
}

# writes outgroup species into a file in the clade folder
sub _write_outgroup {
    my ( $i, $og, $workdir ) = @_;
    my $filename = "$workdir/clade$i/clade$i-outgroup.txt";
    if  ( not -e $filename ) {
        open my $fh, '>', $filename or die $!;
        foreach my $sp ( @{$og} ) {
            print $fh  "$sp\n";
        }
        close $fh;
    }
}

1;
