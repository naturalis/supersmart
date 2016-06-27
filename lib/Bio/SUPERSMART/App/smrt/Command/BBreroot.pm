package Bio::SUPERSMART::App::smrt::Command::BBreroot;

use strict;
use warnings;
use File::Temp 'tempfile';

use Bio::Phylo::Factory;
use Bio::SUPERSMART::Service::TreeService;
use Bio::SUPERSMART::Service::MarkersAndTaxaSelector;
use Bio::Phylo::IO qw(parse parse_tree);

use Bio::SUPERSMART::Service::ParallelService;

use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);


# ABSTRACT: reroots backbone tree

=head1 NAME

BBinfer.pm - reroots backbone tree

=head1 SYNOPSIS

 smrt bbreroot [-h ] [-v ] [-w <dir>] [-l <filename>] -t <file> -b <file> \
	[-g <names>] [-o <filename>] 

=head1 DESCRIPTION

Reroots input backbone tree(s) in newick format (--backbone or -b argument)
ideally with one or more outgroups, which can be specified as comma-separated
values for the --outgroup (or -g) argument.

If no outgroups are specified, a taxa table (as produced by smrt taxize) will
be used to attempt to reconcile the topology with the taxonomy by selecting the
rooting that minimizes the amount of taxonomic non-monophyly on the input tree.
Note that there often are multiple solutions possible and only the first one
will be (arbitrarily) selected.

A taxafile will be used regardless. By default this is the 'species.tsv' file
in the current working directory, but this can be altered using the --taxafile
or -t argument.

=cut

sub options {
	my ($self, $opt, $args) = @_;		
	my $outfile_default = "backbone-rerooted.dnd";
	my $taxa_default = "species.tsv";
	my $tree_default = "backbone.dnd";
	return (
		[
    		 "taxafile|t=s", 
    		 "tsv (tab-seperated value) taxa file as produced by 'smrt taxize'", 
    		 { arg => "file", default => $taxa_default, galaxy_in => 1, galaxy_type => 'data'}
		],
		[
    		 "backbone|b=s", 
    		 "a backbone tree(s) file as produced by 'smrt bbinfer'", 
    		 { arg => "file", default => $tree_default, galaxy_in => 1, galaxy_type => 'data'}
		],	
		[
    		 "outgroup|g=s", 
    		 "one or multiple taxa (names or NCBI identifiers, separated by commata) representing the outgroup at "
    		 ."which the tree is rerooted. Outgroup must be enclosed in quotes.", 
    		 { arg => "taxon,taxon,...", galaxy_in => 1, galaxy_type => 'text' } 
		],
		[
    		 "outgroup_tree|p=s", 
    		 "tree in newick format to extract outgroup from. Outgroup taxa are the terminals of the smallest subtree below the root  ",
    		 { arg => "file" }
		],
		[
		     "outfile|o=s", 
		     "name of the output tree file (in newick format), defaults to '$outfile_default'", 
		     { default => $outfile_default, arg => "file", galaxy_out => 1, galaxy_type => 'data', galaxy_label => 'rerooted backbone'}
		],
		[
		     "smooth|s", 
		     "smooth tip heights left and right of root (i.e. midpointify)",
		     { galaxy_in => 1, galaxy_type => "boolean"}
		],
		[    "ultrametricize|u", 
    		 "adjust terminal branch lengths to yield an ultrametric tree",
			 { galaxy_in => 1, galaxy_type => "boolean"}
		],
	);	
}

sub validate {
	my ($self, $opt, $args) = @_;		

	#  If taxa file or backbone tree is absent or empty, abort  
	my @files = ( $opt->taxafile, $opt->backbone );
	foreach my $file ( @files ){
		$self->usage_error("need backbone and taxafile") if not $file;
		$self->usage_error("file $file does not exist") unless (-e $file);
		$self->usage_error("file $file is empty") unless (-s $file);			
	}
}

sub run {
	my ($self, $opt, $args) = @_;		
	
	# collect command-line arguments
	my $backbone = $opt->backbone;
	my $taxafile = $opt->taxafile;
	my $smooth   = $opt->smooth;
	my $outfile  = $self->outfile;
		
	# instantiate helper objects
	my $ts  = Bio::SUPERSMART::Service::TreeService->new;
	my $mt  = Bio::SUPERSMART::Domain::MarkersAndTaxa->new;
	my $mts = Bio::SUPERSMART::Service::MarkersAndTaxaSelector->new;
	my $log = $self->logger;
	
	# identify outgroup once
	my $outgroup;
	if ( my $csv = $opt->outgroup ) {
		$outgroup = [map { $_->ti } $mts->get_nodes_for_names(split /,/,$csv)];
	}	
	elsif ( my $treefile = $opt->outgroup_tree ) {
		$outgroup = $self->_get_smallest_outgroup( $treefile );
	}
	
	# prepare taxa data
	my @records = $mt->parse_taxa_file($taxafile);
		
	# iterate over trees
	open my $in, '<', $backbone or die $!;
	chomp(my @backbone_trees = <$in>);
	close $in;
	
	# mapping tables for faster id and taxon name lookup
	my %ti_to_name;
	my %name_to_ti;

	# reroot input trees in parallel
	my @rerooted_trees = pmap{	
		my $treestr = $_;

		# read tree
		$log->debug("Reading tree from string");
		my $tree = parse_tree( '-string' => $treestr, 
							   '-format' => 'newick' );	
		
		# create id mapping table
		if ( ! scalar(%ti_to_name) ) {
			%ti_to_name = $ts->make_mapping_table($tree);
			%name_to_ti = reverse(%ti_to_name);
		}
		
		# map identifiers
		$log->debug("Mapping taxon names to taxon IDs");
		$tree = $ts->remap($tree, %name_to_ti);
		$tree->resolve;

		# Perform rerooting at outgroup, if given		
		if ( $outgroup ) {			
			my @ranks = ('forma', 'varietas', 'subspecies', 'species');	
			$log->debug("Rerooting at outgroup $outgroup");
			$ts->outgroup_root(
				'-tree'     => $tree,
				'-ids'      => $outgroup,
				'-ranks'    => \@ranks,
				'-records'  => \@records,
				);		
		} 
		
		# Try to minimize paraphyly if no outgroup given
		else {	
			my $level = $mt->get_highest_informative_level(@records);	
			$log->debug("highest informative taxon level : $level");		
			$tree = $ts->reroot_tree($tree, \@records, [$level]);
		}
		
		# smooth the basal branch, if requested
		if ( $smooth ) {
			$log->debug("smoothing out diff between left and right tip heights");
			$ts->smooth_basal_split($tree);
		}
		
		# clean up labels and map to taxon names
		$tree = $ts->remap($tree, %ti_to_name);
		$ts->remove_internal_names($tree);
        $log->info("Rerooted backbone tree");
		
		$tree->ultrametricize if $opt->ultrametricize;

		return( $tree->to_newick )

	} @backbone_trees;
	
	$log->warn("Number of rerooted trees different than number of input trees") if scalar(@backbone_trees) != scalar(@rerooted_trees);
	
	# write output file
	open my $out, '>', $outfile or die $!;			
	print $out $_  . "\n" for @rerooted_trees;
	close $out;

	$log->info("DONE, results written to $outfile");		
}


# given a tree, returns the terminal names of the smallest of the two subtrees below the
# tree root.
sub _get_smallest_outgroup {
	my ($self, $treefile) = @_;
	my $logger = $self->logger;
	
	# parse tree
	my $tree = parse_tree( '-file' => $treefile, 
						   '-format' => 'newick' );	
	$tree->resolve;

	# get subtrees below root
	my @children = @{$tree->get_root->get_children};
	$logger->warn("Expecting two children of root node") if not  scalar(@children) == 2;
	
	# collect sets of terminal names for subtrees
	my @subtree_terminals = map { my $c = $_; [ map {$_->get_name} @{ $c->get_terminals }] } @children;

	# sort sets of terminal names by size, increasing
	@subtree_terminals = sort { scalar(@$a) <=> scalar(@$b)  } @subtree_terminals;
	
	my $outgroup = $subtree_terminals[0];
	$logger->info("Determined outgroup '" . join( ',', @$outgroup) . "' from tree $treefile");
	
	return $outgroup;

}

1;
