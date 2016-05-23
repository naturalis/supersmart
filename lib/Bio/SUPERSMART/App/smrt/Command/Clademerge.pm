package Bio::SUPERSMART::App::smrt::Command::Clademerge;

use strict;
use warnings;

use List::MoreUtils 'uniq';
use File::Spec;

use Bio::SUPERSMART::Config;
use Bio::Phylo::Factory;
use Bio::Phylo::IO 'parse_matrix';

use Bio::SUPERSMART::Domain::MarkersAndTaxa;
use Bio::SUPERSMART::Service::TreeService;
use Bio::SUPERSMART::Service::ParallelService;
use Bio::SUPERSMART::Service::MarkersAndTaxaSelector;

use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);

# ABSTRACT: merges sets of alignments into input files for clade inference

=head1 NAME

Clademerge.pm - For each decomposed clade, merges the set of alignments assembled for this clade into an input file for tree inference.

=head1 SYNOPSYS

smrt clademerge [-h ] [-v ] [-w <dir>] [-o <format>]

=head1 DESCRIPTION

Given a working directory, traverses it looking for subdirectories of the pattern
clade*. Perusing each of these, it merges the *.fa (FASTA) files in it and
produces a single output file that can be analysed by the subcommand cladeinfer.

=cut


sub options {
    my ($self, $opt, $args) = @_;       
    return (
        [     
			  "outformat|o=s", 
			  "output format for merged clade files (phylip or nexml), defaults to 'nexml'", 
			  { arg=>"format", default=> 'nexml'} 
		],
        [ 
		      "enrich|e", 
		      "enrich the selected markers with additional haplotypes", 
		      { galaxy_in => 1, galaxy_type => 'boolean'} 
		],
    );  
}

sub validate {};

sub run {
    my ($self, $opt, $args) = @_;       
    
    my $workdir   = $self->workdir;
        
    # instantiate helper objects
    my $mt      = Bio::SUPERSMART::Domain::MarkersAndTaxa->new;
	my $config  = Bio::SUPERSMART::Config->new;
    my $log     = $self->logger;


    # collect candidate dirs
    $log->info("Going to look for clade data in $workdir");
    my @dirs;
    opendir my $odh, $workdir or die $!;
    while( my $dir = readdir $odh ) {
        if ( $dir =~ /^clade\d+$/ and -d "${workdir}/${dir}" ) {
            push @dirs, $dir;
        }
    }
	# merge alignments for each clade, don't include this in the pmap below since
	# orthologize_cladedir also uses pmap, hence avoid nested pmap calls
	for my $cladedir ( @dirs ) {
		my $curr_dir = File::Spec->catdir($workdir, $cladedir);
		my $clusterdir = File::Spec->catdir($cladedir, 'clusters');
		$self->make_outputdir( $clusterdir );
		$log->debug("Will write merged alignment clusters to $clusterdir");
		$mt->orthologize_cladedir(
			'dir'=> $curr_dir, 
			'outdir'=>$clusterdir, 
			'maxdist'=>$config->CLADE_MAX_DISTANCE );
	}
	
    # enrich (if requested in argument) and write matrix to file
    my @result = grep { defined $_ and -e $_ } pmap {
		(my $clade) = @_;
		my $dir = $clade;
		
		my $cladedir = File::Spec->catdir($self->workdir, $clade);
		my $clusterdir = File::Spec->catdir( $cladedir, 'clusters' );
		
		# check if we have merged alignments
		if ( scalar( glob ( "${clusterdir}/*.fa") ) ) {
			$mt = Bio::SUPERSMART::Domain::MarkersAndTaxa->new($clusterdir, $config->CLADE_MIN_COVERAGE);
			$mt->write_clade_matrix( 
				'markersfile' => "${cladedir}/${clade}-markers.tsv",
				'outfile' => $opt->outformat eq 'phylip' ? "${cladedir}/${clade}.phy" : "${cladedir}/${clade}.xml",
				'max_markers' => $config->CLADE_MAX_COVERAGE,
				'enrich' => $opt->enrich,
				'format' => $opt->outformat);
		} 
		else {
			$log->warn("Clade $clade is missing merged alignments, possibly due to low coverage. Skipping.");
		}
	} @dirs;
    
    # report result
    $log->info("Wrote outfile $_") for @result;
    $log->info("DONE");
    return 1;
}

1;
