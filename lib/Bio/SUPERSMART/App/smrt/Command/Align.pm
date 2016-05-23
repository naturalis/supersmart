package Bio::SUPERSMART::App::smrt::Command::Align;

use strict;
use warnings;

use Bio::Phylo::Matrices::Matrix;
use Bio::SUPERSMART::Domain::MarkersAndTaxa;
use Bio::SUPERSMART::Service::SequenceGetter;
use Bio::SUPERSMART::Service::MarkersAndTaxaSelector;
use Bio::SUPERSMART::Service::ParallelService;

use Bio::PrimarySeq;
use Bio::AlignIO;
use File::Spec;

use Bio::SUPERSMART::App::SubCommand;
use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);


# ABSTRACT: writes multiple sequence alignments for a list of taxa

=head1 NAME

Align.pm - writes phylogenetically informative clusters for taxa

=head1 SYNOPSYS

smrt align [-h ] [-v ] [-w <dir>] -i <file> [-o <file>] 

=head1 DESCRIPTION

Given an input table of resolved taxa, performs multiple sequence alignment 
of all potentially useful PhyLoTA clusters for the taxa of interest. 
Alignments will be stored in an output directory, given via the '--dirname' argument, 
with default 'alignments/' in the current workind directory

=cut

sub options {
    
    my ($self, $opt, $args) = @_;       
    my $taxa_default = 'species.tsv';
	my $outdir_default = 'alignments';
	
    return (
		[
		 "infile|i=s", 
		 "taxa file (tab-seperated value format) as produced by 'smrt taxize'", 
		 { arg => "file", default => $taxa_default, galaxy_in => 1, galaxy_format => 'tabular', galaxy_type => "data", galaxy_label => 'taxa file' }
		],
		[
		 "outdir|o=s", 
		 "write alignments to specified directory name, defaults to $outdir_default", 
		 { arg => "dir", default => $outdir_default }
		],
		[
		 "zip|z", 
		 "zip output alignment directory",
		]
    );  
}

sub validate {
    my ($self, $opt, $args) = @_;       
    
    #  If the infile is absent or empty, abort  
    my $file = $opt->infile;
    $self->usage_error("no infile argument given") if not $file;
    $self->usage_error("file $file does not exist") unless (-e $file);
    $self->usage_error("file $file is empty") unless (-s $file);    
}

sub run {
    my ($self, $opt, $args) = @_;
    
    # collect command-line arguments
    my $infile  = $opt->infile;
    
    # create directory for alignments if specified in argument
    my $dir = $self->make_outputdir($opt->outdir);  
    
    # instantiate helper objects
    my $log = $self->logger;
    my $mts = Bio::SUPERSMART::Service::MarkersAndTaxaSelector->new;
    my $mt  = Bio::SUPERSMART::Domain::MarkersAndTaxa->new;
    my $sg  = Bio::SUPERSMART::Service::SequenceGetter->new;
    my $ps  = 'Bio::SUPERSMART::Service::ParallelService';
    
    # instantiate taxonomy nodes from infile
    my @nodes = $mts->get_nodes_for_table($mt->parse_taxa_file($infile));     
    $log->info("Found ".scalar(@nodes)." nodes in taxa file $infile");
    
    # fetch all clusters for the focal nodes and organize them into
    # roughly even-sized chunks (they were initially sorted from big
    # to small so the first worker would get all the big ones).
    my @sorted_clusters = $mts->get_clusters_for_nodes(@nodes); 
    my @clusters = $ps->distribute(@sorted_clusters);   
    
    # NOTE: We do the below in sequential mode, which may take some time. This
    # is for the following reason: One sequence (GI) can easily be part of
    # multiple clusters. Since including one sequence into multiple alignments
    # (which would then most likely be merged later anyway ) uses space and
    # resources, we do not include a GI  more than once. To make sure that we
    # skip the GI always in the same alignment, we avoid race conditions in
    # parallel processing. Also see issue #56
    # Collect and filter sequences for all clusters
    my %ti =  map { $_->ti => 1 } @nodes;
    my ( %seqs, %seen, @clinfos );
    $log->info("Going to collect sequences for " . scalar(@clusters) . " clusters");
    for my $cl ( @clusters ) {
    
        # get cluster information
        my $ci      = $cl->{'ci'};
        my $type    = $cl->{'cl_type'};
        my $single  = $sg->single_cluster($cl);
        my $seed_gi = $single->seed_gi;
        my $mrca    = $single->ti_root->ti;
        my $clinfo  = $single->clinfo;
        
        # XXX we should be able to parameterize this more so that we can filter
        # a set down to $max haplotypes per taxon
        my @seqs = $sg->filter_seq_set($sg->get_sequences_for_cluster_object($cl));   
        $log->debug("fetched ".scalar(@seqs)." sequences for cluster $clinfo");

        # filter out sequences that we have processed before or that are
        # of uninteresting taxa
        my @matching = sort {$a->gi <=> $b->gi} grep { $ti{$_->ti} } @seqs;
        @matching    = grep { ! $seen{$_->gi} } @matching;
        
        # let's not keep the ones we can't build alignments from
        if ( scalar @matching > 1 ) {
            $log->info("Will align cluster: $clinfo");
            $seqs{$clinfo} = [
                map {
                    $seen{$_->gi} = 1;          # skip this seq next time                                       
                    $_->to_primary_seq(         # Bio::PrimarySeq
                        'mrca'    => $mrca,
                        'seed_gi' => $seed_gi,
                    );
                } @matching
            ];
            push @clinfos, $clinfo; # to keep optimized distribution for pmap
        }       
    }

    # now make the alignments in parallel mode
    my @result = pmap {             
        my $clinfo = shift;
        
        # align to file
        my $filename = File::Spec->catfile( $dir, $clinfo.'.fa' );      
        $sg->align_to_file( $seqs{$clinfo} => $filename );         
		$filename;
     
    } @clinfos;
    
	my $output = $opt->zip ? $self->zip_outputdir($dir) : $dir;

    $log->info("DONE, results written to $output.");    
    
    return 1;
}
    
1;
