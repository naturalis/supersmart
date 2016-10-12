package Bio::SUPERSMART::App::smrt::Command::Orthologize;

use strict;
use warnings;

use Bio::SUPERSMART::Config;
use Bio::SUPERSMART::Domain::MarkersAndTaxa;
use Bio::SUPERSMART::Service::SequenceGetter;

use Bio::SUPERSMART::App::SubCommand;
use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);

# ABSTRACT: creates orthologous clusters of aligned sequences

=head1 NAME

Align.pm - assesses orthology in different sequence alignments and merges them into orthologous clusters

=head1 SYNOPSYS


=head1 DESCRIPTION

Given a list of aligned candidate clusters, assigns orthology among the
clusters by performing reciprocal BLAST searches on the seed sequences
around which the clusters were assembled. Writes sequence clusters as FASTA files
into a user-specified dirsctory ('./clusters/' by default). Also writes a file 'cluster-names.tsv'
containing the definition lines of the first sequence of each cluster.

=cut

sub options {
	my ($self, $opt, $args) = @_;
	my $indir_default = 'alignments';
	my $outdir_default = 'clusters';
	return (
		[
		 "indir|i=s",
		 "directory or zip file with sequence alignments, as produced by 'smrt align'",
		 { arg => "file", default => $indir_default, galaxy_in => 1, galaxy_format => 'tabular', galaxy_type => "data", galaxy_label => 'alignments' }
		],
		[
		 "outdir|o=s",
		 "write alignments to specified directory name, defaults to $outdir_default",
		 { arg => "dir", default => $outdir_default }
		],
		[
		 "zip|z",
		 "zip output directory containing orthologous sequence alignments",
		],
	);
}

sub validate {
	my ($self, $opt, $args) = @_;

	my $in = $opt->indir;
	$self->usage_error("no indir argument given") if not $in;
}

sub run {
	my ( $self, $opt, $args ) = @_;

	# collect command-line arguments
	my $indir  = $self->process_inputdir( $opt->indir );
	my $outdir = $opt->outdir;
	my $workdir = $self->workdir;

    # create directory for merged alignments
    my $dir = $self->make_outputdir($opt->outdir);  
	
	# instantiate helper objects
	my $sg = Bio::SUPERSMART::Service::SequenceGetter->new;
	my $mt = Bio::SUPERSMART::Domain::MarkersAndTaxa->new;
	my $config  = Bio::SUPERSMART::Config->new;
	my $log     = $self->logger;
  
	# get seed GIs of input alignment files
	my @alnfiles = $mt->parse_aln_dir( $indir );
	my @gis = $mt->extract_seed_gis( @alnfiles );

	# merge alignments to orthologous clusters and store in output directory
	my @clfiles = $sg->merge_alignments( $config->BACKBONE_MAX_DISTANCE, $workdir, $indir, $outdir,  @gis );
	
	# write a file with the cluster names (definitions) for each cluster
	my $names_file = "cluster-names.tsv";
	$sg->write_cluster_definitions( "$outdir/$names_file", @clfiles );

	# cleanup working directory and zip output if required
	$self->cleanup_inputdir( $opt->indir );
	$outdir = $self->zip_outputdir( $outdir ) if $opt->zip;

	$log->info("DONE, results written to $outdir");
	return 1;
}

1;
