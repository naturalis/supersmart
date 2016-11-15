package Bio::SUPERSMART::Service::InferenceService::raxml;
use strict;
use warnings;
use Cwd 'abs_path';
use File::Spec;
use Bio::Tools::Run::Phylo::Raxml;
use Bio::SUPERSMART::Service::InferenceService;
use base 'Bio::SUPERSMART::Service::InferenceService';

=head1 NAME

Bio::SUPERSMART::Service::InferenceService::raxml - Infers phylogenetic trees

=head1 DESCRIPTION

Provides functionality inferring trees and sets of trees by wrapping RAxML.

=over

=item create

Instantiates and configures the wrapper object, returns the instantiated wrapper.

=cut

sub create {
    my $self = shift;
    my $tool = Bio::Tools::Run::Phylo::Raxml->new;

    # configure raxml runner
    $tool->w( abs_path($self->workdir) );

    return $tool;
}

=item configure

Provides the $config object to the wrapper so that settings defined in $config object can
be applied to the wrapper.

=cut

sub configure {
    my ( $self ) = @_;
    my $tool   = $self->wrapper;
    my $config = $self->config;

	# number of threads
	if ( $config->NODES > $self->bootstrap or $self->{'rapid_boot'} ) {
		$tool->T($config->NODES);
	}

	my ($volume,$directories,$id) = File::Spec->splitpath( $self->outfile );
    $tool->outfile_name($id);
    $tool->m($config->RAXML_MODEL);
    $tool->N($config->RAXML_RUNS);
    $tool->p($config->RANDOM_SEED);

	# set starting tree if given
	if ( my $tree = $self->usertree and  ! $self->{'rapid_boot'} and ! $self->constraint_tree ) {
		$self->logger->info("Setting starting tree $tree");
		$tool->t($tree);
	}
	if ( my $ct = $self->constraint_tree ) {
		$self->logger->info("Setting constraint tree $ct");
		$tool->g($ct);
	}

}

=item run

Runs the analysis on the provided supermatrix. Returns a tree file.

=cut

sub run {
    my ( $self, %args ) = @_;
    my $t      = $self->wrapper;
    my $logger = $self->logger;

    # create run ID
    my $id = $t->outfile_name;

    # create output file name
    my $treefile = File::Spec->catfile( $t->w, 'RAxML_bestTree.' . $id );

    # configure rapid bootstrap options
    if ( $self->{'rapid_boot'} ) {

		# Do ML search and boostrapping in one go
		$t->f('a');

        # bootstrap random seed
        $t->x($self->config->RANDOM_SEED);

		# we will return the best ML tree annotated with bootstrap values as node labels
        $treefile = File::Spec->catfile( $t->w, 'RAxML_bootstrap.' . $id );

		# set number of rapid boostrap replicates; override number of runs
		$t->N($self->bootstrap);

		# tell RaXML to put branch lenghts on bootstrap trees
		$t->k(1);
    }
    # run raxml, returns bioperl tree
    my $bptree = $t->run($args{'matrix'});
    $logger->fatal('RAxML inference failed; outfile not present') if not -e $treefile;

    return $treefile;
}

=item cleanup

Cleans up any intermediate files.

=cut

sub cleanup {
    my $self = shift;
    my $dir = $self->wrapper->w;
    my $outstem = $self->wrapper->outfile_name;
    unlink $dir . '/RAxML_info.' . $outstem;
    unlink $dir . '/RAxML_bipartitionsBranchLabels.' . $outstem;
    unlink $dir . '/RAxML_bipartitions' . $outstem;
	unlink $dir . '/RAxML_bestTree.' . $outstem;
	unlink $dir . '/RAxML_bootstrap.' . $outstem;

    opendir my $dh, $dir or die $!;
    while ( my $entry = readdir $dh ) {
        for my $prefix ( qw(parsimonyTree log result) ) {
            if ( $entry =~ /RAxML_${prefix}.${outstem}.RUN.\d+/ ) {
                unlink $dir . '/' . $entry;
            }
        }
    }
}

=item rapid_boot

=cut

sub rapid_boot {
	my $self = shift;
	$self->{'rapid_boot'} = shift if @_;
	return $self->{'rapid_boot'};
}


=back

=cut

1;
