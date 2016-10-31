package Bio::SUPERSMART::Service::InferenceService;
use strict;
use warnings;
use Bio::SUPERSMART::Service;
use Bio::Phylo::Util::Exceptions 'throw';
use Bio::Phylo::Util::CONSTANT '/looks_like/';
use base 'Bio::SUPERSMART::Service';

=head1 NAME

Bio::SUPERSMART::Service::InferenceService - Infers phylogenetic trees

=head1 DESCRIPTION

Provides functionality inferring trees and sets of trees by wrapping external tools.

=over

=item new

Object constructor. Returns the appropriate subclass. Arguments:

	'tool'    => tool name, e.g. examl, raxml, exabayes
	'workdir' => directory name for temp files
	'outfile' => (optional) name for output file
	'boostrap'=> (optional) number of boostrap replicates

=cut

sub new {
	my ( $package, %args ) = @_;
	my $self;

	# invoked as this class->new, construct subclass name and re-invoke
	if ( $args{'tool'} ) {
        	my $subclass = __PACKAGE__ . '::' . lc($args{'tool'});
        	delete $args{'tool'};
		if ( looks_like_class $subclass ) {
			$self = $subclass->new(%args);
		}
	}
	
	# invoked second time, as subclass->new, traverse up
	else {
		$self = $package->SUPER::new(%args);
	}

	# apply other arguments
	for my $key ( keys %args ) {
		$self->$key($args{$key});
	}
	$self->wrapper( $self->create );
	return $self;
}

=item create

Instantiates and configures the wrapper object, returns the instantiated wrapper.
This is an abstract method that needs to be implemented by the child class.

=cut

sub create {
    throw 'NotImplemented' => "missing create method in child class " . ref(shift);
}

=item configure

Provides the $config object to the wrapper so that settings defined in $config object can 
be applied to the wrapper. This is an abstract method that needs to be implemented by the 
child class.

=cut

sub configure {
    throw 'NotImplemented' => "missing confgure method in child class " . ref(shift);
}

=item run

Runs the analysis on the provided supermatrix. can potentially be run multiple
times (in parallel?) to do bootstrapping. returns a tree file. This is an abstract method 
that needs to be implemented by the child class.

=cut

sub run {
    throw 'NotImplemented' => "missing run method in child class " . ref(shift);
}

=item is_bayesian

Indicates whether the inference service is bayesian, in which case we will
not run a bootstrapping analysis. Note that certain wrappers, e.g. exabayes, 
need to override this

=cut

sub is_bayesian { 0 }

=item cleanup

Runs a cleanup operation for any intermediate files. This is an abstract method
that needs to be implemented by the child class.

=cut

sub cleanup {
    throw 'NotImplemented' => "missing cleanup method in child class ". ref(shift);
}

=item replicate

Getter/setter, is used internally to track the current bootstrap replicate

=cut

sub replicate {
	my ( $self, $rep ) = @_;
	if ( defined $rep ) {
		$self->{'replicate'} = $rep;
	}
	return $self->{'replicate'};
}

=item outfile

Getter/setter, is used to store the name of the output file.

=cut

sub outfile {
	my ( $self, $outfile ) = @_;
	if ( $outfile ) {
		$self->{'outfile'} = $outfile;
	}
	return $self->{'outfile'};
}

=item workdir

Getter/setter, is used to store the name of the working directory.

=cut

sub workdir {
	my ( $self, $workdir ) = @_;
	if ( $workdir ) {
		$self->{'workdir'} = $workdir;
	}
	return $self->{'workdir'} || '.';
}

=item workdir

Getter/setter, is used to store the number of bootstrap replicates

=cut

sub bootstrap {
	my ( $self, $bootstrap ) = @_;
   	if ( $bootstrap ) {
		$self->{'bootstrap'} = $bootstrap;
	}
	return $self->{'bootstrap'} || 1;
}


=item wrapper

Getter/setter, is used to store the inference tool wrapper (Bio::Tools::Run::Phylo::*)

=cut

sub wrapper {
	my ( $self, $wrapper ) = @_;
	if ( $wrapper ) {
		$self->{'wrapper'} = $wrapper;
	}
	return $self->{'wrapper'};
}

=item usertree

Getter/setter, is used to store the starting tree.

=cut

sub usertree {
	my ( $self, $tree ) = @_;
	if ( $tree ) {
		$self->{'usertree'} = $tree;
	}
	return $self->{'usertree'};
}

=item constraint_tree

Getter/setter, is used to store the constraint tree

=cut

sub constraint_tree {
	my ( $self, $tree ) = @_;
	if ( $tree ) {
		$self->{'constrainttree'} = $tree;
	}
	return $self->{'constrainttree'};
}

=back 

=cut


1;
