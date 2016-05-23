package Bio::Tools::Run::Phylo::StarBEAST;
use strict;
use version;
use XML::Twig;
use Template;
use File::Temp 'tempfile';
use Bio::AlignIO;
use Bio::Phylo::Factory;
use Bio::Phylo::IO 'parse';
use Bio::Phylo::Forest::Tree;
use Bio::Align::DNAStatistics;
use Bio::Tree::DistanceFactory;
use Bio::Tools::Run::Phylo::PhyloBase;
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::Util::CONSTANT ':objecttypes';
use Bio::Phylo::Util::Exceptions 'throw';
use base qw(Bio::Tools::Run::Phylo::PhyloBase);

our $PROGRAM_NAME = 'beast';
our @beast_PARAMS = (
	'mc3_chains',       # number of chains
	'mc3_delta',        # temperature increment parameter
	'mc3_temperatures', # a comma-separated list of the hot chain temperatures
	'mc3_swap',         # frequency at which chains temperatures will be swapped
	'seed',             # Specify a random number generator seed
	'threshold',        # Full evaluation test threshold (default 1E-6)
	'beagle_order',     # set order of resource use
	'beagle_scaling',   # specify scaling scheme to use
	'beagle_rescale',   # frequency of rescaling (dynamic scaling only)	
	'beagle_instances', # divide site patterns amongst instances	
	'threads',          # set/limit number of threads for likelihood calculation 
);
our @beast_SWITCHES = (
	'verbose',          # Give verbose XML parsing messages
	'warnings',         # Show warning messages about BEAST XML file
	'strict',           # Fail on non-conforming BEAST XML file
	'beagle_CPU',       # use CPU instance		
	'beagle_GPU',       # use GPU instance if available
	'beagle_SSE',       # use SSE extensions if available
	'beagle_cuda',      # use CUDA parallization if available
	'beagle_opencl',    # use OpenCL parallization if available
	'beagle_single',    # use single precision if available
	'beagle_double',    # use double precision if available	
	'overwrite',        # Allow overwriting of log files
);
my $beast_ns  = 'http://beast.bio.ed.ac.uk/BEAST_XML_Reference#';
my $beast_pre = 'beast';
my $fac = Bio::Phylo::Factory->new;
my $log = Bio::Phylo::Util::Logger->new;

=head1 NAME

Bio::Tools::Run::Phylo::StarBEAST - module for interface with the BEAST software 

=head1 SYNOPSYS

 use Bio::Tools::Run::Phylo::StarBEAST;
 use FindBin '$Bin';
 
 my $beast = Bio::Tools::Run::Phylo::StarBEAST->new;
 $beast->chain_length(100);
 $beast->sample_freq(100);
 $beast->run( "$Bin/beast.xml" );


=head1 DESCRIPTION

This module provides an interface to the  *BEAST software
(http://beast.bio.ed.ac.uk) for inference of phylogenetic trees 
using Markov-Chain Monte-Carlo methods.   

=head1 METHODS

In addition to the methods described below, all named arguments described for the
constructor are available as methods as well. For example, the argument C<-verbose>
becomes:

 $beast->verbose(1);

=over

=item new

The constructor takes the following optional named arguments that require a numerical 
value:

	'-mc3_chains',       # number of chains
	'-mc3_delta',        # temperature increment parameter
	'-mc3_temperatures', # a comma-separated list of the hot chain temperatures
	'-mc3_swap',         # frequency at which chains temperatures will be swapped
	'-seed',             # Specify a random number generator seed
	'-threshold',        # Full evaluation test threshold (default 1E-6)
	'-beagle_order',     # set order of resource use
	'-beagle_scaling',   # specify scaling scheme to use
	'-beagle_rescale',   # frequency of rescaling (dynamic scaling only)	
	'-beagle_instances', # divide site patterns amongst instances	
    '-threads',          # threads for likelihood calculation

In addition, the following named arguments that require a boolean are available:

	'-verbose',          # Give verbose XML parsing messages
	'-warnings',         # Show warning messages about BEAST XML file
	'-strict',           # Fail on non-conforming BEAST XML file
	'-beagle_CPU',       # use CPU instance		
	'-beagle_GPU',       # use GPU instance if available
	'-beagle_SSE',       # use SSE extensions if available
	'-beagle_cuda',      # use CUDA parallization if available
	'-beagle_opencl',    # use OpenCL parallization if available
	'-beagle_single',    # use single precision if available
	'-beagle_double',    # use double precision if available	
	'-overwrite',        # Allow overwriting of log files


=cut

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);
    $self->_set_from_args(
        \@args,
        '-methods' => [ @beast_PARAMS, @beast_SWITCHES ],
        '-create'  => 1,
    );
    my ($out) = $self->SUPER::_rearrange( [qw(OUTFILE_NAME)], @args );
    $self->outfile_name( $out || '' );
    return $self;
}

sub _setparams {
    my ( $self, $file ) = @_;
    my $param_string = '';

	# iterate over parameters and switches
    for my $attr (@beast_PARAMS) {
        my $value = $self->$attr();
        next unless defined $value;
        $param_string .= ' -' . $attr . ' ' . $value;
    }
    for my $attr (@beast_SWITCHES) {
        my $value = $self->$attr();
        next unless $value;
        $param_string .= ' -' . $attr;
    }

    # Set default output file if no explicit output file has been given
    if ( ! $self->outfile_name ) {
        my ( $tfh, $outfile ) = $self->io->tempfile( '-dir' => $self->work_dir );
        close $tfh;
        undef $tfh;
        $self->outfile_name($outfile);
    }
    
    # set input file name
	$param_string .= ' ' . $file;
    
    # hide stderr
    my $null = File::Spec->devnull;
    $param_string .= " > $null 2> $null" if $self->quiet() || $self->verbose < 0;

    return $param_string;
}

=item logfile_name

Getter/setter for the location of a log file that can be loaded in tracer
to assess convergence and ESS

=cut

sub logfile_name {
	my $self = shift;
	$self->{'_logfile_name'} = shift if @_;
	return $self->{'_logfile_name'};
}

=item beastfile_name

Getter/setter for the location of the xml file that
is generated as input for *BEAST

=cut

sub beastfile_name {
	my $self = shift;
	$self->{'_beastfile_name'} = shift if @_;			
	
	# set filename to temporary one if not present
	if ( ! $self->{'_beastfile_name'} ) {
		my ($fh, $fn) = tempfile();
		$self->{'_beastfile_name'} = $fn; 
	}	
	return $self->{'_beastfile_name'};
}

=item log_freq

Getter/setter for the logging frequency.

=cut


sub log_freq {
	my $self = shift;
	$self->{'_log_freq'} = shift if @_;
	return $self->{'_log_freq'};

}

=item chain_length

Returns the number of steps for the Markov chain for running *BEAST.
If chain length is provided it is set to the given value.

=cut

sub chain_length {
	my $self = shift;
	$self->{'_chain_length'} = shift if @_;
	return $self->{'_chain_length'};
}

=item outgroup

Getter/Setter for outgroup

=cut

sub outgroup {
	my $self = shift;
	$self->{'_outgroup'} = shift if @_;
	return $self->{'_outgroup'};
}

=item sample_freq

Returns the number determining how often the Markov Chain is sampled.
If the sample frequence is provided it is set to the given value.

=cut

sub sample_freq {
	my $self = shift;
	$self->{'_sample_freq'} = shift if @_;
	return $self->{'_sample_freq'};
}

=item collapse_species

Getter/Setter for flag if taxa of rank subspecies, varietas and formas of one 
species should be collapsed in the inferred species tree

=cut

sub collapse_species {
	my $self = shift;
	$self->{'_collapse_species'} = shift if @_;
	return $self->{'_collapse_species'};
}

=item rebuild

Getter/setter for flag to specify whether to rebuild previously existing 
input files (which may have been edited by hand)

=cut

sub rebuild {
	my $self = shift;
	$self->{'_rebuild'} = shift if @_;
	return $self->{'_rebuild'};
}

=item template

Getter/setter for a L<Template> file to be used for generating BEAST
input

=cut

sub template {
	my $self = shift;
	$self->{'_template'} = shift if @_;
	return $self->{'_template'};
}

=item program_name

Returns the name of the executable.

=cut

sub program_name { $PROGRAM_NAME }

=item program_dir

Not yet implemented.

=cut

sub program_dir { undef }

=item run

Runs *Beast using the data in the provided NeXML file. 
Returns the output file generated from the *Beast run.

=cut

sub run {
	my ($self,$nexml) = @_;
	my $filename = $self->beastfile_name();
	
	# let's assume people may want to finetune their BEAST runs
	if ( -e $filename and not $self->rebuild ) {
		$log->warn("not overwriting previous input file $filename")
	}
	else {
	
		# parse the input file
		$log->info("going to read nexml file $nexml");	
		$self->_alignment( parse(
			'-format' => 'nexml',
			'-file'   => $nexml,
			'-as_project' => 1,
		) );

		# create BEAST xml
		if ( my $template = $self->template ) {
			# interpolate BEAST xml template
			my $tt = Template->new({ 'ABSOLUTE' => 1 });
			
			
			my $args = { 'data'         => $self->_alignment,
						 'chain_length' => $self->chain_length,
						 'sample_freq'  => $self->sample_freq,
						 'outfile_name' => $self->outfile_name,
						 'logfile_name' => $self->logfile_name,
			};
			
			# set outgroup if given to constrain monophyly of ingroup
			if ( $self->outgroup ) {
				my %og = map { $_=>1 } @{$self->outgroup};
				$args->{'outgroup'} = \%og;
			}

			$tt->process( $template, $args, $filename ) or throw 'API' => $tt->error();
		}
		else {
			# nope, need a template
			$log->fatal("no template provided!");
			throw 'FileError' => "*BEAST xml template file not provided";			
		}
	}
	
	# create the invocation and run the command	
	my $command = $self->executable . $self->_setparams($filename);
	$log->info("going to execute '$command'");	
    	my $status  = system $command;
    
    	# fetch the output file
    	my $outfile = $self->outfile_name();
    	if ( !-e $outfile || -z $outfile ) {
        	$self->warn("*BEAST call had status of $status: $? [command $command]\n");
        	return undef;
	}
	return $outfile;
}

=item version

Returns the version number of *Beast used.

=cut

sub version {
    my ($self) = @_;
    my $exe;
    return undef unless $exe = $self->executable;
    my $string = `$exe -version 2>&1`;
    $string =~ /BEAST (v\d+\.\d+\.\d+)/;
    return version->parse($1) || undef;
}

sub _alignment {
	my ( $self, $thing, $format ) = @_;
	if ( $thing ) {
		if ( -e $thing ) {
			$self->{'_alignment'} = parse(
				'-format'     => $format,
				'-file'       => $thing,
				'-as_project' => 1,
			);
			#$self->_validate;
		}
		elsif ( ref $thing ) {
			$self->{'_alignment'} = $thing;
			#$self->_validate;
		}
	}
	return $self->{'_alignment'};
}

sub _validate {
	my $self = shift;
	my $project = $self->{'_alignment'};
	my ($taxa) = @{ $project->get_items(_TAXA_) };
	my %taxon = map { $_->get_id => $_ } @{ $taxa->get_entities };
	my $ntax = $taxa->get_ntax;
	
	# add missing rows to sparse matrices
	MATRIX: for my $matrix ( @{ $project->get_items(_MATRIX_) } ) {
		my @rows = @{ $matrix->get_entities };
		#next MATRIX if $ntax == scalar @rows;
		my $nchar = $matrix->get_nchar;
		my %seen;
		$seen{ $_->get_taxon->get_id }++ for @rows;
		my @missing = grep { ! $seen{$_} } keys %taxon;
		my $to = $matrix->get_type_object;
		for my $tid ( @missing ) {
			$matrix->insert( $fac->create_datum( 
				'-type_object' => $to,
				'-name'        => $taxon{$tid}->get_name,
				'-taxon'       => $taxon{$tid},
				'-char'        => ( 'N' x $nchar ),
			) );
		}
	}
}

sub _escape {
	my $string = shift;
	$string =~ s/ /_/g;
	return $string;
}

=back

=cut

1;
