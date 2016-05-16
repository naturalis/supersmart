package Bio::SUPERSMART::App::SubCommand;

use Cwd;
use Term::ANSIColor ':constants';
use Bio::SUPERSMART::Config;
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::Util::Exceptions 'throw';

=head1 NAME

Bio::SUPERSMART::App::SubCommand Superclass for C<smrt> subcommands

=head1 DESCRIPTION

This class provides global options (e.g. verbosity, working directory) that are
inherited by all subcommand classes.

=cut

my $logfh;

=head1 METHODS

All methods below are inherited by child classes.

=over

=item logger

Getter/setter for the subcommand's L<Bio::Phylo::Util::Logger> object, 
verbosity can be de- or increased by  the argument -v for all subcommands. 

=cut

sub logger {
	my $self = shift;
	if ( @_ ) {
		$self->{'logger'} = shift;		
	}
	return $self->{'logger'};
}

=item outfile

Getter/Setter the output file name of the subcommand. Returns absolute path.

=cut

sub outfile {
	my $self = shift;
	if ( @_ ) {
		$self->{'outfile'} = $self->absolute_path(shift);		
	}
	return $self->{'outfile'};
}

=item workdir

Getter/setter for the working directory for the subcommand

=cut

sub workdir {
	my $self = shift;

	if ( @_ ) {
		my $wd = shift;
		$self->logger->debug ("Setting working directory to $wd" );
		$self->{'workdir'} = $wd;		
		if ( !  ( -e $wd && -w $wd ) ) {
			#$self->logger->fatal ("Working directory $wd does not exist or is not writable") ;		
			#exit(1);
			mkdir $wd or die $!;
		}
	}
	return $self->{'workdir'};
}

=item config

Getter/setter for the L<Bio::SUPERSMART::Config> object

=cut

sub config {
	my $self = shift;
	if ( @_ ) {
		$self->{'config'} = shift;
	}
	return $self->{'config'};
}

=item absolute_path

creates an absolute path to a file location, if file is not yet
given as absolute path. This is done by prepending the working directory
to the file name.

=cut

sub absolute_path {
	my $self = shift;
	my $filename = shift;
	if (not $filename =~ /\//) {
		if ( my $wd = $self->workdir ) {	
			$filename = $filename =~ /^\// ? $filename : $wd . "/" . $filename;
		} else {
			$self->logger->warn("no working directory specified, using relative paths ")
		}
	}
	return $filename; 
}

=item execute

Overrides the default classes' 'execute' method such that some global properties 
for all subcommands (e.g. verbosity level, ...) can be set in this method. The
method then calls the 'run' subroutine, which must be implemented for all child classes

=cut

sub execute {
	my ($class, $opt, $args) = @_;
	my $c = Bio::SUPERSMART::Config->new;
	$class->logger->info("This is SUPERSMART release " . $c->RELEASE);
	
	my $result = $class->run( $opt, $args );	
	close $logfh;
	return $result;
}

=item init

The init subroutine initializes objects that are shared by all subcommands (as for instance the working 
directory, output file and a L<Bio::Phylo::Util::Logger> instance). Also, we check for command line options that require a file argument
and, if necessary, turn their paths into absolute paths.

=cut

sub init {
	my ($self, $opt, $args) = @_;
	my $verbosity = INFO;
	$verbosity += $opt->verbose ? $opt->verbose : 0;
 
 	# create logger object with user-defined verbosity
	$self->logger( Bio::Phylo::Util::Logger->new(
		'-level' => $verbosity,
		'-style' => $opt->logstyle,
		'-class' => [ 
		    ref( $self ), 
		    'Bio::SUPERSMART::App::SubCommand', 
		    'Bio::SUPERSMART::Domain::MarkersAndTaxa', 
		    'Bio::SUPERSMART::Service::ParallelService', 
		    'Bio::SUPERSMART::Service::TreeService', 
		    'Bio::SUPERSMART::Service::CalibrationService',
		    'Bio::SUPERSMART::Service::SequenceGetter',
		    'Bio::SUPERSMART::Service::MarkersAndTaxaSelector',
		    'Bio::SUPERSMART::Service::DecorationService',
		    'Bio::Tools::Run::Phylo::ExaML',
		    'Bio::Tools::Run::Phylo::ExaBayes',
		    'Bio::Tools::Run::Phylo::StarBEAST',
		],		
	));
   
	# set working directory
	my $wd = $opt->workdir || '.';
	$self->workdir($wd);

	# loop through options to see which ones are file options; 
	# set absolute path for all filenames given
	my %file_opts = map { (my $optname= $_->[0])=~s/\|.+//g; $_->[2]->{'arg'} eq 'file' ? ($optname=>1) : () } $self->options;
	$file_opts{'logfile'} = 1;
	for my $given_opt ( keys %$opt ) {
		if ( exists $file_opts{$given_opt} ) {
			my $new_filename = 
			$opt->{$given_opt} = $self->absolute_path($opt->{$given_opt});
		}		
	}
	
 	# set outfile name
	if ( my $of = eval { $opt->outfile } ) {
		$self->outfile($of);
	}
     
    # create config object
    $self->config( Bio::SUPERSMART::Config->new );
    
    # redirect logger output to file if specified by command call
    if ( my $logfile = $opt->logfile ) {
    	open $logfh, '>', $logfile or die $!;
		$self->logger->set_listeners(sub{$logfh->print(shift)});
    } 
}

=item opt_spec

Speciefies the global options for the smrt command. This subroutine will
also call the 'options' for the child classes to capture also the
subcommand specific options.

=cut 

sub opt_spec {
	my ($class, $app) = @_;
	( my $subcommand = lc $class ) =~ s/.*:://g ;
	return (
		$class->options($app),
		[ "help|h", "display help screen", { urgency => 'super' } ],
		[ "verbose|v+", "increase verbosity level", { urgency => 'super' } ],
		[ "workdir|w=s", "directory in which results and intermediate files are stored", { arg => "dir", urgency => 'super' } ],
		[ "logfile|l=s", "write run-time information to logfile", 
		  { arg => "file", urgency => 'super', galaxy_out => 1, galaxy_label => "logfile-${subcommand}", galaxy_type => 'text'}],	
		[ "logstyle|y=s", "toggles logging style between 'simple' and 'detailed'", { default => "simple", urgency => 'super' }],
	);	
}

=item description

Returns the description for a specific subcommand. This method is invoked by 
the child classes when the global option 'help' is used, e.g. in 
$ smrt help subcommand
The actual description of the subcommand is parsed from its DESCRIPTION field
in the POD documentation.

=cut

sub description {
	my ($class, $opt, $args) = @_;

	# classname to filename
  	my $pm_file = $class;
  	$pm_file =~ s/::/\//g; 
  	$pm_file =~ s/=.*//g;   	
  	$pm_file .= '.pm';
  	$pm_file = $INC{$pm_file} or return "(unknown)";
	
  	# if the pm file exists, open it and parse it
  	open my $fh, "<", $pm_file or return "(unknown)";
	my $seen;
	my $description = '';
	LINE: while( <$fh> ) {
#		chomp;
		if ( /DESCRIPTION/ ) {
			$seen++;
			next LINE;
		}
		if ( /^=cut/ and $seen ) {
			$seen--;
		}
		if ( $seen and not $description and not /\S/ ) {
			next LINE;
		}
		if ( $seen and not $descripton and /\S/ ) {
			$description .= $_;
			next LINE;
		}
		$description .= $_ if $seen and $description;
	}
	return $description;
}

=item validate_args

Overrides the default method mainly for one reason: There are two ways of 
displaying the help screen: One is '$ smrt help command', the other is
'$ smrt command -h'. Here we make the latter to behave exactly as
the former.

=cut

sub validate_args {
	my ($class, $opt, $args) = @_;		
	
	# first init the properties of the superclass
	$class->init($opt, $args);
			
	if ($opt->help){

		# This is a small hack to make the help screen called with the 
		#  option (as here, e.g. "$ smrt command -h" ) look the same
		#  as when we call "$ smrt help command"  				
		my ($cmd) = $class->app->prepare_command("help"); 
		my @name = ($class->command_names);
		$cmd->execute($opt, \@name);
		exit;
	} 
	else {
		$class->validate($opt, $args);
	}
}

=item usage_desc

Overrides the method from L<App::Cmd> to compile a custom 'usage' 
string for a given subcommand with options and arguments 
of the form 'smrt subcommand option1 <arg1> [option2 <arg2>] [option3]'.  
Options that require arguments followed by the argument name enclosed in '<>'.
Options that are not mandatory are enclosed in []. 

=cut 

sub usage_desc {
	my ($class, $opt, $args) = @_;	
	my $cmd = $class->command_names;	
	my @opts = $class->opt_spec;			
	my $usage = "%c $cmd ";
		
	# build custom usage string (default was "%c $cmd %o" )
	my @args = ( '' );
	@opts = sort { ( $a->[2]->{urgency} eq 'super' ) <=> ( $b->[2]->{urgency} eq 'super' ) } @opts;
	for my $opt (@opts){
		my $s = $opt->[0];
		my %h = %{ $opt->[2] };
		my ($short_opt) = ($s =~ /\|([a-z]+)/);
		my $required = $h{'mandatory'};
		my $arg      = $h{'arg'};
		my $opt_str  = "-$short_opt";
		$opt_str    .= $arg ? " <$arg>"  :"";
		$opt_str     = !$required ? "[$opt_str] " : " $opt_str ";		

		# encode colored
		if ( $h{'urgency'} =~ /super/ ) {
			$args[-1] .= "\033[0;37m" . $opt_str . "\033[0m";
		}
		elsif ( $required ) {
			$opt_str =~ s/^\s*//;
			$opt_str =~ s/\s*$//;
			$args[-1] .= " \033[4m" . $opt_str . "\033[0m ";
		}
		else {
			$args[-1] .= $opt_str;
		}
		if ( length($args[-1]) > 60 ) {
			push @args, '';
		}
	}
	return $usage . join( "\\\n\t", grep { /\S/ } @args );	
}

=item usage_error

Subclasses the method from L<App::Cmd> to display the error message using
the logger before delegating to the super class.

=cut

sub usage_error {
	my ( $self, $message ) = @_;
	$self->logger->error($message);
	die "Usage: " . $self->_usage_text . "\n";
}

=item options

This is a wrapper for L<App::Cmd>'s 'opt_spec' method and is implemente by the child classes to
specify command-line options for a specific subcommand.

=cut

sub options {
	throw 'NotImplemented' => "subroutine 'options' not implemented by " . ref(shift);
}

=item validate

This is a wrapper for L<App::Cmd>'s 'validate_args' method and is implemente by the child classes to
validate command-line options for a specific subcommand.

=cut

sub validate {
	throw 'NotImplemented' => "subroutine 'validate' not implemented by " . ref(shift);

}

=item run

This is a wrapper for L<App::Cmd>'s 'execute' method and is implemente by the child classes to
run a specific subcommand
=cut

sub run {
	throw 'NotImplemented' => "subroutine 'run' not implemented by " . ref(shift);
}

=back

=cut

1;
