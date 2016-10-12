package Bio::Tools::Run::Phylo::ExaBayes;
use strict;
use version;
use Cwd;
use File::Spec;
use File::Temp 'tempfile';
use Bio::Phylo::Generator;
use Bio::Phylo::IO 'parse';
use Bio::Tools::Run::Phylo::PhyloBase;
use Bio::Phylo::Util::CONSTANT ':objecttypes';
use Bio::Phylo::Util::Logger;
use Bio::SUPERSMART::Service::TreeService;

use base qw(Bio::Tools::Run::Phylo::PhyloBase);

my $treeservice = Bio::SUPERSMART::Service::TreeService->new;
my $log = Bio::Phylo::Util::Logger->new;

# program name required by PhyloBase
our $PROGRAM_NAME = 'exabayes'; 

our @ExaBayes_PARAMS = (
        'f',        # alnFile:    alignment file (binary and created by parser or plain-text phylip)
        's',        # seed:       master seed for the MCMC
        'n',        # ruid:       a run id
        'r',        # id:         restart from checkpoint id
        'q',        # modelfile:  RAxML-style model file
        't',        # treeFile:   file containing starting trees (newick format) for chains
        'm',        # model:      data type for single partition non-binary alignment file (DNA | PROT)
        'c',        # confFile:   file configuring your ExaBayes run
        'w',        # dir:        working directory for output files
        'R',        # num:        number of runs (i.e., independent chains) to be executed in parallel
        'C',        # num:        number of chains (i.e., coupled chains) to be executed in parallel
        'M',        # mode:       memory versus runtime trade-off, value 0 (fastest) to 3 (most memory efficient)        
	    'T',        # x           start x threads per MPI process. If you do not use MPI, simply start x threads. 
    );

our @ExaBayes_CONFIGFILE_PARAMS = (
        # RUN parameters
        'numRuns',           # number of independant runs
        'numGen',            # number of generations to run. Defaults to 1e6
        'parsimonyStart',    # if 1, start from most parsimoneous tree
        #MCMC parameters
        'numCoupledChains'   # number of chains per independent run
    );

our @Consense_PARAMS = (
        # The following command line parameters are not for "ExaBayes" per se 
        # but for the tool "Consense" which is also handled in this class.
        't',        # thresh:     threshold for the consenus tree 
                    #             values between 50 (majority rule) and 100 (strict) or MRE (the greedily refined MR consensus).  
                    #             Default: MRE
        'b',        # relBurnin:  proportion of trees to discard as burn-in (from start). Default: 0.25        
    );

our @Misc_PARAMS = (
        # The following parameters are artificial, not a 'real' parameter of ExaBayes nor Consense
        'outfile_name',        # file name of the inferred consensus trees
        'outfile_format',       # format of output file. Either 'nexus' or 'newick' 
        'consense_bin',         # location of consense binary 
    );


our @ExaBayes_SWITCHES = (
        'v',        # version
        'z',        # quiet mode
        'd',        # dry run
        'Q',        # per-partition data distribution
        'S'         # try to save memory using the SEV-technique for gap columns       
    );


# Create function aliases for (one letter) command line arguments and switches

=over

=item work_dir

Getter/setter. Intermediate files will be written here.

=cut

*work_dir = *w;

=item run_id

Getter/setter for an ID string for this run. by default based on the PID. this only has
to be unique during the runtime of this process, all files that have the run ID in them
are cleaned up afterwards.

=cut

*run_id = *n;

=item seed

Getter/setter for master seed for MCMC.

=cut

*seed = *s;

=item burnin_fraction

Getter/setter for fraction of burn in (used in the consense program).

=cut

*burnin_fraction = *b;

=item config_file

Getter/setter for a config file for the ExaBayes run.

=cut

*config_file = *c;

=item model

Getter/setter for Exabayes model (DNA | Protein)

=cut

*model = *m;

=item mode

Getter/setter for Exabayes Mode (performance vs memory usage, value between 0 and 3)

=cut

*mode = *M;

=item alnFile

Getter/Setter for phylip alignment file 

=cut

*alnFile = *f;

=item treeFile

Getter/setter for file with starting tree

=cut

*treeFile = *t;

=item threads

Getter/setter for number of threads

=cut

*threads = *T;

=item dryRun

Getter/setter for whether or not just a dry run should be performed

=cut

*dryRun = *d;

=item restart

Getter/setter to restart from checkpoint (set previous run ID)

=cut

*restart = *r;

##sub program_name { $PROGRAM_NAME }

=item program_dir

(no-op)

=cut

sub program_dir { undef }


=item version

Returns the version number of the ExaBayes executable.

=cut

sub version {
    my ($self) = @_;
    my $exe;
    return undef unless $exe = $self->executable;
    my $string = `$exe -v 2>&1`;
    $string =~ /ExaBayes, version (\d+\.\d+\.\d+)/;
    return version->parse($1) || undef;
}

=item new

This constructor takes the following mandatory arguments: 

    -f alnFile       a alignment file (either binary and created by parser or plain-text phylip)
    -s seed          a master seed for the MCMC
    -n ruid          a run id
    -r id            restart from checkpoint. Just specify the id of the previous run (-n) here. 
                       Make sure, that all files created by the previous run are in the working directory.
                       This option is not mandatory for the start-up, seed (via -s) will be ignored.
    -q modelfile     a RAxML-style model file (see manual) for multi-partition alignments. Not needed 
                       with binary files.
    -t treeFile      a file containing starting trees (in Newick format) for chains. If the file provides less
                       starting trees than chains to be initialized, parsimony/random trees will be used for
                       remaining chains. If a tree contains branch lengths, these branch lengths will be used
                       as initial values.
    -m model         indicates the type of data for a single partition non-binary alignment file
                       (valid values: DNA or PROT)

Optional arguments are:

    -d               execute a dry-run. Procesess the input, but does not execute any sampling.
    -c confFile      a file configuring your ExaBayes run. For a template see the examples/ folder
    -w dir           specify a working directory for output files
    -R num           the number of runs (i.e., independent chains) to be executed in parallel
    -C num           number of chains (i.e., coupled chains) to be executed in parallel
    -Q               per-partition data distribution (use this only with many partitions, check manual
                       for detailed explanation)
    -S               try to save memory using the SEV-technique for gap columns on large gappy alignments
                       Please refer to  http://www.biomedcentral.com/1471-2105/12/470
                       On very gappy alignments this option yields considerable runtime improvements. 
    -M mode          specifies the memory versus runtime trade-off (see manual for detailed discussion).
                       <mode> is a value between 0 (fastest, highest memory consumption) and 3 (slowest,
                       least memory consumption)

=cut

sub new {
        my ( $class, @args ) = @_;
        my $self = $class->SUPER::new(@args);
        $self->_set_from_args(
                \@args,
                '-methods' => [ @ExaBayes_PARAMS, @ExaBayes_SWITCHES, @ExaBayes_CONFIGFILE_PARAMS, @Misc_PARAMS, @Consense_PARAMS ],
                '-create'  => 1,
            );
        my ($out) = $self->SUPER::_rearrange( [qw(OUTFILE_NAME)], @args );
        $self->outfile_name( $out || '' );
        return $self;
}

=item run

Runs an analysis with ExaBayes. Input is an alignment file in phylip format.
Returns the file location of the consense tree.

=cut

sub run {
	my $self = shift;
	my $ret;
	my %args = @_;
	my $phylip  = $args{'-phylip'} || die "Need -phylip arg";
	
	my $binary = $self->run_id . '-dat' ;
	$binary = $treeservice->make_phylip_binary( $phylip, $binary, $self->parser, $self->work_dir );
	$binary = File::Spec->catfile($self->work_dir, $binary);
	
	$self->alnFile($binary);
	
	# compose argument string: add MPI commands, if any
	my $string;
	if ( $self->mpirun && $self->nodes ) {
		$string = sprintf '%s -np %i ', $self->mpirun, $self->nodes;
	}
	
	$string .= $self->executable . $self->_setparams;       
	
	my $current_dir = getcwd; 
	if ( $self->work_dir ) {
		chdir($self->work_dir);
	}
	
	# run exabayes
	$log->info("going to run '$string'");        
	system($string) and $self->warn("Couldn't execute '$string': $! (errno: $?)");
	
	$ret = $self->_run_consense;
	
	chdir $current_dir;
	
	#return $self->_cleanup;
	# return outfile name or 1 if it was a dry run
	return $ret || $self->dryRun;
}

sub _run_consense {
        my $self = shift;
        my $ret = "";
        
        # add all topology files created by Exabayes. CAUTION! Some 'magic words' in here
        my $exabayesout = $self->work_dir . "/ExaBayes_topologies." . $self->run_id . ".*";
        my @arr = glob $exabayesout; 
        if (scalar( @arr ) < 1 ) {
                $log->warn("cannot run consense since ExaBayes did not produce output. Did you perform a dry run?");
        }
        else {
                # build command string for running 'consense'
                my $string = $self->consense_bin . ( $self->_set_consense_params || '' );
                $string .= " -f " . $exabayesout . " -n " . $self->run_id;
                
                # run consense
                $log->info("Going to run $string");        
                my $output = system($string) and $self->warn("Couldn't run Consense: $?");
                
                # get format for output file, if not given, chose 'newick'
                my $format = $self->outfile_format || "newick";
                $format = ucfirst($format);
                
                # consense output file name is not controllable, therefore rename consense outfile
                my $filename = "ExaBayes_Consensus*".$format . "." . $self->run_id;
                
                if ( $self->outfile_name ) {
                        my $mvstr = "mv " . $filename . " " . $self->outfile_name;
                        system ($mvstr) and $self->warn("Couldn't run command $mvstr : $?");
                        $ret = $self->outfile_name;
                } 
                else {
                        $ret = $filename;
                }
        }
        return $ret;
}


=item mpirun

Getter/setter for the location of the mpirun program for parallelized runs. When unset,
no parallelization is attempted.

=cut

sub mpirun {
        my ( $self, $mpirun ) = @_;
        if ( $mpirun ) {
		$self->{'_mpirun'} = $mpirun;
	}
	return $self->{'_mpirun'};
}

=item nodes

Getter/setter for number of mpi nodes

=cut

sub nodes {
	my ( $self, $nodes ) = @_;
	if ( $nodes ) {
		$self->{'_nodes'} = $nodes;
	}
	return $self->{'_nodes'};
}

=item parser

Getter/setter for the location of the parser program (which creates a compressed, binary 
representation of the input alignment) that comes with ExaBayes. 

=cut

sub parser {
	my ( $self, $parser ) = @_;
	if ( $parser ) {
		$self->{'_parser'} = $parser;
	}
	return $self->{'_parser'} || 'parser';
}


sub _write_config_file {
        my $self = shift;
        my $conffile = File::Spec->catfile( $self->work_dir, $self->run_id . '.nex' );
	open my $conffh, '>', $conffile or die $!;
        print $conffh "#NEXUS\n";
        print $conffh "begin run;\n";
        for my $attr (@ExaBayes_CONFIGFILE_PARAMS) {
                my $value = $self->$attr();
                next unless defined $value;
                print $conffh $attr . "\t" . $value . "\n";
        }
        print $conffh "end;\n";
        close $conffh;
        $self->config_file($conffile);
        return $self->config_file;
}

sub _setparams {
    my $self = shift;
    my $param_string = '';

    # add config file to parameters if not already existant
    if (! $self->config_file){
            $self->_write_config_file;
    }
    
    # iterate over parameters and switches
    for my $attr (@ExaBayes_PARAMS) {
        my $value = $self->$attr();
        next unless defined $value;
        $param_string .= ' -' . $attr . ' ' . $value;
    }
    for my $attr (@ExaBayes_SWITCHES) {
        my $value = $self->$attr();
        next unless $value;
        $param_string .= ' -' . $attr;
    }
    
    # hide stderr
    my $null = File::Spec->devnull;
    $param_string .= " > $null 2> $null" if $self->quiet() || $self->verbose < 0;

    return $param_string;
}

sub _set_consense_params {
        my $self = shift;
        my $param_string;
        for my $attr (@Consense_PARAMS) {
                        my $value = $self->$attr();
                        next unless defined $value;
                        $param_string .= ' -' . $attr . ' ' . $value;                        
        }
        # hide stderr
        my $null = File::Spec->devnull;
        $param_string .= " > $null 2> $null" if $self->quiet() || $self->verbose < 0;
        return $param_string;
}

=back

=cut

1;
