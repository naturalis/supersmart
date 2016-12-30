package Bio::SUPERSMART::Service::ParallelService;
use strict;
use warnings;
use Carp 'croak';
use POSIX 'ceil';
use base 'Bio::SUPERSMART::Service';
use Bio::Phylo::Util::Exceptions 'throw';

my $mode;
my $num_workers;
my $logger = Bio::Phylo::Util::Logger->new;

=over

=item END

MPI_Finalize needs to be invoked exactly once (per mpirun execution) after all parallel 
jobs have finished, therefore this is done in the END block.

=cut

END {
	if ( $mode ) {
		if ( $mode eq 'mpi' ) {
			MPI_Finalize();
		}
	}
}

=item import

When modules are 'use'd, the 'import' sub is automatically called as a package method,
i.e. ParallelService->import. If the 'use' statement has additional arguments, these
are passed into the import statement. For example:

 use ParallelService 'mpi';

Results in the import sub receiving 1. package name, 2. a string determining the mode. 
The mode can either be 'mpi', 'pthreads', or 'pfm'. If none is given, 'pthreads' is the 
default. If perl doesn't have threading support, 'pfm' becomes the default.

When importing this module, the subroutines 'pmap' and 'sequential' are exported to the 
caller's namespace, including thir prototypes (&@ or &). Pmap can be used as a drop-in 
replacement of the built-in perl function map. The method 'sequential' accepts a code
block which is evaluated only on the master node. Note that 'sequential' does not need 
to be used when in 'pthread' mode, but it is still present for convenience. 
 
The number of processes is determined automatically in this class, if e.g. the script 
that uses the ParallelService class is invoked with 'mpirun -np 5', the number of 
processes is set to 4 (4 worker nodes and a master node). In 'pthread'
mode, the number of processes is set to the number of cores present on the machine.

=cut

sub import {
	my $package = shift;
	($mode) = @_;
	
	# no mode specified
	if ( !$mode ) {
		$logger->info("no parallelization mode specified");
		eval { require threads };
		if ( $@ ) {
			$mode = 'pfm';
			$logger->info("no threading support, using fork manager");
		}
		else {
			$mode = 'pthreads';
			$logger->info("threading support available, using pthreads");
		}
	}
	
	# mode is mpi, initialize
	if ( $mode eq 'mpi' ) {
		$logger->info("initializing mpi");
		require Parallel::MPI::Simple;
		Parallel::MPI::Simple->import;
		MPI_Init();
		$num_workers = MPI_Comm_size( MPI_COMM_WORLD() ) - 1;
	}
		
	elsif ( $mode eq 'pthreads' || $mode eq 'pfm' ) {

		# load threads
		if ( $mode eq 'pthreads' ) {
			require threads;
			require threads::shared;
		}
		
		# load fork manager
		else {
			require Parallel::ForkManager;
		}

		# get number of processes from config file
		require Bio::SUPERSMART::Config;
		my $config  = Bio::SUPERSMART::Config->new;
		$num_workers = $config->NODES || 1; 
	} 
	else {
		$num_workers = 1;
	}
	$logger->info("Initializing $mode ParallelService with $num_workers workers");

	my ($caller) = caller();
	eval "sub ${caller}::pmap(\&\@);";
	eval "*${caller}::pmap = \\&${package}::pmap;";

	eval "sub ${caller}::sequential(\&);";
	eval "*${caller}::sequential = \\&${package}::sequential;";

	eval "sub ${caller}::num_workers;";
	eval "*${caller}::num_workers = \\&${package}::num_workers;";

}

=item pmap_pthreads

Does a parallel version of map using pthreads. The function pmap deletages to this 
function if ParallelService was called with 'pthreads' as its argument.

=cut

sub pmap_pthreads (&@) {
	my ( $func, @data ) = @_;
	my $counter = 0;
	my $size    = scalar @data;
	my @threads;
	my @result;
	my $inc = ceil( $size / $num_workers );
	$logger->debug("pmap pthreads has $num_workers nodes available");
	my $thread = 0;

	for ( my $i = 0 ; $i < $size ; $i += $inc ) {
		++$thread;
		my $max = ( $i + $inc - 1 ) >= $size ? $size - 1 : $i + $inc - 1;
		my @subset = @data[ $i .. $max ];
		$logger->debug("Detaching " . scalar(@subset) . " items to thread # " . $thread );
		eval {
			push @threads, threads->create(
				sub {
					map {
						$counter++;
						$logger->debug(
							"Thread $thread is processing item # $counter / "
							. scalar(@subset) );

						# execute code block given in $func with argument
						my $ret = $func->($_);
						$logger->debug("Thread $thread finished processing item # $counter / "
									  . scalar(@subset) );
						$ret;
					} @subset;
				}
				);
		};
		if ( $@ ) {
			$logger->warn("Error in thread $thread");
			throw 'API' => $@;
		}
	}
	for my $i ( 1..scalar(@threads) ) {
		my $thread = $threads[$i-1];
		my @thread_results = $thread->join;
		$logger->debug("Collecting " . scalar(@thread_results) . " results for thread $i");
		push @result, @thread_results;
	}
	return @result;
}

=item pmap_mpi

Does a parallel version of map using MPI. The function pmap deletages to this 
function if ParallelService was called with 'mpi' as its argument. 

=cut

sub pmap_mpi (&@) {
	my ( $func, @data ) = @_;
	my $counter = 0;
	my $WORLD   = MPI_COMM_WORLD();
	my $BOSS    = 0;
	my $JOB     = 1;
	my $RESULT  = 2;
	my $rank    = MPI_Comm_rank($WORLD);

	if ( num_workers() == 0 ) {
		return map {
			$counter++;
			$logger->info( "Worker $rank is processing item # $counter / "
				  . scalar(@data) );
			$func->($_);
		} @data;
	}

	# this is the boss node
	if ( $rank == 0 ) {

		# submit the data in chunks
		my @chunks = map      { [] } 1 .. $num_workers;
		my @idx    = sort map { $_ % $num_workers } 0 .. $#data;

		for my $i ( 0 .. $#data ) {
			my $worker = $idx[ $i - 1 ] + 1;
			push @{ $chunks[ $idx[$i] ] }, $data[$i];
		}

		# send subset of data to worker
		for my $worker ( 1 .. $num_workers ) {
			my $subset = $chunks[ $worker - 1 ];
			$logger->info(
				"Sending " . scalar(@$subset) . " items to worker # $worker" );
			MPI_Send( $subset, $worker, $JOB, $WORLD );
		}

		# receive the result
		my @result;
		for my $worker ( 1 .. $num_workers ) {
			my $subset = $chunks[ $worker - 1 ];
			my $result_subset = MPI_Recv( $worker, $RESULT, $WORLD );
			$logger->info("Received results from worker # $worker");
			push @result, @{$result_subset};
		}
		return @result;
	}
	else {

		# receive my job and process it
		my $subset = MPI_Recv( $BOSS, $JOB, $WORLD );
		my @result = map {
			$counter++;
			$logger->info( "Worker $rank is processing item # $counter / "
				  . scalar( @{$subset} ) );
			$func->($_);
		} @{$subset};

		# send the result back to the boss
		MPI_Send( \@result, $BOSS, $RESULT, $WORLD );
	}
}

=item pmap_pfm

Does a parallel version of map using the parallel fork manager (PFM). The function pmap deletages to this 
function if ParallelService was called with 'pfm' as its argument. 

=cut

sub pmap_pfm {
	my ( $func, @all_data ) = @_;

	my $pm = Parallel::ForkManager->new( num_workers() );

	my $counter = 0;

	# store the results coming from the single threads
	my @all_results;
	my $total = scalar(@all_data);
	# tell the fork manager to return values (references) to the main process
	$pm->run_on_finish(
		sub {
			my ( $pid, $exit_code, $ident, $exit_signal, $core_dump,
				$data_structure_reference )
			  = @_;
			push @all_results, @$data_structure_reference;
		}
	);

	foreach my $data (@all_data) {		
		$counter++;
		
		my $pid = $pm->start and next;
		my @res = ();

		$logger->info( "Processing item $counter of $total");			
		@res = &$func( $data );
		$pm->finish( 0, \@res );    # Terminates the child process		
	}

	$pm->wait_all_children;

	return @all_results;

}

=item sequential

This wrapper allows code to be executed in sequential mode. This is needed
if a script is called using 'mpirun' but one would like to be able to
run some code sequentially and not as many times as passed to mpirun's '-np' 
argument. Code can be simply run sequentially as follows:

my $something = sequential { some code here };

When using 'pthreads' for parallelization, the script is called
with the normal 'perl' command and thus the code in it would run sequencially
anyways unless something is sent to the ParallelService (via pmap). The function
'sequential' is therefore not needed. Thus, in mode 'pthreads', this prototype
simply evaluates the code given as argument.

=cut

sub sequential (&) {
	my ($func) = @_;
	croak "Need mode argument!" if not $mode;
	if ( $mode eq 'mpi' ) {

		# only allow master node to perform operation
		my $WORLD = MPI_COMM_WORLD();
		my $rank  = MPI_Comm_rank($WORLD);
		if ( $rank == 0 ) {
			$func->();
		}
	} else {
		# simply call the function
		$func->();
	}
}

=item pmap
    
This is a wrapper around whichever implementation the user has indicated. Note the
usage of 'goto'. Since the seminal work of Dutch computer scientist Edsger Dijkstra,
we all know that "go to" statements should be considered harmful
(for the original paper from 1968: http://www.u.arizona.edu/~rubinson/copyright_violations/Go_To_Considered_Harmful.html),
however, in this case it should be permissible because it is obviously labelled where
we delegate to, and it serves the purpose of erasing pmap() from the call stack, 
directly forwarding @_ to the actual handler.

=cut

sub pmap (&@) {
	croak "Need mode argument!" if not $mode;
	if ( $mode eq 'pthreads' ) {
		goto &pmap_pthreads;
	}
	elsif ( $mode eq 'mpi' ) {
		goto &pmap_mpi;
	}
	elsif ( $mode eq 'pfm' ) {
		goto &pmap_pfm;
	} 
	elsif ( $mode eq 'mock' ) {
		goto &pmap_mock;
	}
}

=item num_workers

Getter/setter for the number of worker nodes.

=cut

sub num_workers {
	$num_workers = shift if @_;
	return $num_workers;
}

=item mode

Getter/setter for the parallel mode

=cut

sub mode {
	$mode = shift if @_;
	return $mode;	
}

=item distribute

Given an array of data which is assumed to be sorted in
descending "work load", redistributes this work load into
an array of data optimized for parallel processing using
the available number of worker nodes.

=cut

sub distribute {
    my ( $self, @data ) = @_;
	my @result;
    my $nworkers = num_workers();
    if ( scalar @data >= $nworkers ) {
        my @subset; 
        for my $i ( 0 .. $#data ) {
            my $j = $i % $nworkers;
            $subset[$j] = [] if not $subset[$j];
            push @{ $subset[$j] }, $data[$i];
        }       
        push @result, @{ $subset[$_] } for 0 .. ( $nworkers -1 );
    }
    else {
        @result = @data;
    }
	return @result;
}

=back

=cut

1;
