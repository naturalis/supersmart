package Bio::SUPERSMART::Service::InferenceService::exabayes;
use Bio::Tools::Run::Phylo::ExaBayes;
use Bio::SUPERSMART::Service::InferenceService;
use base 'Bio::SUPERSMART::Service::InferenceService';

=head1 NAME

Bio::SUPERSMART::Service::InferenceService::exabayes - Infers phylogenetic trees

=head1 DESCRIPTION

Provides functionality inferring trees and sets of trees by wrapping ExaBayes.

=over

=item configure

Provides the $config object to the wrapper so that settings defined in $config object can 
be applied to the wrapper.

=cut

sub configure {
    my ( $self, %args ) = @_;
    my $tool    = $self->wrapper;
    my $config  = $self->config;
    my $logger  = $self->logger;
    my $outfile = $self->outfile;
        
	# set checkpoint ID from previous run, if given
	if ( my $file = $args{'checkpoint'} ) {
		$logger->info("Restarting from checkpoint file $file");
		# extract previous run ID from checkpoint file named ExaBayes_checkpoint.RUNID
		(my $checkpoint = $file) =~ s/ExaBayes_checkpoint\.//g;
		$tool->restart($checkpoint);
	}

    # set outfile name
    $logger->info("going to create output file $outfile");
    $tool->outfile_name($outfile);    

    # set mpirun location
    $logger->info("going to use mpirun executable ".$config->MPIRUN_BIN);
	$tool->mpirun($config->MPIRUN_BIN);
    
    # set exabayes location
    $logger->info("going to use executable ".$config->EXABAYES_BIN);
    $tool->executable($config->EXABAYES_BIN);
        
    # set parser location
    $logger->info("going to use parser executable ".$config->EXABAYES_PARSER_BIN);
    $tool->parser($config->EXABAYES_PARSER_BIN);
        
    # set number of independent runs
    $logger->info("going to perform ".$config->EXABAYES_NUMRUNS." independent runs");
    $tool->numRuns($config->EXABAYES_NUMRUNS);

	# set number of runs to run in parallel
	$logger->info("running ".$config->EXABAYES_NUMRUNS." independent runs in parallel");
    $tool->R($config->EXABAYES_NUMRUNS);

    # set number of coupled chains
    $logger->info("setting coupled chains to ".$config->EXABAYES_NUMCHAINS);
    $tool->numCoupledChains($config->EXABAYES_NUMCHAINS);

	# set number of coupled chains to run in parallel
	$logger->info("running ".$config->EXABAYES_NUMCHAINS." coupled chains in parallel");
	$tool->C($config->EXABAYES_NUMCHAINS);
	
    # set number of mpi processes; increase when numruns * numchains too large
	# to prevent ExaBayes error
	my $nodes = $config->NODES;
	if ( $config->EXABAYES_NUMRUNS * $config->EXABAYES_NUMCHAINS > $nodes ) {
		$nodes = $config->EXABAYES_NUMRUNS * $config->EXABAYES_NUMCHAINS;
	}
	$logger->info("setting number of MPI nodes to $nodes");
	$tool->nodes($nodes);

    # set number of generations
    $logger->info("setting number of generations to ".$config->EXABAYES_NUMGENS);
    $tool->numGen($config->EXABAYES_NUMGENS);
            
    # set random seed to config value
    $tool->seed($config->RANDOM_SEED); 
    
    # set binary of consense program
    $logger->info("setting consense binary");
    $tool->consense_bin($config->EXABAYES_CONSENSE_BIN);
        
    # set burnin fraction
    $logger->info("setting burnin fraction to ".$config->BURNIN);
    $tool->burnin_fraction($config->BURNIN); 

	# set starting tree
	if ( my $tree = $self->usertree ) {
		# user defined
		$logger->info("Setting starting tree $tree");
		$tool->treeFile($tree);
	}
	else {
		# parsimony start
		$logger->info("setting ExaBayes to start from parsimony tree");
		$tool->parsimonyStart('true');
	}
	
}

=item create

Instantiates and configures the wrapper object, returns the instantiated wrapper.

=cut

sub create {
    my $self = shift;

    # instantiate helper objects and variables
    my $logger  = $self->logger;
    my $workdir = $self->workdir;
    my $tool    = Bio::Tools::Run::Phylo::ExaBayes->new;
    
    # set working directory
    $logger->info("going to use working directory $workdir");
    $tool->work_dir($workdir);        
	
    # set run id to current process id
    $tool->run_id("infer_backbone_".$$);
        
    # set oufile format
    $logger->info("setting outfile format to 'newick'");
    $tool->outfile_format("newick");
        
    # set mode to zero (faster, but takes more memory)
    $logger->info("setting exabayes 'mode' argument to zero ");
    $tool->mode(0);             
    return $tool;
}

=item run

Runs the analysis on the provided supermatrix. Returns a tree file.

=cut

sub run {
    my ( $self, %args ) = @_;
    my $consensus = $self->wrapper->run( '-phylip' => $args{'matrix'} );
    return $self->interdigitate;
}

=item interdigitate

Combines the set of topology (nexus) files into a single newick file such 
that the trees from the different chains are interdigitated so that the
correct burnin is discarded should a consense operation be applied to them.

=cut

sub interdigitate {
    my $self = shift;
    my $runid = $self->wrapper->run_id;
    my $dir = $self->workdir;    
    my $log = $self->logger;
    
    # assemble set of topology files
    my @files;
    opendir my $dh, $dir or die $!;
    while( my $entry = readdir $dh ) {
        if ( $entry =~ /ExaBayes_topologies.${runid}.\d+/ ) {
            push @files, "${dir}/${entry}";
        }
    }
	
    # iterate over files
    my @tree_sets;	
    for my $file ( @files ) {
        my %translate;
        my @trees;
        my $tflag = 0;		
        open my $fh, '<', $file or die $!;
        LINE: while(<$fh>) {
            chomp;
			
            # capture translation table
            $tflag++ and next LINE if /translate/;
            if ( $tflag and /^\s+(\d+)\s+(\d+)([,;])/ ) {
                my ( $id, $name, $token ) = ( $1, $2, $3 );
                $translate{$id} = $name;
                if ( $token eq ';' ) {
                    $tflag--;
                }
            }
			
            # capture trees
            if ( /tree (.+?) = \[&U\] (.+?;)/ ) {
                my ( $name, $tree ) = ( $1, $2 );
				
                # capture the starting positions of all IDs in the tree string
                my %pos;
                while( $tree =~ /(\d+):/g ) {
                    my $id = $1;
                    $pos{$id} = pos($tree) - ( length($id) + 1 );
                }
				
                # replace IDs with names, working from right to left
                for my $id ( sort { $pos{$b} <=> $pos{$a} } keys %pos ) {
                    substr $tree, $pos{$id}, length($id), $translate{$id};
                }
                push @trees, $tree;
            }
        }
        push @tree_sets, \@trees;	
    }
	
    # merge
    my @interdigitated;		
    my $i = 0;
    TREE: while(1) {
        my $have_trees = 0;
        for my $set ( @tree_sets ) {
            if ( $set->[$i] ) {
                $have_trees++;
                push @interdigitated, $set->[$i];
            }
        }
	$i++;
        last TREE unless $have_trees;
    }

    # for some reason, two spurious initial trees (starting tree?)
    # are still included, so we end up with ngens/samplefreq + 2 
    # trees (e.g. 402). these trees are so bad that rerooting and
    # calibrating often fails, so let's just omit them here
    shift @interdigitated;
    shift @interdigitated;
	
    # write to file
    open my $out, '>', $self->outfile or die $!;
    print $out join "\n", @interdigitated;
    return $self->outfile;
}

=item is_bayesian

Indicates whether the inference service is bayesian, in which case we will
not run a bootstrapping analysis. Returns true.

=cut

sub is_bayesian { 1 }

=item cleanup

Cleans up any intermediate files.

=cut

sub cleanup {
    my $self = shift;
    my $runid = $self->wrapper->run_id;
    my $dir = $self->workdir;    

    # remove the input matrices
    unlink "${dir}/${runid}-dat.binary";
    unlink "${dir}/${runid}.nex";

    # remove single intermediate files
    for my $prefix ( qw(diagnostics checkpoint info prevCheckpointBackup ConsensusExtendedMajorityRuleNexus) ) {
        unlink "${dir}/ExaBayes_${prefix}.${runid}";
    }

    # remove files from chains
    opendir my $dh, $dir or die $!;
    while( my $entry = readdir $dh ) {
        for my $prefix ( qw(parameters topologies) ) {
            if ( $entry =~ /ExaBayes_${prefix}.${runid}.\d+/ ) {
                unlink "${dir}/${entry}";
            }
        }
    }
}

=back 

=cut

1;
