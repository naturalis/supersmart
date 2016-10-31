package Bio::SUPERSMART::App::smrt::Command::BBinfer;
use strict;
use warnings;
use File::Copy qw(copy move);
use File::Temp 'tempfile';
use Bio::Phylo::IO qw(parse parse_tree);
use Bio::Phylo::Factory;
use Bio::Phylo::Util::CONSTANT ':namespaces';
use Bio::SUPERSMART::Config;
use Bio::SUPERSMART::Domain::MarkersAndTaxa;
use Bio::SUPERSMART::Service::TreeService;
use Bio::SUPERSMART::Service::SequenceGetter;
use Bio::SUPERSMART::Service::InferenceService;
use Bio::SUPERSMART::Service::MarkersAndTaxaSelector;
use Bio::SUPERSMART::Service::ParallelService;
use Bio::Phylo::Util::Exceptions 'throw';
use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);

# ABSTRACT: inference of a genus-level backbone tree

=head1 NAME

BBinfer.pm - inference of genus-level backbone tree

=head1 SYNOPSIS

smrt bbinfer [-h ] [-v ] [-w <dir>] -s <file> -t <file> [-i <tool>] [-o <file>] [-r ]

=head1 DESCRIPTION

Given an input supermatrix (in interleaved PHYLIP format), infers a backbone tree or set
of trees. In the case of maximum likelihood inference, the set of trees are bootstrapped
replicates (this applies to the tools ExaML and RAxML, currently). For large amounts of data, 
RAxML's fast bootstrap option provides a feasible way to obtain abackbone phylogeny. For ExaBayes, this
the set is a sample of the posterior distribution of trees. The tree(s) resulting from this command 
are written to the specified output file.

=cut

# define additional command line arguments
sub options {
    my ($self, $opt, $args) = @_;       
    my $outfile_default = "backbone.dnd";
    my $tool_default    = "ExaML";
    my $matrix_default  = "supermatrix.phy";
    my $taxa_default    = "species.tsv";
    my $boot_default    = 1;
	my @tools = qw(RAxML ExaML PhyML ExaBayes);
    return (
        [
		     "supermatrix|s=s", 
		     "matrix of concatenated multiple sequence alignments as produced by 'smrt bbmerge'", 
		     { arg => "file", default => $matrix_default, galaxy_in => 1, galaxy_type => 'data' }
		],  
        [
		     "inferencetool|i=s", 
    		 "software tool for backbone inference, defaults to $tool_default. Available tools: " . join(', ', @tools), 
    		 {default => $tool_default, arg => "tool", galaxy_in => 1, galaxy_type => 'select', galaxy_value => $tool_default, galaxy_options => \@tools }
		],
		[
		     "taxafile|t=s", 
			 "[ExaML inference only] taxa file to generate starting tree. Defaults to $taxa_default", 
		     { arg => "file", default => $taxa_default, galaxy_in => 1, galaxy_type => 'data', galaxy_format => 'tabular' }
		],
		[
		     "constrain_genera|c",
		     "[RAxML inference only] constrain genera to be monophyletic, won't work with --random start option",
     		 { galaxy_in => 1, galaxy_type => 'boolean', galaxy_condition => { 'inferencetool' => ['ExaML'] } }
		],
		[
		     "starttree|k=s", 
			 "[ExaML inference only] starting tree for maximum-likelihood inference", 
		     { galaxy_in => 1, galaxy_type => 'data', galaxy_condition => { 'inferencetool' => 'ExaML' } }
		],
		[
		     "random_start|m", 
			 "[ExaML inference only] start inference from random starting tree", 
		     { galaxy_in => 1, galaxy_type => 'boolean', galaxy_condition => { 'inferencetool' => ['ExaML'] } }
		],
        [
		     "bootstrap|b=i", 
    		 "[ExaML and RaXML inference only] number of bootstrap replicates. Will add the support values to the backbone tree.", 
    		 { 
				 default => $boot_default, 
				 galaxy_in => 1, 
				 galaxy_type => 'integer', 
				 galaxy_value => $boot_default, 
				 galaxy_condition => { 'inferencetool' => ["RAxML", "ExaML", "PhyML"] } 
			 }
		],
		[
		     "rapid_boot|r", 
    		 "[RaXML inference only] use RAxML's rapid bootstrap algorithm. Returns single consensus tree with annotated nodes.", 
		     { galaxy_in => 1, galaxy_type => 'boolean', galaxy_condition => { 'inferencetool' => 'RAxML' } }
		],
        [
    		 "ids|n", 
    		 "return tree with NCBI identifiers instead of taxon names", 
		     { galaxy_in => 1, galaxy_type => 'boolean' }
		],
        [
		     "outfile|o=s", 
		     "name of the output tree file (in newick format), defaults to '$outfile_default'", 
		     {default => $outfile_default, arg => "file", galaxy_out => 1, galaxy_type => "data", galaxy_label => "backbone" }
		],
		[    
			 "append|a=s",
			 "[ExaBayes inference only] append output to previous run, need to specify ExaBayes checkpoint file (ExaBayes_checkpoint.infer_backbone_XXXX)",
			 {}
		],
        [    "cleanup|x", 
			 "if set, cleans up all intermediate files", 
			 {}
		],
    );  
}

# validate command line arguments
sub validate {
    my ($self, $opt, $args) = @_;

    my $sm = $opt->supermatrix;
	my $tax = $opt->taxafile;
    my $tool = $opt->inferencetool;
    
    $self->usage_error("Need supermatrix") if not $sm;
    $self->usage_error("File $sm does not exist") unless -e $sm;
    $self->usage_error("File $sm is empty") unless -s $sm;	
	$self->usage_error("Rapid boostrap only supported when 'inferencetool' argument is RAxML.") if ! (lc $tool eq 'raxml') and $opt->rapid_boot;
	$self->usage_error("starttree only supported when 'inferencetool' argument is ExaML or RaXML.") if ! ($tool =~ /^examl$|^raxml$|^phyml$/i) and $opt->starttree;
	if ( lc ($tool)  eq 'examl' ) {
		$self->usage_error("Need taxa file for ExaML inference") unless $opt->taxafile;
		$self->usage_error("File $tax does not exist") unless -e $tax;
		$self->usage_error("File $tax is empty") unless -s $tax;
	}	
}

# run the analysis
sub run {
    my ( $self, $opt, $args ) = @_;     
        
    # collect command-line arguments
    my $supermatrix = $opt->supermatrix;
    my $bootstrap   = $opt->bootstrap;
	
    # instantiate and configure helper objects
    my $ts = Bio::SUPERSMART::Service::TreeService->new;     
    my $ss = Bio::SUPERSMART::Service::SequenceGetter->new;   
    my $is = Bio::SUPERSMART::Service::InferenceService->new(
        'tool'      => lc( $opt->inferencetool ),
        'workdir'   => $opt->workdir,
		'bootstrap' => $opt->bootstrap
    );
   
	# need starting tree for examl inference
	if ( $opt->inferencetool =~ /^examl$|^raxml$|^phyml$/i ) {
		$self->_set_usertree( $is, $supermatrix, $opt->taxafile, $opt->random_start, $opt->starttree );
	}
	
	# constrain genera to be monophyletic if specified by user
	if ( $opt->constrain_genera and $opt->inferencetool eq 'raxml' ) {
		my $cs = 'constraint.dnd';
		$is->constraint_tree($ts->make_constraint_tree($supermatrix, $cs));
	}
	
	# For RaXML's rapid bootstrap, do not create bootstrap matrices since it's taken care
	#  of by RaXML
	if ( $opt->rapid_boot ) {
		$bootstrap = 1;
		$is->rapid_boot(1);
	}
	
    # run the analysis, process results
    my $base = $self->outfile;	
	pmap {
		my $i = $_;
		
		# assign input matrix
		my $matrix;
		if ( $is->is_bayesian or $bootstrap == 1 ) {
			$matrix = $supermatrix;
			}
		else {
			$matrix = $ss->bootstrap( $supermatrix, $i );
		}
		
		# create replicate settings
		$is->outfile( "${base}.${i}" );
		$is->replicate( $i );
		$is->configure('checkpoint' => $opt->append);
		
		# run
		my $backbone = $is->run( 'matrix' => $matrix );  
		$self->_append( $backbone, $base . '.trees' );
		
		# cleanup, if requested
		if ( $opt->cleanup ) {
			unlink $backbone;
			$is->cleanup;
			if ( not $is->is_bayesian and $bootstrap != 1 ) {
				unlink $matrix;
			}
		}
	} ( 1 .. $bootstrap );
	
	# finalize
	$self->_process_result( 
		$base . '.trees',                       # the set of trees (file)
		!$opt->ids,                             # whether to remap
		( $bootstrap - 1 || $is->is_bayesian ), # whether to consense
		$supermatrix,                           # the supermatrix (file)
		$is,                                    # inference service
		);		

    unlink $self->workdir . '/user.dnd' if $opt->cleanup;
    return 1;
}

# append bootstrap replicate results
sub _append {
    my ( $self, $in, $out ) = @_;
    $self->logger->info("Appending contents of $in to $out");

    # read the in file
    open my $infh, '<', $in or die $!;
    my $content = do { local $/; <$infh> };
    close $infh;
    
    # write to the out file
    open my $outfh, '>>', $out or die $!;
    print $outfh $content;
    close $outfh;
}

# process the inferred tree, e.g. by mapping identifiers
sub _process_result {
    my ( $self, $backbone, $remap, $consense, $matrix, $service ) = @_;

    my $ts = Bio::SUPERSMART::Service::TreeService->new;
    my $outfile = $self->outfile;

    # map ID to name
    if ( $remap ) {
		my $ms  = Bio::SUPERSMART::Service::MarkersAndTaxaSelector->new;
		my $mt  = Bio::SUPERSMART::Domain::MarkersAndTaxa->new;		

		my %map = map { $_ => $ms->find_node($_)->taxon_name } $mt->get_supermatrix_taxa($matrix);
        $map{$_} =~ s/ /_/g for keys %map;		
		my ( $fh, $filename ) = tempfile();
		close $fh;
		
		$ts->remap_newick( $backbone => $filename, %map );
		move( $filename => $outfile );
        unlink $backbone;
    }
    
    $self->logger->info("DONE, results written to $outfile");
}

# set usertree to inference service object if applicable.
# currently this can be done when inference tool is examl
sub _set_usertree {
	my ( $self, $service, $supermatrix, $taxafile, $random_start, $treefile ) = @_;

    my $ts = Bio::SUPERSMART::Service::TreeService->new;     
	my $log = $self->logger; 
		
	# starting tree is either a given one, an NCBI classification tree, 
	# or a random one (if starttree is empty argument to 'make_usertree')
	my $starttree;
	if ( $treefile ) {
		$log->info("Starting inference from tree $treefile");
		($starttree) = $ts->read_tree( '-file' => $treefile );
		$starttree = $ts->remap_to_ti( $starttree );
	}
	else {
		if ( ! $random_start ) {
			$log->info("Starting inference from NCBI classification tree");
			my $mt = Bio::SUPERSMART::Domain::MarkersAndTaxa->new;
			my @taxatable = $mt->parse_taxa_file( $taxafile );
			$starttree = $ts->make_ncbi_tree( @taxatable );
		} 
		else {
			$log->info("Starting inference from random tree");
		}
	}
	my $usertree = $ts->make_usertree( $supermatrix, $starttree, $self->workdir.'/user.dnd' ); 
	$service->usertree( $usertree );
	return $service;
}

1;
