package Bio::SUPERSMART::App::smrtutils::Command::Replicate;

use strict;
use warnings;

use File::Spec;
use Data::Dumper;
use List::Util qw(max);

use Bio::Phylo::IO qw(parse unparse parse_matrix);
use Bio::Phylo::Util::CONSTANT ':objecttypes';
use Bio::Phylo::Models::Substitution::Dna;

use Bio::SUPERSMART::Service::MarkersAndTaxaSelector;
use Bio::SUPERSMART::Domain::MarkersAndTaxa;
use Bio::SUPERSMART::Config;
use Bio::SUPERSMART::Service::TreeService;
use Bio::SUPERSMART::Service::ParallelService 'pfm';

use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrtutils qw(-command);

# ABSTRACT: replicate trees and alignments with simulations


=head1 NAME

Replicate - make simulated replicates of trees and alignments

=head1 SYNOPSYS

This subcommand provides functionality to generate a complete synthetic dataset from a supersmart run.
The generated synthetic dataset consists of a replicated (final) species tree and replicated alignments,
which resemble the input alignments.
The replicated dataset can then be used as input for the SUPERSMART pipeline to validate the method and
parameters chosen for the initial run.

=head1 DESCRIPTION

=cut

my $config = Bio::SUPERSMART::Config->new;

sub options {
	my ($self, $opt, $args) = @_;
	my $aln_outfile_default = 'aligned-replicated.txt';
	my $tree_outfile_default = 'tree-replicated.dnd';
	my $taxa_outfile_default = 'taxa-replicated.tsv';
	return (
		['tree|t=s', 'file name of input tree, must be ultrametric', { arg => 'file', mandatory => 1 } ],
		['alignments|a=s', "list of alignment files to replicate, as produced by 'smrt align'", { arg => 'file' } ],
		["aln_outfile|o=s", "name of output file listing the simulated alignments, defaults to '$aln_outfile_default'", {default => $aln_outfile_default, arg => "file"}],
		["tree_outfile|b=s", "name of the output tree file (newick format), defaults to '$tree_outfile_default'", {default => $tree_outfile_default, arg => "file"}],
		['replicated_tree|r=s', 'Replicated tree from former replication', { arg => 'file' } ],
		["taxa_outfile|c=s", "name the output taxa file", {default => $taxa_outfile_default, arg => "file"}],
		["no_marker_sim|n", "do not simulate marker presence with binary simulation, take taxa from original alignments instead", { default=> 0}],
		["random_rootseq|d", "simulate alignments from random sequence; otherwise the median of the sequences will be the root sequence", { default=> 0}],
		["ids|i", "return NCBI identifiers in remapped tree instead of taxon names", { default=> 0}],
		["phylogenetically_informative|p", "Discard simulated alignments in which > 50% sequences are identical, and simulate again  ", { default=> 0}],
	    );
}

sub validate {
	my ($self, $opt, $args) = @_;

	# We only have to check the 'infile' argument.
	#  If the infile is absent or empty, abort
	my $file = $opt->tree;
	$self->usage_error('no tree argument given') if not $file;
	$self->usage_error("file $file does not exist") unless (-e $file);
	$self->usage_error("file $file is empty") unless (-s $file);
}

sub run {
	my ($self, $opt, $args) = @_;

	my $logger = $self->logger;
	my $treefile = $opt->tree;
	my $aln_outfile = $opt->aln_outfile;
	my $tree_outfile = $opt->tree_outfile;
	my $taxa_outfile = $opt->taxa_outfile;   
	my $ts = Bio::SUPERSMART::Service::TreeService->new;
	my $mt  = Bio::SUPERSMART::Domain::MarkersAndTaxa->new;
	my $mts    = Bio::SUPERSMART::Service::MarkersAndTaxaSelector->new;

	# read tree
	(my $tree) = $ts->read_tree( '-file' => $treefile );
	
	# prune negative branches from tree, if any
	$tree = $self->_prune_negative_branches($tree);

	# replicate tree and write to file, if not replicated already
	my $tree_replicated;
	if ( my $filename_rep = $opt->replicated_tree ) {
		#load replicate tree from file
		$logger->info("Reading replicated tree from $filename_rep");
		($tree_replicated) = $ts->read_tree( '-file'   => $filename_rep );
	}
	else {
		$logger->info("Replicating tree from $treefile");
		$tree_replicated = $self->_replicate_tree($tree)->first;
		$ts->to_file( '-tree' => $tree_replicated, '-file' => $tree_outfile );
		$logger->info("wrote replicated tree to $tree_outfile");
	}

	$self->_write_taxafile($tree_replicated, $taxa_outfile) if not -e $taxa_outfile;
	$logger->info("wrote taxa file to $taxa_outfile");
	my @records = $mt->parse_taxa_file( $taxa_outfile );

	if ( my $aln = $opt->alignments ) {
		
		$ts->remap_to_ti( $tree, @records );
		$ts->remap_to_ti( $tree_replicated, @records );

		# read list of alignments
		my @alnfiles = $mt->parse_aln_file( $aln );

		# replicate all alignments fiven in input align file,
		# write alignments to file and also create a list with all
		# newly written alignments
	    open my $outfh, '>', $aln_outfile or die $!;

		my @replicated = pmap {
			my ($filename_orig) = @_;

			$logger->info("Attempting to replicate alignment $filename_orig");

			# set random seed to prevent issue with forking and chosing tempfile names
			my $seed = 0;
			$seed += $_ for  map {ord $_} split(//, $filename_orig);
			srand($seed);

			# file name for replicated alignment
			my ( $volume, $directories, $filename_rep ) = File::Spec->splitpath( $filename_orig );
			$filename_rep =~ s/\.fa$/-replicated\.fa/g;
			$filename_rep = $self->workdir . '/' . $filename_rep;
			
			# replicate the alignment
			my $rep_aln = $self->_replicate_alignment( 'alnfile' => $filename_orig, 
													   'tree_replicated' => $tree_replicated, 
													   'tree_original' => $tree, 
													   'simulate_markers' => ! $opt->no_marker_sim, 
													   'random_rootseq' => $opt->random_rootseq,
				                                       'informative' => $opt->phylogenetically_informative);
			
			if ( $rep_aln ) {
				# simulated alignment will have the same file name plus added '-simulated'
				$logger->info("Writing alignment $filename_rep");
				unparse ( -phylo => $rep_aln, -file => $filename_rep, -format=>'fasta' );
				# output average distances in alignment
				my $dist_orig = $self->_mean_dist($filename_orig);
				my $dist_rep = $self->_mean_dist($filename_rep);
				my $ident_orig = $self->_num_identical_seqs($filename_orig);
				my $ident_rep = $self->_num_identical_seqs($filename_rep);
				$logger->info("Average distance in alignments: original : $dist_orig, replicated : $dist_rep");
				$logger->info("Pairs of identical seqs in alignments: original : $ident_orig, replicated : $ident_rep");
				
			}
			else {
				$logger->warn("Could not write replicated alignment to file: alignment not replicated ");
				return 0;
			}
			# write filename to alignment list
			print $outfh "$filename_rep\n";

			return $filename_rep;

		} @alnfiles;

		close $outfh;
		$logger->info("Replicated " . scalar(grep {$_} @replicated) . " of " . scalar(@alnfiles) . " alignments from $aln");
		$logger->info("Alignment written to $aln_outfile");
	}

	$logger->info("DONE. Tree in $tree_outfile, taxa in $taxa_outfile" );
	return 1;
}

# This produces a taxa file givan a tree
sub _write_taxafile {
	my ($self, $tree, $filename) = @_;

	my $logger = $self->logger;
	my $mts    = Bio::SUPERSMART::Service::MarkersAndTaxaSelector->new;
	my $mt  = Bio::SUPERSMART::Domain::MarkersAndTaxa->new;

	my @names = map { s/_/ /g; $_;} map {$_->get_name} @{$tree->get_terminals};

	# identify terminals with artificial names and
	my %artificial = map {$_=>1} grep {! $mts->find_node({taxon_name=>$_})} @names;
	@names = grep {! $artificial{$_}} @names;

	# write initial taxa table to file
	$logger->info("Going to make taxa table");
	my @taxa_table = $mts->make_taxa_table( '-taxon_names' => \@names );
	$logger->info("Writing taxa table to $filename");
	$mts->write_taxa_file( '-file' => $filename, '-table' => \@taxa_table );

	# for the species with artificial taxon names, do as follows: first: assign an artificial
	# ti to the species. Then iteratively select the children of the father of the focal node, and
	# check if any of its children have a classification, until we found a (maybe far away) sister
	# node that is classified. Then copy the genus, family etc. and append the entries to the table
	my $new_ti = $mts->max_ti + 1;

	my @levels = reverse $mts->get_taxonomic_ranks;

	my @records = $mt->parse_taxa_file($filename);

	open my $fh, '>>', $filename or die $!;
	# iterate over all artificial species, assign taxon id and append entry to taxa table
	for my $an ( keys %artificial ) {
		(my $name = $an) =~ s/ /_/g;
		$logger->debug("Looking for closest non-artificial relative of artificial species $name");
		my ($node)  = grep {$_->get_name eq $name} @{$tree->get_terminals};
		my $relative;

		# search for the closest 'real' relative
		while ( ! $relative and ! $node->is_root ) {
			my @relatives;
			$node->visit_depth_first(-in => sub {
				my $n = shift;
				(my $current_name = $n->get_name) =~ s/_/ /g;
				if ( $current_name and (! $artificial{$current_name} )) {
					push @relatives, $current_name;
				}

						 });
			$relative = $relatives[0] if scalar(@relatives);
			# if the whole clade consisted of artificial species, go one up
			$node = $node->get_parent;
		}
		$logger->debug("Closest non-artificial relative: $relative");

		# get all taxon IDs for the closest nonartificial relative, sorted from low to high ranks
		my $tid = $mts->find_node({taxon_name=>$relative})->ti;
		my ($all_ids) = grep { my %h = %{$_}; grep{ $tid } values%h } @records;

		# assign taxon ID to artificial species and take all other taxon ids from closest real relative
		$all_ids->{'species'} = $new_ti++;
		my @ids = @$all_ids{@levels};
		$logger->debug("Assinging taxon id $new_ti to artificial species $name");

		# append information to taxon table
		$name =~ s/_/ /g;
		print $fh $name . "\t" . join("\t", @ids) . "\n";
	}
	close $fh;
	$logger->info("Added " . scalar(keys %artificial) . " artificial species to taxa file $filename") if keys %artificial;
}

sub _replicate_tree {
	my ($self, $tree) = @_;

	$tree->generize(
		'-delim'     => '_',
		'-monotypic' => 1,
		'-polypara'  => 1,
	    );

	my $rep = $tree->replicate( '-genera' => 1);

	return $rep;
}

# Replicate an alignment. Arguments:
#   'alnfile'
#	'tree_replicated' 
#	'tree_original'
#	'simulate_markers'
#	'random_rootseq'

sub _replicate_alignment {
	my ($self, %args ) = @_; #$fastafile, $tree, $original_tree, $simulate_presence) = @_;
	
	my $alnfile = $args{'alnfile'};
	my $tree_replicated = $args{'tree_replicated'};
	my $tree_original = $args{'tree_original'};
	
	my $logger = $self->logger;
	
	# create matrix object from FASTA
	$logger->fatal("Alignment file $alnfile does not exist") if  ( ! ( -e $alnfile and -s $alnfile ) );

	my $project = parse(
		'-format'     => 'fasta',
		'-type'       => 'dna',
		'-file'     => $alnfile,
		'-as_project' => 1,
	    );
	my ($matrix) = @{ $project->get_items(_MATRIX_) };
	$matrix = $self->_clean_fasta_defline( $matrix );

	if ( scalar (@{$matrix->get_entities}) < 3 ) {
		$logger->warn("Cannot replicate alignment $alnfile with less than three sequences. Skipping.");
		return 0;
	}

	# determine for which taxa we want replicated sequences
	my @rep_taxa;

	if ( $args{ 'simulate_markers' } ) {
		# do a binary charachter simulation
		$logger->info("Determining taxa in alignment with binary simulation");
		@rep_taxa = $self->_simulate_marker_presence( '-matrix'=>$matrix, '-tree'=>$tree_original, '-replace'=>0 );
	}
	else {
		# take the taxa that are in the original alignment
		$logger->info("Using taxa from original alignment");
		@rep_taxa = map { $_->get_name } @{$matrix->get_entities};
	}
	$logger->debug("Taxa to replicate sequences for : " . Dumper(\@rep_taxa) );

	if ( scalar(@rep_taxa) < 3 ) {
		$logger->warn( "Less than three taxa predicted to have marker. Skipping replication of $alnfile." );
		return 0;
	}

	# determine substitution model for given alignment
	my $timeout = 1800; # set to 30 min
	my $model = 'Bio::Phylo::Models::Substitution::Dna'->modeltest( '-matrix' => $matrix, '-timeout' => $timeout, '-tree' => $tree_original );

	# timeout in model test will most surely also timeout in simulation; therefore do not simulate
	if ( ! $model ) {
		$logger->warn("No substitution model determined. Skipping.");
		return 0;
	}

	# prune tree for faster sequence simulations
	my $pruned = parse('-format'=>'newick', '-string'=>$tree_replicated->to_newick)->first;
	$pruned->keep_tips( \@rep_taxa );
	$logger->debug("Input tree for simulation : " . $pruned->to_newick);

	# simulate sequences
	#  In case the 'informative' option is set: if the number of identical sequences in the alignments is greater than 50%, we  
	#  try 10 times to get a phylogenetically informative alignment by re-simulating.   
	my $rep;
	my $frac_ident = 1;
	my $tries = 0;
	while ( $frac_ident > 0.5) {
		$tries++;
		$logger->info("Try # $tries to simulate phylogenetically informative alignment.") if $args{'informative'};
		if ( $tries >= 10 ) {
			$logger->warn("Exceeded 10 tries for simulating phylogenetically informative alignment. Skipping.");
			return 0;
		}
		
		eval { $rep = $matrix->replicate('-tree'=>$pruned, '-seed'=>$config->RANDOM_SEED, '-model'=>$model, '-random_rootseq'=> $args{'random_rootseq'} ); };
		if ( $@ ) {
			$logger->warn("Problem replicating matrix, " . $@);
			return 0;
		}
		
		my %rp = map {$_=>1} @rep_taxa;
		for my $r ( @{$rep->get_entities} ) {
			if ( ! $rp{$r->get_name} ) {
				$logger->info("Removing taxon " . $r->get_name . " from replicated matrix");
				$rep->delete($r);
			}
		}
		$frac_ident = $self->_get_frac_identical_seqs( $rep );
		# never try to re-simulate if '--phylogenetically_informative' flag not set
		last if not $args{'informative'};
	}
	
	$logger->info("Number of seqs in original alignment: " . scalar(@{$matrix->get_entities}) . ", number of seqs in rep alignment: " . scalar(@{$rep->get_entities}));
	# If we had less than two simulated marker presences, the replicated alignment is not an alignment, therefore skip
	if ( @{ $rep->get_entities } < 2 ) {
		$logger->warn("Replication produced alignment with less than 2 sequences, skipping");
		return 0;
	}
	$rep = $self->_fix_fasta_defline( $rep );

	return $rep;
}

# given an alignment as matrix object, returns the highest fraction
#  of sequences that are identical in the alignment 
sub _get_frac_identical_seqs {
	my ($self, $matrix) = @_;	

	my @rows = @{$matrix->get_entities};
	my %ident_by_seq;
	$ident_by_seq{$_->get_char}++ for @rows;
	my $frac = max (values (%ident_by_seq)) / scalar(@rows);
	
	return $frac;
}

# Given an alignment as matrix object and a tree, simulates the presence in the alignment
# for each taxon in the tree using a birth-death simulation on a binary character
sub _simulate_marker_presence {
	my ($self,%args) = @_;

	my $matrix = $args{'-matrix'};
	my $tree = $args{'-tree'};

	my $logger = $self->logger;

	# collect all taxa from alignment that are also present in the tree
	my %taxa_in_aln = map{ $_->get_name=>1 }  @{ $matrix->get_entities };
	my %taxa_in_tree = map { $_->get_name=>1 } @{ $tree->get_terminals };
	my %marker_taxa = map {$_=>1} grep { exists $taxa_in_tree {$_} } keys(%taxa_in_aln);

	# if desired, taxa not present in the tree are substituted by random taxa that are in the tree
	if ( $args{'-replace'} ) {
		if ( scalar(keys %marker_taxa) < scalar(keys %taxa_in_aln) ) {
			$logger->info("Some taxa in alignment are not in the tree. Adding random taxa from tree to taxon set.");
			%marker_taxa =  $self->_add_random_taxa( \%marker_taxa, $tree, scalar(keys %taxa_in_aln) );
		}
	}

	# make binary matrix with occurrences of taxa in the alignment
	my $bin = [];
	for my $tax ( keys %taxa_in_tree ) {
		my $present = $marker_taxa{$tax} ? '1':'0';
		push @$bin, [$tax=>$present];
	}
	my $fac = Bio::Phylo::Factory->new;
	my $marker_matrix = $fac->create_matrix( '-matrix' => $bin);

	# run binary simulation to determine taxa that are simulated to be present
	my $rep_marker_matrix = $marker_matrix->replicate('-tree'=>$tree, '-seed'=>$config->RANDOM_SEED);

	# get taxa from replicate that have a simulated marker presence
	my %simulated_taxa =  map {$_->[0]=>1} grep {$_->[1] == 1} @{$rep_marker_matrix->get_raw};

	return ( keys %simulated_taxa );
}

# given a set of taxon id's and a tree, randomly add taxa from the tree to the set until it
# reaches user supplied size
sub _add_random_taxa {
	my ($self, $set, $tree, $size) = @_;

	my %taxa = %{$set};
	my @tree_taxa = map { $_->get_name } @{ $tree->get_terminals };

	while ( scalar(@tree_taxa) > scalar(keys %taxa) and scalar(keys %taxa) < $size ) {
		my $id = @tree_taxa[rand @tree_taxa];
		$self->logger->debug("attempting to add taxon $id to set of taxa");
		$taxa{$id} = 1;
	}

	return %taxa;
}

# removes excess '>' from fasta def line
sub _fix_fasta_defline {
	my ($self, $mat) = @_;
	for my $seq ( @{ $mat->get_entities }) {
		# need to fix definition line to prevent a ">>" in FASTA definition line (issue #21 in Bio::Phylo)
		my $defline = $seq->get_generic('fasta_def_line');
		$defline =~ s/>//g;
		$seq->set_generic('fasta_def_line', $defline);
	}
	return $mat;
}

# Change the FASTA definition line from the SUPERSMART style,
# e.g. ">gi|443268840|seed_gi|339036134|taxon|9534|mrca|314294/1-1140"
# to contain only the taxon ID, e.g. 9534
sub _clean_fasta_defline {
	my ($self, $matrix) = @_;

	for my $seq ( @{ $matrix->get_entities } ) {
		my $seqname = $seq->get_name;
		if ( $seqname =~ m/taxon\|([0-9]+)/ ) {
			$seqname = $1;
			$self->logger->debug('Changing FASTA definition line from ' . $seq->get_name . " to $seqname");
			$seq->set_name($seqname);
		}
		$seq->set_generic('fasta_def_line'=>$seqname);
	}

	return $matrix;
}

sub _mean_dist {
	my ($self, $filename) = @_;
	my $mt    = Bio::SUPERSMART::Domain::MarkersAndTaxa->new;
	open my $fh, '<', $filename;
	read $fh, my $string, -s $fh;
	close $fh;
	my $dist = $mt->calc_mean_distance($string);
	return $dist;
}

sub _num_identical_seqs {
	my ($self, $file) = @_;

	$self->logger->info("Checking for identical sequences in alignment $file");
	my $mts    = Bio::SUPERSMART::Service::MarkersAndTaxaSelector->new;

	my $matrix = parse_matrix(
		'-file'       => $file,
		'-format'     => 'fasta',
		'-type'       => 'dna',
		);

	my @rows = @{$matrix->get_entities};

	my $identical = 0;
	for (my $i=0; $i<=$#rows; $i++) {
		my $row_i = $rows[$i];
		for (my $j=$i+1; $j<=$#rows; $j++) {
			my $row_j = $rows[$j];

			if ($row_i->get_char eq $row_j->get_char) {

				my $name_i = $row_i->get_name;
				my $name_j = $row_j->get_name;
				if ( $name_i=~/.*taxon\|([0-9]+)/) {
					$name_i =$1;
				}
				if ( $name_j=~/.*taxon\|([0-9]+)/) {
					$name_j =$1;
				}

				$name_i = $mts->find_node($name_i)->taxon_name;
				$name_j = $mts->find_node($name_j)->taxon_name;

				$identical++;
				$self->logger->info("seqs " . $name_i . " and " . $name_j  . " are identical !");
			}
		}
	}
	return $identical;
}

sub _prune_negative_branches {
	my ($self, $tree) = @_;

	my $logger = $self->logger;
	my $count = 0;
	$tree->visit(sub{
		my $node = shift;
		my $bl = $node->get_branch_length || 0;
		if ( $bl < 0 ) {
			my @terminals = map {$_->get_name} @{$node->get_terminals};
			$count += scalar @terminals;
			$logger->info("Found negative branch length $bl, following tips will be removed: " . join(', ', @terminals));
			$node->get_parent->prune_child( $node );
		}
				 });
	$tree->remove_unbranched_internals;
	$logger->info("Removed $count terminals pruning all negative branches");
	return $tree;
}

1;
