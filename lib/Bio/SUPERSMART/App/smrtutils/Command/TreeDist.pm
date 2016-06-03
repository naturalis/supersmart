package Bio::SUPERSMART::App::smrtutils::Command::TreeDist;

use strict;
use warnings;

use Bio::Phylo::IO qw(parse_tree);
use Bio::SUPERSMART::Service::TreeService; 

use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrtutils qw(-command);

# ABSTRACT: calculates the distance between two phylogenetic trees

=head1 NAME

Treedist.pm - Calculates the distance between two phylogenetic trees

=head1 SYNOPSYS

smrt-utils treedist [-t <file>] [-u <file>] [-f <format>] [-d] [-h] [-v] [-w <dir>] [-l <file>] [-y] 

=head1 DESCRIPTION

Calculates the distance between two phylogenetic trees. Currently supported distance measures: Robinson Foulds symmetric distance,
squared euclidean branch lengths distance (see Kuhner & Felsenstein, 1994).

=cut

sub options {    
	my ($self, $opt, $args) = @_;
	my $format_default = 'newick';
	my $measure_default = 'rf';
	return (
		['tree1|t=s', "tree file", { arg => 'file' }],		
		['tree2|u=s', "tree file to compare to", { arg => 'file' }],		
		['distance_measure|d=s', "distance measure, currntly supported: rf (Robinson Foulds),  eb (squared euclidean branch length distance, Kuhner and Felsenstein), defaults to $measure_default", { default => $measure_default}],
	);	
}

sub validate {
	my ($self, $opt, $args) = @_;			
	$self->usage_error('need two tree files as argument') if not ($opt->tree1 and $opt->tree2);
	$self->usage_error('tree file(s) not found') if not (-e $opt->tree1 and -e $opt->tree2);
	$self->usage_error('tree file(s) empty') if not (-s $opt->tree1 and -s $opt->tree2);
	
	my %measures = ( 'rf' => 1, 'eb' => 1 );
	$self->usage_error('distance measure ' . $opt->distance_measure . ' not supported') if not ($measures{$opt->distance_measure});
}

sub run {
	my ($self, $opt, $args) = @_;    

	my $logger = $self->logger;      	
	my $ts = Bio::SUPERSMART::Service::TreeService->new; 	

	(my $tree1) = $ts->read_tree( '-file'   => $opt->tree1 );
	
	(my $tree2) = $ts->read_tree( '-file'   => $opt->tree2 );

	# prune tips such that both trees have the same taxa.
	$logger->info("Pruning tips to have the same taxa in both trees");
	my @trees = $ts->intersect_trees($tree1, $tree2);
	$tree1 = $trees[0];
	$tree2 = $trees[1];

	# calculate robinson-foulds distance and a normalized robinson foulds-distance
	my $diff;		
	my $ndiff;
	if ( $opt->distance_measure eq 'rf' ) { 
		$diff = $tree1->calc_symdiff($tree2);
		$ndiff = $tree1->calc_symdiff($tree2, 1);
	}
	elsif ( $opt->distance_measure eq 'eb' ) { 
		$diff = $tree1->calc_branch_length_score($tree2);
		$ndiff = $tree1->calc_branch_length_score($tree2, 1);				
	}

	$logger->info("Distance (measure " . $opt->distance_measure .  ") between " . $opt->tree1 . " and " . $opt->tree2 . " : $diff, normalized by splits: $ndiff");	
	$logger->info("DONE");
	
}

1;
