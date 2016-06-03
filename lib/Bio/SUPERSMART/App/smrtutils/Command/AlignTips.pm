package Bio::SUPERSMART::App::smrtutils::Command::AlignTips;

use strict;
use warnings;

use Bio::Phylo::IO qw(parse_tree);
use Bio::SUPERSMART::Service::TreeService; 

use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrtutils qw(-command);

# ABSTRACT: reorder tips to match between two phylogenetic trees

=head1 NAME

AlignTips.pm - reorder tiplabels to match between two trees

=head1 SYNOPSYS

smrt-utils aligntips [-t <file>] [-u <file>] [-f <format>] [-o <file>] [-q <file>] [-h] [-v] [-w <dir>] [-l <file>] [-y] 

=head1 DESCRIPTION

Given two phylogenetic trees, prunes reorders their tips by node order and thus aligns the tips of both trees.

=cut

sub options {    
	my ($self, $opt, $args) = @_;
	my $format_default = 'figtree';
	return (
		['tree1|t=s', "tree file", { arg => 'file' }],		
		['tree2|u=s', "tree file to compare to", { arg => 'file' }],				
		['format|f=s', "file format of output trees, defaults to $format_default", { default => $format_default, arg => "format" }],
		['outtree1|o=s', "outfile 1: filename to export tree1 with reordered tip labels, defaults to <tree1>-sorted", { arg => 'file' }],
		['outtree2|q=s', "outfile 2: filename to export tree2 with reordered tip labels, defaults to <tree2>-sorted", { arg => 'file' }],		
		['prune|p', "prune tips such that both trees have the same taxa" ],
		);	
}

sub validate {
	my ($self, $opt, $args) = @_;			
	$self->usage_error('need two tree files as argument') if not ($opt->tree1 and $opt->tree2);
	$self->usage_error('tree file(s) not found') if not (-e $opt->tree1 and -e $opt->tree2);
	$self->usage_error('tree file(s) empty') if not (-s $opt->tree1 and -s $opt->tree2);	
    if ( $opt->format !~ /^(?:newick|nexus|figtree)$/i ) {
        $self->usage_error("only newick and nexus format are supported");
    }
}

sub run {
	my ($self, $opt, $args) = @_;    

	my $logger = $self->logger;      	
	my $ts = Bio::SUPERSMART::Service::TreeService->new; 	
   
	(my $tree1) = $ts->read_tree( '-file'   => $opt->tree1 );
	(my $tree2) = $ts->read_tree( '-file'   => $opt->tree2 );
  
	if ( $opt->prune ) {
		# prune tips such that both trees have the same taxa.
		$logger->info("Pruning tips to have the same taxa in both trees");
		my @trees = $ts->intersect_trees($tree1, $tree2);
		$tree1 = $trees[0];
		$tree2 = $trees[1];
	}
	
	# get tip labels, sorted by distance to root
	# collect labels of the larger tree
	my $tree = $tree1;
	$tree = $tree2 if scalar( @{$tree2->get_terminals} ) > scalar( @{$tree1->get_terminals} );
	
	my @labels;
	$tree->visit_depth_first(
        # pre-order
        '-pre' => sub {
			my $n = shift;
			if ( $n->is_terminal ) {
				my $name = $n->get_name;
				push @labels, $name;
			}
        }
        );
	
	$tree1->sort_tips(\@labels);
	$tree2->sort_tips(\@labels);
	
	my ($of1, $of2) = $self->_get_outfilenames($opt);

	$ts->to_file( '-tree'=>$tree1, '-format'=> $opt->format, '-file'=>$of1 );
	$ts->to_file( '-tree'=>$tree2, '-format'=> $opt->format, '-file'=>$of2 );
	
	$logger->info("DONE, trees written to $of1, $of2");
	
}

# construct names for output files, if not given.
# output file names default to <infilename>-sorted
sub _get_outfilenames {
	my ( $self, $opt ) = @_;

	my $out1 = $opt->outtree1;
	my $out2 = $opt->outtree2;
	my $suffix = $opt->format eq 'newick' ? '.dnd' : '.nex';
	
	if ( ! $out1 ) {
		($out1 = $opt->tree1) =~ s/\.[^.]+$//;
		$out1 .= '-sorted' . $suffix; 
	}
	if ( ! $out2 ) {
		($out2 = $opt->tree2) =~ s/\.[^.]+$//;
		$out2 .= '-sorted' . $suffix;
	}

	return ($out1, $out2);
}

1;
