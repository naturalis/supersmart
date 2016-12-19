package Bio::SUPERSMART::App::smrtutils::Command::TreeInfo;

use strict;
use warnings;

use Bio::Phylo::IO qw(parse parse_tree unparse);
use Bio::SUPERSMART::Service::MarkersAndTaxaSelector;

use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrtutils qw(-command);

# ABSTRACT: information on a phylogenetic tree


=head1 NAME

TreeInfo

=head1 SYNOPSYS

aa

=head1 DESCRIPTION

=cut

sub options {
	my ($self, $opt, $args) = @_;
	my $format_default = 'figtree';
	return (
		['tree|t=s', 'file name of input tree', { arg => 'file', mandatory => 1 } ],
		['tree_format|f=s', "format of tree input file (newick or nexus), defaults to $format_default", { default => $format_default } ],
	    );
}

sub validate {
	my ($self, $opt, $args) = @_;

	# We only have to check the 'infile' argument.
	#  If the infile is absent or empty, abort
	my $file = $opt->tree;
	$self->usage_error('no tree argument given') if not $file;
}

sub run {
	my ($self, $opt, $args) = @_;

	my $logger = $self->logger;
	my $treefile = $opt->tree;

	# read tree
	my $tree = parse_tree(
		'-file'   => $treefile,
		'-format' => $opt->tree_format,
	    );

	# calculate average posterior support value
	my $sum;
	my $cnt;
	for my $node ( @{$tree->get_internals} ) {
		if ( my $posterior = $node->get_meta_object('fig:posterior') || $node->get_meta_object('fig:bootstrap')) {
			$sum += $posterior;
			$cnt++;
		}
	}

	if ( $cnt ) {
		my $avg_post = sprintf("%.3f", $sum/$cnt);
		$logger->info("Average posterior or bootstrap value: $avg_post");
	}

	# get number of terminals
	my $num_term = scalar(@{$tree->get_terminals});

	$logger->info("Number of terminals : $num_term");

	return 1;
}

1;

