# this is an object oriented perl module

package Bio::Phylo::PhyLoTA::DAO::Result::Node;

# Created by DBIx::Class::Schema::Loader
# DO NOT MODIFY THE FIRST PART OF THIS FILE

=head1 NAME

Bio::Phylo::PhyLoTA::DAO::Result::Node

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<nodes>

=cut

__PACKAGE__->table("nodes");

=head1 ACCESSORS

=head2 ti

  data_type: 'int'
  is_nullable: 0
  size: 10

=head2 ti_anc

  data_type: 'int'
  default_value: null
  is_nullable: 1
  size: 10

=head2 terminal_flag

  data_type: 'tinyint'
  default_value: null
  is_nullable: 1
  size: 1

=head2 rank_flag

  data_type: 'tinyint'
  default_value: null
  is_nullable: 1
  size: 1

=head2 model

  data_type: 'tinyint'
  default_value: null
  is_nullable: 1
  size: 1

=head2 taxon_name

  data_type: 'varchar'
  default_value: null
  is_nullable: 1
  size: 128

=head2 common_name

  data_type: 'varchar'
  default_value: null
  is_nullable: 1
  size: 128

=head2 rank

  data_type: 'varchar'
  default_value: null
  is_nullable: 1
  size: 64

=head2 n_gi_node

  data_type: 'int'
  default_value: null
  is_nullable: 1
  size: 10

=head2 n_gi_sub_nonmodel

  data_type: 'int'
  default_value: null
  is_nullable: 1
  size: 10

=head2 n_gi_sub_model

  data_type: 'int'
  default_value: null
  is_nullable: 1
  size: 10

=head2 n_clust_node

  data_type: 'int'
  default_value: null
  is_nullable: 1
  size: 10

=head2 n_clust_sub

  data_type: 'int'
  default_value: null
  is_nullable: 1
  size: 10

=head2 n_piclust_sub

  data_type: 'int'
  default_value: null
  is_nullable: 1
  size: 10

=head2 n_sp_desc

  data_type: 'int'
  default_value: null
  is_nullable: 1
  size: 10

=head2 n_sp_model

  data_type: 'int'
  default_value: null
  is_nullable: 1
  size: 10

=head2 n_leaf_desc

  data_type: 'int'
  default_value: null
  is_nullable: 1
  size: 10

=head2 n_otu_desc

  data_type: 'int'
  default_value: null
  is_nullable: 1
  size: 10

=head2 ti_genus

  data_type: 'int'
  default_value: null
  is_nullable: 1
  size: 10

=head2 n_genera

  data_type: 'int'
  default_value: null
  is_nullable: 1
  size: 10

=cut

__PACKAGE__->add_columns(
  "ti",
  { data_type => "int", is_nullable => 0, size => 10 },
  "ti_anc",
  { data_type => "int", default_value => "null", is_nullable => 1, size => 10 },
  "terminal_flag",
  {
    data_type => "tinyint",
    default_value => "null",
    is_nullable => 1,
    size => 1,
  },
  "rank_flag",
  {
    data_type => "tinyint",
    default_value => "null",
    is_nullable => 1,
    size => 1,
  },
  "model",
  {
    data_type => "tinyint",
    default_value => "null",
    is_nullable => 1,
    size => 1,
  },
  "taxon_name",
  {
    data_type => "varchar",
    default_value => "null",
    is_nullable => 1,
    size => 128,
  },
  "common_name",
  {
    data_type => "varchar",
    default_value => "null",
    is_nullable => 1,
    size => 128,
  },
  "rank",
  {
    data_type => "varchar",
    default_value => "null",
    is_nullable => 1,
    size => 64,
  },
  "n_gi_node",
  { data_type => "int", default_value => "null", is_nullable => 1, size => 10 },
  "n_gi_sub_nonmodel",
  { data_type => "int", default_value => "null", is_nullable => 1, size => 10 },
  "n_gi_sub_model",
  { data_type => "int", default_value => "null", is_nullable => 1, size => 10 },
  "n_clust_node",
  { data_type => "int", default_value => "null", is_nullable => 1, size => 10 },
  "n_clust_sub",
  { data_type => "int", default_value => "null", is_nullable => 1, size => 10 },
  "n_piclust_sub",
  { data_type => "int", default_value => "null", is_nullable => 1, size => 10 },
  "n_sp_desc",
  { data_type => "int", default_value => "null", is_nullable => 1, size => 10 },
  "n_sp_model",
  { data_type => "int", default_value => "null", is_nullable => 1, size => 10 },
  "n_leaf_desc",
  { data_type => "int", default_value => "null", is_nullable => 1, size => 10 },
  "n_otu_desc",
  { data_type => "int", default_value => "null", is_nullable => 1, size => 10 },
  "ti_genus",
  { data_type => "int", default_value => "null", is_nullable => 1, size => 10 },
  "n_genera",
  { data_type => "int", default_value => "null", is_nullable => 1, size => 10 },
);

=head1 PRIMARY KEY

=over 4

=item * L</ti>

=back

=cut

__PACKAGE__->set_primary_key("ti");


# Created by DBIx::Class::Schema::Loader v0.07046 @ 2016-10-31 13:06:27
# DO NOT MODIFY THIS OR ANYTHING ABOVE! md5sum:zX7fUaRr4REgoMAjT0YvmA

# You can replace this text with custom code or comments, and it will be preserved on regeneration

use Bio::SUPERSMART::Config;
use Bio::Phylo::Forest::NodeRole;
use Bio::SUPERSMART::Service::MarkersAndTaxaSelector;
push @Bio::Phylo::PhyLoTA::DAO::Result::Node::ISA, 'Bio::Phylo::Forest::NodeRole';

my $mts = Bio::SUPERSMART::Service::MarkersAndTaxaSelector->new;
my %tree;

=head2 table

Getter/setter that maps this ORM class onto the correct version (184) of the underlying
database table.

=cut

sub table {
	my $class = shift;
	my $table = shift;
	my $release = Bio::SUPERSMART::Config->new->currentGBRelease;
	$class->SUPER::table( $table . '_' . $release );
}

=head2 get_parent

Returns direct parent node.

=cut

sub get_parent {
	my $self = shift;
	if ( $self->get_generic('root') ) {
		return;
	}
	if ( my $parent_ti = $self->ti_anc ) {
		return $mts->find_node($parent_ti);
	}
	return;
}

=head2 set_parent

This is a no-op: the tree structure is immutable.

=cut

sub set_parent { return shift }

=head2 get_children

Returns array ref of direct children.

=cut

sub get_children {
	my $self = shift;
	my $ti = $self->ti;
	my @children = $mts->search_node( { ti_anc => $ti } )->all;
	return \@children;
}

=head2 get_siblings

Returns array ref of siblings

=cut

sub get_siblings {
	my $self = shift;

	return [ grep { $_->ti != $self->ti } @{ $self->get_parent->get_children } ];
}

=head2 get_descendants_at_rank

Given a taxonomic rank, gets all the decendants of the
invocant which are of that rank.

=cut

sub get_descendants_at_rank {
	my ( $self, $rank ) = @_;
	
	my @result;
	my @queue = ($self);
	while (@queue) {
		my $current = shift (@queue);
		if ( $current->rank eq $rank ) {
			push @result, $current;
		}
		if ( my @ch = @{$current->get_children} ) {
			push @queue, @ch;
		} 
	}
	return \@result;
} 

=head2 get_branch_length 

Returns nothing: in this implementation (i.e. a taxonomy) there are no branch lengths

=cut

sub get_branch_length { return }

=head2 set_branch_length

This is a no-op: the tree structure is immutable.

=cut

sub set_branch_length { return shift }

=head2 get_id

Alias for C<ti>.

=cut

sub get_id { shift->ti }

=head2 set_tree

Stores a reference to the containing tree, if any.

=cut

sub set_tree {
	my ( $self, $tree ) = @_;
	$tree{ $self->get_id } = $tree;
	return $self;
}

=head2 get_tree

Returns a reference to the containing tree, if any.

=cut

sub get_tree { $tree{ shift->get_id } }

=head2 get_name

Alias for taxon_name

=cut

sub get_name { shift->taxon_name }

=head2 get_genus

Alias for ti_genus

=cut

sub get_genus {shift->ti_genus}

1;
