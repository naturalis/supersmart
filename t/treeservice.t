#!/usr/bin/perl
use strict;
use warnings;
use FindBin '$Bin';
use Test::More 'no_plan';
use Bio::Phylo::IO 'parse_tree';
use Bio::SUPERSMART::Domain::MarkersAndTaxa 'parse_taxa_file';
use Bio::SUPERSMART::Config;

# load the package
BEGIN { use_ok('Bio::SUPERSMART::Service::TreeService'); }

# create a new instance
my $ts = new_ok('Bio::SUPERSMART::Service::TreeService');
my $conf = Bio::SUPERSMART::Config->new;


# test reading trees

# newick format, explicitly specified and autodetect
my $tf = $Bin . '/testdata/testtree.dnd';
(my $tree) = $ts->read_tree( '-file' => $tf, '-format' => 'newick' );
isa_ok ( $tree, 'Bio::Phylo::Forest::Tree');
($tree) = $ts->read_tree( '-file' => $tf );
isa_ok ( $tree, 'Bio::Phylo::Forest::Tree');

# nexus (pure) format, explicitly specified and autodetect
$tf = $Bin . '/testdata/nexus.nex';
($tree) = $ts->read_tree( '-file' => $tf, '-format' => 'nexus' );
isa_ok ( $tree, 'Bio::Phylo::Forest::Tree');
($tree) = $ts->read_tree( '-file' => $tf );
isa_ok ( $tree, 'Bio::Phylo::Forest::Tree');

# figtree nexus format, explicitly specified and autodetect
$tf = $Bin . '/testdata/consensus.nex';
($tree) = $ts->read_tree( '-file' => $tf, '-format' => 'figtree' );
isa_ok ( $tree, 'Bio::Phylo::Forest::Tree');
($tree) = $ts->read_tree( '-file' => $tf );
isa_ok ( $tree, 'Bio::Phylo::Forest::Tree');

# parse multiple trees to array
$tf = $Bin . '/testdata/piperaceae-backbone.dnd';
my @trees = $ts->read_tree( '-file' => $tf, '-format' => 'newick' );
for ( @trees ) {
	isa_ok( $_,  'Bio::Phylo::Forest::Tree');
}

# parse taxon mapping
my $taxa = $Bin . '/testdata/species.tsv';
my $mt = 'Bio::SUPERSMART::Domain::MarkersAndTaxa';
my @records = $mt->parse_taxa_file($taxa);

$tf = $Bin . '/testdata/tree-not-rerooted.dnd';
($tree) = $ts->read_tree( '-file'   => $tf, '-format' => 'newick');

# test rerooting a tree
my $rerooted = $ts->reroot_tree($tree, \@records, ["suborder"]);
isa_ok ($rerooted, 'Bio::Phylo::Forest::Tree');

# rerooted tree should be different
ok(! ($tree->to_newick eq $rerooted->to_newick), "rerooted tree differs from initial one");

# parse phylip files for tree tips
my $phylip_small = "$Bin/testdata/testmatrix-small.phy";
my @names_small = $ts->read_tipnames($phylip_small);
is(scalar @names_small, 5, "read_tipnames from phylip file");

my $phylip_large = "$Bin/testdata/testmatrix.phy";
my @names_large = $ts->read_tipnames($phylip_large);
is($names_large[0], "Orthogeomys_heterodus", "read_tipnames from phylip file");

# build a consensus clade tree
my $treefile  = "$Bin/testdata/testclade.nex";
my $cladetree = $ts->consense_trees( '-infile' => $treefile, '-heights' => $conf->NODE_HEIGHTS );
isa_ok( $cladetree,'Bio::Phylo::Forest::Tree' );

# parse the backbone tree
my $bbfile = "$Bin/testdata/backbone-consensus.nex";
(my $bbtree) = $ts->read_figtree( '-file' => $bbfile );

# graft clade onto backbone tree
my $nof_terminals = scalar(@{$bbtree->get_terminals});
$bbtree = $ts->remap_to_ti($bbtree);
my $grafted = $ts->graft_tree( $bbtree, $cladetree );

isa_ok($grafted,'Bio::Phylo::Forest::Tree');
ok (scalar(@{$grafted->get_terminals}) > $nof_terminals, "tree has more terminals after grafting");

# test remapping a tree back to taxon identifiers 
my $newick = "((('Echinochloa_crus-galli',Echinochloa_colona),Echinochloa_stagnina),(Panicum_turgidum, \"x_Brassolaeliocattleya_'Sung_Ya_Green'\"));";
($tree) = $ts->read_tree( '-string' => $newick, '-format' => 'newick' );
my $remapped = $ts->remap_to_ti($tree);
isa_ok ($remapped, 'Bio::Phylo::Forest::Tree');
