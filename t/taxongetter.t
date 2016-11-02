#!/usr/bin/perl
use strict;
use warnings;
use FindBin '$Bin';
use Test::More 'no_plan';
use Bio::SUPERSMART::Config;
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::Matrices::Matrix;
use Bio::SUPERSMART::Service::SequenceGetter;
use Bio::SUPERSMART::Domain::MarkersAndTaxa;

use Data::Dumper;

# BEWARE: this is a lengthy test to run!

my $config = Bio::SUPERSMART::Config->new;

# instantiate logger
my $log = Bio::Phylo::Util::Logger->new( '-level' => WARN, '-class' => 'main' );

# the first tests: can we use and instantiate the MarkersAndTaxaSelector
BEGIN { use_ok('Bio::SUPERSMART::Service::MarkersAndTaxaSelector'); }
my $mts = new_ok('Bio::SUPERSMART::Service::MarkersAndTaxaSelector');
$log->VERBOSE(
	'-level' => ERROR,
	'-class' => 'Bio::SUPERSMART::Service::MarkersAndTaxaSelector',
);

# get species for a given root taxon
my @root_taxa_names = ("Gentianales", "Rubiaceae");

my $expand_rank = "species";
my @taxa_names = $mts->expand_taxa(\@root_taxa_names, $expand_rank);
cmp_ok (scalar @taxa_names, '>', 7000, "expand for root taxa");

@root_taxa_names = ("Primates");
@taxa_names = $mts->expand_taxa(\@root_taxa_names, $expand_rank);
cmp_ok (scalar @taxa_names, '>', 300, "expand for root taxa");

# 'Otolemur garnettii' has a subspecies, which should be 
#  discarded if lowest level is "Species"
@root_taxa_names = ("Otolemur garnettii");
@taxa_names = $mts->expand_taxa(\@root_taxa_names, $expand_rank);

# should be 'Otolemur garnettii' and NOT 'Otolemur garnettii garnettii'
cmp_ok (scalar @taxa_names, '==', 1, "single taxon");
cmp_ok ( $taxa_names[0], 'eq', "Otolemur garnettii", "test subspecies");


# we are going to read the text files in the results/specieslists dir,
my $file = $Bin . '/testdata/testnames.txt';
$log->info("going to read species names from $file");

# read names from file, clean line breaks
open my $fh, '<', $file or die $!;
my @names = <$fh>;
chomp(@names);
$log->info("read species names from $file");

# this will take some time to do the taxonomic name resolution in the
# database and with webservices
my @nodes = $mts->get_nodes_for_names(@names);
ok( @nodes, "found nodes" );

# try taxon Id for 'Acrocomia crispa'
my ($n) = $mts->get_nodes_for_names('Acrocomia crispa');
my $ti = $n->ti;
ok ( $ti, 'found ti for Acrocomia crispa');

# search for taxon that cannot be found anywhere
my @no_nodes = $mts->get_nodes_for_names("xy");
is ( scalar @no_nodes, 0, "nonsense query yields no taxon");

# test getting the nodes for a whole taxon table 
my $mt =  Bio::SUPERSMART::Domain::MarkersAndTaxa->new;
my $taxa = $Bin . '/testdata/species.tsv';
my @records = $mt->parse_taxa_file($taxa);
@nodes = $mts->get_nodes_for_table(@records);
ok( @nodes, "getting nodes for taxa table" );

# test getting nodes with a different taxon table
$taxa = $Bin . '/testdata/species-allranks.tsv';
@records = $mt->parse_taxa_file($taxa);
@nodes = $mts->get_nodes_for_table(@records);
ok( @nodes, "getting nodes for taxa table (many ranks)" );

# test function get_rank_for_taxon
my $tid = 3656;
my $rank = $mts->get_rank_for_taxon($tid);
is ($rank, "species", "taxon $tid is of rank species");

$tid = 2759; 
$rank = $mts->get_rank_for_taxon($tid);
is ($rank, "superkingdom", "taxon $tid is of rank superkingdom");

# test function get_outgroup_taxa
my @ingroup = (319220, 869827, 396991, 442820);
my $tf = $Bin . '/testdata/commontree-cucumis.dnd';
my $tree = parse_tree(
	'-file'   => $tf,
	'-format' => 'newick',
	'-as_project' => 1,
    );
my @outgroup = $mts->get_outgroup_taxa($tree, \@ingroup);

cmp_ok ( scalar @outgroup, ">", 0, "found outgroup species for list of taxon ids");

