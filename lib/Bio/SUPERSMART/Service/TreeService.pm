package Bio::SUPERSMART::Service::TreeService;
use strict;
use warnings;
use Cwd;
use File::Temp 'tempfile';
use Bio::Phylo::IO qw'parse parse_tree unparse';
use Bio::Phylo::Factory;
use Bio::Phylo::Util::CONSTANT ':objecttypes';
use Bio::SUPERSMART::Config;
use Bio::SUPERSMART::Service;
use Bio::Phylo::Util::CONSTANT ':namespaces';
use Bio::Phylo::Util::Exceptions 'throw';
use base 'Bio::SUPERSMART::Service';
use Bio::SUPERSMART::Domain::MarkersAndTaxa;
use Bio::SUPERSMART::Service::MarkersAndTaxaSelector;
use List::MoreUtils 'uniq';
use List::Util qw'min sum';
use Data::Dumper;

my $config = Bio::SUPERSMART::Config->new;
my $log = Bio::SUPERSMART::Service->logger->new;
my $fac = Bio::Phylo::Factory->new;

=head1 NAME

Bio::SUPERSMART::Service::TreeService - Operations on phylogenetic Trees

=head1 DESCRIPTION

Provides functionality for handling trees, e.g. decomposing, rerooting and grafting
clade trees onto a backbone tree.

=over

=item to_file

Writes a tree object to file.

 -file => filename
 -tree => tree object or arrayref of tree objects
 -format => output format, newick, nexus or figtree, defaults to newick

=cut

sub to_file {
	my ( $self, %args ) = @_;
	
	my $file = $args{'-file'} or throw 'BadArgs' => "Need -file argument";
	my $string = $args{'-string'};
	my $tree = $args{'-tree'};
	my $format = $args{'-format'} || 'newick';

	throw 'BadArgs' => "Need -tree or -string argument" if not $tree and not $string;
	$string = $self->to_string( $tree, $format ) if not $string;
   
	open my $fh, '>', $file or die $!;
	print $fh $string;
	close $fh;
	
	$log->info("Wrote tree in $format format to file $file.");
}


=item to_string

Returns the string representation of a tree in the specified format

=cut

sub to_string {
	my ( $self, $tree, $format ) = @_;

	my $fac = Bio::Phylo::Factory->new;
	my $forest = $fac->create_forest;
	
	# either one tree is provided or an arrayref of trees
	if ( ref($tree) ) {
		$forest->insert($_) for @{$tree};
	}
	else {
		$forest->insert($tree);	
	}
	
	# create output
	my $project = $fac->create_project;
	my $taxa = $forest->make_taxa;
	$project->insert($taxa);
	$project->insert($forest);

	my $string = unparse(
		'-format' => $format,
		'-phylo'  => $project,
	    );

	return $string;
}

=item read_tree

Reads phylogenetic tree(s) from file, detects format if not given as argument,
returns an array with trees

 -file => filename
 -string => tree string
 -format => format

=cut 

sub read_tree {
	my ( $self, %args ) = @_;
	
	my $file = $args{'-file'}; 
	my $format = $args{'-format'} || $self->_detect_treeformat( $file );
	my $string = $args{'-string'};
	
	throw 'BadArgs' => "Need -file or -string argument" if not $file and not $string;
	
	my %parse_args = %args;
	$parse_args{'-as_project'} = 1;
	$parse_args{'-format'} = $format;
	my $project = parse(%parse_args);

	my @trees = @{ $project->get_items( _TREE_ ) };
	
	if ( ! scalar( @trees) ) {
		$log->fatal("Could not parse tree(s) in file $file");
	}
	return @trees;
}

sub _detect_treeformat {
	my ($self, $file) = @_;
	
	my $format;   
	$log->info("Detecting format of tree $file");
	open my $fh, '<', $file;
	
	if ( grep {/#NEXUS/} <$fh> ) {
		if ( grep {/FigTree/} <$fh> ) {
			$format = 'figtree';
		}
		else {
			$format = 'nexus';
		}
		$log->info("Detected $format format");
	}
	else {
		$format = 'newick';
		$log->info("File not in nexus or figtree format; assuming newick format");
	}

	return $format;
}
	
=item read_figtree

Reads a tree in figtree-NEXUS format

 -file   => filename with figtree tree
 -string => tree in string representation

=cut

sub read_figtree {
	my ( $self, %args ) = @_;

	my $logger = $self->logger;

	my %parse_args = ( '-format' => 'figtree' );
	if ( my $file    = $args{'-file'} ) {
		# parse tree from fiel
		$logger->debug("going to read $file as figtree/NEXUS");
		$parse_args{ '-file' } = $file;
	}
	elsif ( my $string    = $args{'-string'}) {
		$logger->debug("going to read tree $string");
		$parse_args{ '-string' } = $string;
	}
	else {
		$logger->warn('Could not read figtree tree, need -string or -file argument');
	}

	my $tree = parse_tree( %parse_args );

	return $tree;
}

=item reroot_tree

Reroots a backbone tree. All possible rerootings are evaluated and the tree is returned
that minimizes the paraphyly with respect to the exemplar species of all genera.

=cut

sub reroot_tree {
	my ( $self, $tree, $taxatable, $levels ) = @_;

	my @records = @{$taxatable};

	# taxonomic ranks to consider for determining paraphyly
	my @levels = @{$levels};

	# store the scores (number of paraphyletic species per rerooting) for each level and
	# node index at which we reroot
	my %scores = ();
	@scores{ @levels } = ();

	# store the tree objects for the rerooting at each node
	my %rerooted_trees = ();

	my $num_internals = scalar @{ $tree->get_internals };
	my $newi          = $tree->to_newick(nodelabels=>1);

	# iterate over all possible rerootings
	for my $i (0 .. $num_internals - 1 ) {

		# get fresh tree, so the node order etc. won't be messed up
		my $current_tree = parse(
			'-string' => $newi,
			'-format' => 'newick',
		)->first;

		my @internals = @{ $current_tree->get_internals };

		# reroot the tree
		if ( my $node = $internals[$i] ) {
			$log->debug("rerooting tree at internal node # " . ($i+1) . "/" . $num_internals);

			$node->set_root_below();

			# store the tree for later
			$rerooted_trees{$i} = $current_tree;

			# calculate the 'scores' : count the amount of paraphyly at each taxonomic level
			foreach my $level (@levels) {
				my $para_count =
					$self->_count_paraphyletic_species( $current_tree, $level, @records );
				$log->debug("Paraphyletic species (counted from level $level) : $para_count");
				$scores{$level}{$i} = $para_count;
			}
		}
	}

        # traverse all scores and just keep the node indices for the min scores for each level
	my %best_node_idx = ();
	foreach my $level (@levels) {
		$best_node_idx{$level} = ();
		my $minscore = min( values %{$scores{$level}} );
		foreach my $k ( keys %{$scores{$level}} ) {
			if ( $scores{$level}{$k} == $minscore ) {
				push @{ $best_node_idx{$level} }, $k;
			}
		}
	}

	# now we take the intersection of all indices for all the levels
	my @best_indices =
	  @{ $best_node_idx{ $levels[0] } };   #initialize with set from first level

	my $count = @levels;
	for my $ii ( 1 .. $count - 1 ) {
		@best_indices =
		  _intersect( \@best_indices, $best_node_idx{ $levels[$ii] }  );
	}

	if ( scalar @best_indices > 1 ) {
		$log->warn ("Found more than one optimal rerooted tree. Returning the first one");
	}
	if ( scalar @best_indices == 0 ) {
		$log->warn ("Found no optimal rerooted tree. Returning the original one");
		return $tree;
	}
	$rerooted_trees{ $best_indices[0] }->get_root->set_branch_length(undef);
	return $rerooted_trees{ $best_indices[0] };
}

sub _count_paraphyletic_species {
	my ( $self, $current_tree, $level, @records ) = @_;

	# get all distinct members of the NCBI tree at the given level
	my @members = uniq( map { $_->{$level} } @records );

	# sum up the amount of paraphyletic species for all members at that level
	my $para_count = 0;

	for my $m (@members) {

	   # get the subset of terminals in the NCBI tree that belong to this family
		my @ncbi_terminal_names =
		  map { $_->{'species'} } grep { $_->{$level} eq $m } @records;

		# get node objects for species in our tree
		my @nodes = map{ $current_tree->get_by_name( $_ )  } @ncbi_terminal_names;

		# get most recent common ancestor of the species in this family in our tree
		my $mrca = $current_tree->get_mrca( \@nodes );
		if ($mrca) {

			# get descendants of the mrca
			my @terminals = @{ $mrca->get_terminals };
			my @terminal_names = map { $_->get_name } @terminals;

			# count how many terminals of the mrca in our tree are in the
			#  terminals in the NCBI tree, paraphyletic species are then the
			#  ones that are descendents in our tree but not in the NCBI tree
			my @matching = _intersect( \@ncbi_terminal_names, \@terminal_names );
			my @paraspecies = _array_minus( \@terminal_names, \@ncbi_terminal_names );
			$para_count += scalar @paraspecies;
		}
	}
	return $para_count;
}

=item intersect_trees

	Given two trees, prunes the trees such that both will have the 
	same set of tip labels

=cut

sub intersect_trees {
	my ($self, $tree1, $tree2) = @_;
	
	# get leaves shared by both trees
	my @terminals1 = map { $_->get_name } @{$tree1->get_terminals};
	my @terminals2 = map { $_->get_name } @{$tree2->get_terminals};	
	my @intersect = _intersect(\@terminals1, \@terminals2);

	# prune trees
	$log->debug("Pruning trees to tips: " . Dumper(\@intersect));
	$tree1->keep_tips(\@intersect);
	$tree2->keep_tips(\@intersect);
	
	return ( $tree1, $tree2 );
	
}

# get intersection between elements of two array references
sub _intersect {
	my ($a1, $a2) = @_;
	my @arr1 = @{$a1};
	my @arr2 = @{$a2};
	my %arr1 = map{$_=>1} @arr1;
	my %arr2 = map{$_=>1} @arr2;
	return grep( $arr1{$_}, @arr2 );
}

# for two arrays A and B, returns the array elements that only exist in array A and not in array B
sub _array_minus {
	my ($a, $b) = @_;
	my @a = @{$a};
	my @b = @{$b};
	my %a = map{$_=>1} @a;
	my %b = map{$_=>1} @b;
	return grep ( ! defined $b{$_}, @a );
}


=item smooth_basal_split

After rerooting, one side of the root may have received all of the length of
the branch on which the root was placed. This method smooths that out as much
as possible.

=cut

sub smooth_basal_split {
	my ( $self, $tree ) = @_;
	my $root = $tree->get_root;

	# compute paths left and right of the root
	my ( $left, $right ) = @{ $root->get_children };
	my ( @hl, @hr, %h );
	for my $tip ( @{ $left->get_terminals } ) {
		push @hl, $tip->calc_path_to_root;
	}
	for my $tip ( @{ $right->get_terminals } ) {
		push @hr, $tip->calc_path_to_root;
	}

	# compute averages
	my $lm = sum(@hl)/scalar(@hl);
	$log->debug("mean height left: $lm");
	my $rm = sum(@hr)/scalar(@hr);
	$log->debug("mean height right: $rm");
	my $diff = abs($rm-$lm) / 2;
	my $r_length = $right->get_branch_length;
	my $l_length = $left->get_branch_length;

	if ( $lm < $rm ) {

		# don't want negative branches
		if ( $r_length < $diff ) {
			$diff = $r_length;
		}

		# adjust branch lengths
		$left->set_branch_length( $l_length + $diff );
		$right->set_branch_length($r_length - $diff );
		$log->info("Stretched left, shrunk right, by $diff");
	}
	else {

		# don't want negative branches
		if ( $l_length < $diff ) {
			$diff = $l_length;
		}

		# adjust branch lengths
		$left->set_branch_length( $l_length - $diff );
		$right->set_branch_length($r_length + $diff );
		$log->info("Stretched right, shrunk left, by $diff");
	}
	$root->set_branch_length(undef);
}

=item outgroup_root

Roots a tree on a list of outgroup taxon names. Arguments:

 -tree     => input tree object
 -outgroup => [ list of names ]
 ( optional: -ids => [ list of outgroup IDs ] )
 -ranks    => [ list of taxonomic ranks of interest ]
 -records  => [ list of taxa table records ]

=cut

sub outgroup_root {
	my ( $self, %args ) = @_;
	my $tree    = $args{'-tree'};
	my @records = @{ $args{'-records'} };
	my @ranks   = $args{'-ranks'} ? @{ $args{'-ranks'} } : qw[forma varietas subspecies species];
	my @ids     = $args{'-ids'} ? @{ $args{'-ids'} } : ();
	my @names   = $args{'-outgroup'} ? @{ $args{'-outgroup'} } : ();
	my $mt = Bio::SUPERSMART::Domain::MarkersAndTaxa->new;
	my $mts = Bio::SUPERSMART::Service::MarkersAndTaxaSelector->new;

	# remap outgroup species names to taxon identifiers
	if ( not @ids ) {
		$log->info("Resolving outgroup name(s): @names");
		@ids = map { $_->ti } $mts->get_nodes_for_names(@names);
	}
	my @all_ids = $mt->query_taxa_table(\@ids, \@ranks, @records);

	# fetch outgroup nodes and mrca
	my %og = map{ $_=>1 } @all_ids;
	my @ognodes = grep { exists($og{$_->get_name}) } @{$tree->get_terminals};

	if ( ! scalar @ognodes ) {
		my @specnames = map{$self->find_node($_)->taxon_name} @all_ids;
		$log->warn("Cannot reroot at outgroup! None of the following names are in the tree: " . join(', ', @specnames));
		return;
	}

	my $mrca = $tree->get_mrca(\@ognodes);
	if ($mrca->is_root) {
		my @ignodes = grep { ! $og{$_->get_name} } @{$tree->get_terminals};
		$mrca = $tree->get_mrca(\@ignodes);
		$log->info("Rooting below MRCA of ingroup");
	}
	else {
		$log->info("Rooting below MRCA of outgroup");
	}

	# set previous root edge to zero
	$tree->get_root->set_branch_length(0.00);

	# reroot under mrca of outgroup
	$mrca->set_root_below;
}

=item remove_internal_names

Removes non-numeric labels of internal nodes in a tree, e.g. "root" or
"e1", which may have been introduced by rerooting or resolving the tree

=cut

sub remove_internal_names {
	my ($self, $tree) = @_;
	$tree->visit(sub{
		my $n = shift;
		$n->set_name('') if $n->is_internal and $n->get_name =~ /[a-zA-Z]/;
	});
	return $tree;
}

=item remap_to_name

Given an object of class L<Bio::Phylo::Forest::Tree>,
changes the names of the terminal nodes from NCBI taxonomy identifiers to their respective taxon names.

=cut


sub remap_to_name {
       	my ($self, $tree) = @_;

	my $mts = Bio::SUPERSMART::Service::MarkersAndTaxaSelector->new;

        $tree->visit(sub{
                my $n = shift;
                if ( $n->is_terminal ) {
                        my $id = $n->get_name;
                        my $dbnode = $self->find_node($id);
                        $log->fatal("Could not find name for taxon id $id in database!") if not $dbnode;
                        my $name = $dbnode->taxon_name;
                        $name = $mts->encode_taxon_name($name);
                        $n->set_name( $name );
                }
                     });
        return $tree;
}

=item remap_newick

Given a newick input tree file, an output file name and a mapping hash, maps all tip names
to hash values. This method is intended for larger sets of trees that should not be read
in memory. Each tree should be on a single line.

=cut

sub remap_newick {
	my ( $self, $infile, $outfile, %map ) = @_;

	open my $outfh, '>', $outfile or die $!;
	open my $infh, '<', $infile or die $!;
	while(<$infh>) {
		chomp;
		if ( /(\(.+;)/ ) {
			my $newick = $1;
			my %pos;
			while ( $newick =~ /[\(,]([^\(,]+?):/g ) {
				my $key = $1;
				$pos{$key} = pos($newick) - ( 1 + length($key) );
				$log->warn("Unknown label '$key' at position ".$pos{$key}) if !$map{$key};
			}
			$log->error("Didn't parse any identifiers from $infile") if not scalar keys %pos;
			for my $key ( sort { $pos{$b} <=> $pos{$a} } keys %pos ) {
				substr $newick, $pos{$key}, length($key), $map{$key};
			}
			print $outfh $newick, "\n";
		}
	}
}

=item remap_to_ti

Given an object of class L<Bio::Phylo::Forest::Tree>,
changes the names of all terminal nodes from taxon names to their respective identifiers
as given in the NCBI taxonomy database. The records argument is optional.

=cut

sub remap_to_ti {
    my ($self, $tree, @records) = @_;

	my $mts = Bio::SUPERSMART::Service::MarkersAndTaxaSelector->new;

	# no taxa table given, we will query the database
	$tree->visit(sub{
		my $n = shift;
		if ( $n->is_terminal and my $name = $n->get_name ) {
			$self->logger->debug("remapping taxon name $name");
			my $ti;
			$name = $mts->decode_taxon_name($name);
			$self->logger->debug("decoded name: $name");
			# if taxa table is given, get ids from there
			if (@records) {
				# valid ranks for tip labels
				my @ranks =  ('forma', 'varietas', 'subspecies', 'species');
				RECORD: for my $rec (@records) {
					my $curr_name = $rec->{"name"};
					if ( $curr_name and $curr_name eq $name ) {
						for my $rank (@ranks) {
							if ( $rec->{$rank} and $rec->{$rank} ne 'NA' ) {
								$ti =  $rec->{$rank};
								last RECORD;
							}
						}
					}
				}
			}
			# no taxa table given or taxon not found in table:  search in database
			if ( ! $ti ) {
				my $dbnode = $self->find_node({taxon_name=>$name});
				die "could not find database entry for taxon name $name " if not $dbnode;
				$ti = $dbnode->ti;
			}
			# set taxon name
			$n->set_name( $ti );
			$log->debug("Remapped name $name to $ti ");
		}
	});

    return $tree;
}


=item make_mapping_table

Given a tree with taxon names, creates a hash mapping taxon ids
to taxon names

=cut

sub make_mapping_table {
	my ( $self, $tree ) = @_;

	my $mts = Bio::SUPERSMART::Service::MarkersAndTaxaSelector->new;
	my %mapping;
	for my $t ( @{$tree->get_terminals} ) {
		my $name = $t->get_name;
		$self->logger->debug("remapping taxon name $name");
		$name = $mts->decode_taxon_name($name);
		$self->logger->debug("decoded name: $name");
		my $dbnode = $mts->find_node({taxon_name=>$name});
		die "could not find database entry for taxon name $name " if not $dbnode;
		my $ti = $dbnode->ti;
		$mapping{$ti} = $t->get_name;
	}

	return %mapping;
}


=item remap

Given a hash with an identifier/taxonname mapping,
and a tree, remaps the leave names of the tree

=cut

sub remap {
	my ( $self, $tree, %mapping ) = @_;

	for my $t ( @{$tree->get_terminals} ) {
		my $mapped = $mapping{$t->get_name};
		$t->set_name($mapped);
	}
	return $tree;
}

=item newick2nexus

Given an input newick file and an output nexus file name, writes trees
from input file as a trees block (without translation table) to the output file.
Returns number of trees so written.

=cut

sub newick2nexus {
    my ( $self, $infile, $outfile ) = @_;
    my $date = localtime();
    open my $in, '<', $infile or die $!;
    open my $out, '>', $outfile or die $!;
    print $out <<'HEADER';
#NEXUS
[! written on $date by $self from $infile ]
BEGIN TREES;
HEADER
    my $i = 0;
    while(<$in>) {
        chomp;
        if ( /^\s*(\(.+;)/ ) {
            my $tree = $1;
            print $out "\t", 'tree tree', ++$i, ' = ', $tree, "\n";
        }
    }
    print $out 'END;';
    return $i;
}

=item consense_trees

Given a file with multiple trees (in newick format),
creates a consensus tree using the tool specified by
"TREEANOTATOR_BIN" in the configuration file.
Returns a single newick tree as a string. Arguments:

 -infile (required):  a nexus tree file to be read by treeannotator
 -burnin (optional):  a number indicating the fraction of trees to discard,
                      default is provided by the BURNIN parameter in the
                      configuration file
 -heights (optional): how to process node heights, an option of 'keep',
                      'median', 'mean' or 'ca', default is provided by the
                      NODE_HEIGHTS parameter in the configuration file
 -limit (optional):   minimum support for a node to be retained
 -format (optional):  input file format. Default is nexus; if newick, performs
                      conversion via newick2nexus

=cut

sub consense_trees {
	my ( $self, %args ) = @_;
	my $burnin  = $args{'-burnin'}  || $config->BURNIN;
	my $infile  = $args{'-infile'}  || die "Need -infile argument!";
	my $heights = $args{'-heights'} || $config->NODE_HEIGHTS;
	my $limit   = $args{'-limit'}   || 0.0;
	my $format  = $args{'-format'}  || 'nexus';

	# do conversion if needed
	if ( lc($format) eq 'newick' ) {
		my ( $fh, $filename ) = tempfile();
		close $fh;
		$self->newick2nexus( $infile => $filename );
		$infile = $filename;
	}
	# we first need to count the absolute number of trees in the file
	# and multiply that by burnin to get the absolute number of trees to skip
	my $counter = 0;
	open my $fh, '<', $infile or die $!;
	while(<$fh>) {
		$counter++ if /^\s*tree\s/i; #*
	}
	my $babs = int( $burnin * $counter );

	# create temporary outfile name
	my ( $outfh, $outfile ) = tempfile();
	close $outfh;

	# execute command
	my $tmpl = '%s -burnin %i -heights %s -limit %f %s %s';
	my $command = sprintf $tmpl, $config->TREEANNOTATOR_BIN, $babs, $heights, $limit, $infile, $outfile;
	$log->debug("running command $command");
	system($command) and die "Error building consensus: $?";

	# parse result
	my $tree = parse_tree(
		'-format' => 'figtree',
		'-file'   => $outfile,
	    );
	unlink $outfile;
	return $tree;
}

=item parse_newick_from_nexus

Given the file name of a phylogenetic tree in Nexus format, creates a tree in newick
representation which is returned as a string. If information about the posterior value
of a node is given (as present in the Nexus output of *BEAST), the posterior value is
assigned to the respective inner nodes in the newick tree.

=cut

sub parse_newick_from_nexus {
        my ($self, $nexusfile, %args ) = @_;
        my $newick;
        open my $fh, '<', $nexusfile or die $!;
        my @lines = <$fh>;
        close $fh;
        my @rev = reverse @lines;

        foreach ( @lines ) {
		if ( /^\s*tree/i ) {
			$newick = $_;
			last;
		}
	}
	die("no tree block found in nexus file $nexusfile") if not $newick;
	$newick =~ s/^\s*tree.+\[\&R\]\s+//gi;

        # get ids for taxon labels from 'Translate' section in nexus file
        my @sub;
        foreach (@lines) {
                push(@sub, split) if (/Translate/ .. /\;/);
        }
        shift @sub;
        pop @sub;
        my %id_map = @sub;

        # remove trailing commas from idenifiers, if any
        for ( values %id_map ) { s/,$//g };

        # parse comments [things between square brackets] to
        # set the posterior as node name
        my @comments = ( $newick =~ m/(\[.+?\])/g );

	# iterate through comments and parse out posterior value, if given
        foreach my $c (@comments){
       	my @matches = ( $c =~ m/posterior=([0-9.]+)/g);
  			$log->warn("Found more than one posterior in comment tag") if scalar(@matches) > 1;
  			my $posterior = $matches[0] || "";
  			my $quoted = quotemeta $c;

	                # round posterior value to two digits after the comma and append maximum value '/1'
	                $posterior = sprintf("%.2f", $posterior) . '/1' if $posterior;
  			# substitute in newick tree string. keep comments so that downstream
  			# parsers may be able to do something with them
  			$newick =~ s/$quoted/$posterior$quoted/g;
        }

        # substitute nexus taxon ids with real taxon labels
        # note that there is possible trouble if node names for posteriors (e.g. 1)
        # overlap with nexus identifier. However, BEAST seems to write all posteriors
        # as proper decimals
	my $newicktree = parse(
            '-string' => $newick,
            '-format' => 'newick',
            %args
        )->first;

  	    $newicktree->visit( sub{
                my $n = shift;
                $n->set_name( $id_map{$n->get_name} ) if exists $id_map{$n->get_name};
	});

        return $args{'-id_map'} ? ( $newicktree, %id_map ) : $newicktree;

}

=item process_commontree

Reconcile the common tree with the taxa in the supermatrix and prepare it
for usage (remove unbranched internal nodes, randomly resolve polytomies, deroot)

=cut

sub process_commontree {
    my ( $self, $commontree, @tipnames ) = @_;

    # only retain tips in supermatrix

	# filter out tipnames that are not in the tree	
	my %tnames = map { $_->get_name=>1 } @{ $commontree->get_terminals };

	$log->debug("Taxa in current supermatrix:");
	$log->debug(Dumper(\@tipnames));

	$log->debug("Taxa in current commontree:");
	$log->debug(Dumper(\%tnames));
	@tipnames = grep { exists $tnames{$_} } @tipnames;
	
	$log->debug("Taxa to keep in commontree:");
	$log->debug(Dumper(\@tipnames));
	
    $commontree->keep_tips( \@tipnames );
	
    # it can occur that a taxon that has been chosen as an exemplar is
    # not in the classification tree. For example, if a species has one subspecies,
    # but the species and not the subspecies is an exemplar. Then this node
    # is an unbranched internal and not in the classification tree. We therefore
    # add these node(s) to the starting tree
    my @terminals = @{ $commontree->get_terminals };

    if ( @terminals != @tipnames ) {
        $log->warn("Tips mismatch: ".scalar(@tipnames)."!=".scalar(@terminals));

        # insert unseen nodes into starting tree
        my %seen = map { $_->get_name => 1 } @terminals;
        my ($p)  = @{ $commontree->get_internals };
        for my $d ( grep { ! $seen{$_} } @tipnames ) {
            $log->warn("Adding node $d (" . $self->find_node($d)->get_name . ") to starting tree");
            my $node = $fac->create_node( '-name' => $d );
            $node->set_parent($p);
            $commontree->insert($node);
        }
    }

    # finalize the tree
    $commontree->resolve->remove_unbranched_internals->deroot;
    return $commontree;
}

=item make_random_starttree

Given an array reference of tip names, simulates an equiprobable starting tree. Returns
tree object.

=cut

sub make_random_starttree {
    my ( $self, $tipnames) = @_;
    my $taxa = $fac->create_taxa;
    my $tree = $fac->create_tree;
    my $rootnode = $fac->create_node( '-name' => 'root' );
    $tree->insert($rootnode);
    for my $t (@{$tipnames}) {
        my $node = $fac->create_node( '-name' => $t, '-branch_length' => 0 );
        $node->set_parent($rootnode);
        $tree->insert($node);
    }
    $tree->resolve;
    $log->info($tree->to_newick);
    return $tree;
}

=item make_ncbi_tree

Given a taxa table, creates a tree with classifications according
to the NCBI taxonomy database

=cut

sub make_ncbi_tree {
	my ( $self, @taxatable) = @_;

	my $mts = Bio::SUPERSMART::Service::MarkersAndTaxaSelector->new;

    # instantiate nodes from infile
	my @nodes = $mts->get_nodes_for_table( @taxatable );

	# compute classification tree
	my $tree = $mts->get_tree_for_nodes(@nodes);

	# create node labels with taxon names
	$tree->visit(sub{
		my $node = shift;
		my $label = $node->get_name;
		$label =~ s/_/\|/g;
		$node->set_name( $label );
	});
	$tree->remove_unbranched_internals;
	return ($tree);
}

=item make_usertree

Given a supermatrix and (optionally) a classification tree, makes a starting tree either
by simulation or from input tree. Requires name of output file to write.

=cut

sub make_usertree {
    my ( $self, $supermatrix, $commontree, $outfile ) = @_;
    my @tipnames = $self->read_tipnames($supermatrix);

    # this will be the tree object to write to file
    my $tree;

    # process the common tree to the right set of tips, resolve it,
    # remove unbranched internals (due to taxonomy) and remove the root
    if ( $commontree ) {
        $log->info("Going to prepare starttree for usage");
        $tree = $self->process_commontree($commontree,@tipnames);
    }

    # simulate a random tree
    else {
        $log->info("No starttree given as argument. Generating random starting tree.");
        $tree = $self->make_random_starttree(\@tipnames);
    }

    # write to file
    open my $fh, '>', $outfile or die $!;
    print $fh $tree->to_newick();
    return $outfile;
}

=item graft_tree

Grafts a clade tree into a backbone tree, returns the altered backbone.

=cut

sub graft_tree {
	my ( $self, $backbone, $clade, $squish ) = @_;

	# retrieve the tips from the clade tree that also occur in the backbone, i.e. the
	# exemplars, and locate their MRCA on the backbone
	my @ids = map { $_->get_name } @{ $clade->get_terminals };
	my ( @exemplars, @clade_exemplars );
	for my $clade_tip ( @{ $clade->get_terminals } ) {
		if ( my $bb_tip = $backbone->get_by_name($clade_tip->get_name) ) {
			push @exemplars, $bb_tip;
			push @clade_exemplars, $clade_tip;
			$log->info("Found shared taxon: ".$clade_tip->get_name);
		}
		else {
			$log->debug("Only in clade: ".$clade_tip->get_name);
		}
	}
	
	# it is possible that none or only one exemplars are found in the backbone tree,
	# e.g. when there is no marker overlap of backbone exemplars with other clade taxa
	# or if an exemplar does have less markers than CLADE_TAXON_MON_MARKERS;
	# in this case, do not try to graft because mrca depth will be zero   
	if ( scalar (@clade_exemplars) < 2 ) {
		$log->error("Less than 2 exemplars found in clade tree, probably too few markers for exemplars available. Skipping clade.");
		return $backbone;
	}
	else {
		$log->info("Found ".scalar(@exemplars). " backbone exemplars");
		$log->info("Found ".scalar(@clade_exemplars). " clade exemplars")
	}

	# find MRCA of exemplars on backbone and clade
	my $bmrca = $backbone->get_mrca(\@exemplars);
	my $cmrca = $clade->get_mrca(\@clade_exemplars);
	if ( $bmrca->is_root ){
		$log->fatal("MRCA of exemplar species " . join(',', map{$_->id} @exemplars) . " in backbone is the backbone root!");
		return $backbone;
	}
	if ( my $min = $bmrca->get_meta_object('fig:fossil_age_min') and my $max = $bmrca->get_meta_object('fig:fossil_age_max') ) {
		$log->info("Backbone MRCA is calibration point ($min..$max)");
		$cmrca = $clade->get_root;
	}

	# calculate the respective depths, scale the clade by the ratios of depths
	my $cmrca_depth = $cmrca->calc_max_path_to_tips;
	my $bmrca_depth = $bmrca->calc_max_path_to_tips;
	$log->info("Clade depth: $cmrca_depth " . $cmrca->to_newick);
	$log->info("Backbone MRCA depth: $bmrca_depth " . $bmrca->to_newick);
	
	if ( $cmrca_depth > 0 ) {
		$self->_rescale( $clade, $bmrca_depth / $cmrca_depth );
		$log->debug("multiplied clade tree by ". ($bmrca_depth/$cmrca_depth));
	}

	# calculate the depth of the root in the clade, adjust backbone depth accordingly
	my $crd  = $clade->get_root->calc_max_path_to_tips;
	my $diff = $bmrca_depth - $crd;
	my $branch_length = $bmrca->get_branch_length || 0;
 	$bmrca->set_branch_length( $branch_length + $diff ); # XXX still gives negative branches
	$log->debug("adjusted backbone MRCA depth by $bmrca_depth - $crd");

	# copy the clade annotation, if any
	if ( my $value = $clade->get_root->get_meta_object('fig:clade') ) {
		$bmrca->set_namespaces( 'fig' => _NS_FIGTREE_ );
		$bmrca->set_meta_object( 'fig:clade' => $value );
	}

	# now graft!
 	$bmrca->clear();
 	$clade->visit(sub {
		my $node = shift;
		my $name = $node->get_internal_name;
		if ( my $p = $node->get_parent ) {
			if ( $p->is_root ) {
				$log->info("Grafted node $node onto backbone MRCA");
				$node->set_parent($bmrca);
			}
		}
		if ( $node->is_root ) {
			$log->info("Replacing root $cmrca with backbone MRCA");
		}
		else {
			$backbone->insert($node);
		}
	});
	return $backbone;
}

sub _rescale {
	my ( $self, $tree, $ratio ) = @_;
	$tree->visit(sub{
		my $node  = shift;

		# select figtree annotations
		my $re = qr/^fig:.*(?:height|length).*/;
		my %annos = map { $_->get_predicate => $_ }
		           grep { $_->get_predicate =~ $re }
		               @{ $node->get_meta };

		# iterate over predicates
		for my $predicate ( keys %annos ) {
			my $val = $annos{$predicate}->get_object;
			my $newval;
			if ( $predicate =~ /(.+?)_(min|max)$/ ) {
				my ( $stem, $range ) = ( $1, $2 );
				my $other = $stem . ( $range eq 'min' ? '_max' : '_min' );
				my $otherval = $annos{$other}->get_object;
				my $mid = ( $val + $otherval ) / 2;
				my $deviation = $mid > 0 ? $val / $mid : 0;
				$newval = ( $mid * $ratio ) * $deviation;
			}
			else {
				$newval = $val * $ratio;
			}
			$annos{$predicate}->set_triple( $predicate => $newval );
			$self->logger->debug("$predicate: $val => $newval");
		}

		# adjust branch length
		my $length = $node->get_branch_length || 0;
		$length *= $ratio;
		$node->set_branch_length($length);
	});
}

=item heights_to_lengths

Given node heights as defined by fig:height annotations, computes the branch lengths

=cut

sub heights_to_lengths {
	my ( $self, $tree ) = @_;
	$tree->visit(sub{
		my $node = shift;
		my $height = $node->get_meta_object('fig:height');
		if ( not defined $height ) {
			$log->warn("No fig:height annotation on $node");
		}
		else {
			$node->set_generic( 'age' => $height );
		}
	});
	$tree->agetobl;
	return $tree;
}

=item make_phylip_binary

Given the location of a file in phylip format, creates a binary representation of
the phylip file and writes it to the specified file.

=cut

sub make_phylip_binary {
	my ( $self, $phylip, $binfilename, $parser, $work_dir) = @_;
	$log->info("Going to make binary representation of $phylip => $binfilename");
	$log->debug("using parser $parser");
	my ( $phylipvolume, $phylipdirectories, $phylipbase ) = File::Spec->splitpath( $phylip );
	my $curdir = getcwd;
	chdir $work_dir;

        # check if filename exists in working directory, if not take the full path
        my $phylipfile = -e $phylipbase ? $phylipbase : $phylip;

        my @command = ( $parser,
		'-m' => 'DNA',
		'-s' => $phylipfile,
		'-n' => $binfilename,
		'>'  => File::Spec->devnull,
		'2>' => File::Spec->devnull,
	);
	my $string = join ' ', @command;
	$log->debug("going to run '$string' inside " . $work_dir );
	system($string) and $log->warn("Couldn't execute command '$string': $! (errno: $?)");
	chdir $curdir;
	return "${binfilename}.binary";
}

=item make_phylip_from_matrix

Given an object of class L<Bio::Phylo::Matrices::Datum>, writes the data into a phylip
file to the specified file name.

=cut

sub make_phylip_from_matrix {
	my ( $self, $taxa, $phylipfile, @matrix ) = @_;

	# create phylip file for parser
	open my $phylipfh, '>', $phylipfile or die $!;
	my %nchar_for_matrix;
	my $ntax  = $taxa->get_ntax;
	my $nchar = 0;
	for my $m ( @matrix ) {
		my $mid = $m->get_id;
		$nchar_for_matrix{$mid} = $m->get_nchar;
		$nchar += $nchar_for_matrix{$mid};
	}
	print $phylipfh $ntax, ' ', $nchar, "\n";
	$taxa->visit(sub{
		my $t = shift;
		my @d = @{ $t->get_data };
		print $phylipfh $t->get_name, ' ';
		for my $m ( @matrix ) {
			my $mid = $m->get_id;
			my ($row) = grep { $_->get_matrix->get_id == $mid } @d;
			if ( $row ) {
				print $phylipfh $row->get_char;
			}
			else {
				print $phylipfh '?' x $nchar_for_matrix{$mid};
			}
		}
		print $phylipfh "\n";
	});
	return $phylipfile;
}


=item read_tipnames

Reads the supermatrix (phylip format) and returns the tip names from it

=cut

sub read_tipnames {
	my ($self, $supermatrix) = @_;
	my $ntax = 0;
	my $line = 0;
	my @result;
        open my $fh, '<', $supermatrix or die $!;
	LINE: while(<$fh>) {
		chomp;
		my $word;
		if ( /^\s?(\S+)/ ) {
			$word = $1;
        }
		if ( not $ntax ) {
			$ntax = $word;
			$log->debug("$supermatrix has $ntax taxa");
                        next LINE;
                }
		if ( $word ) {
                        push @result, $word;
                };
		$log->debug("adding taxon $word");
		last LINE if ++$line == $ntax;
	}
	return @result;
}

=item extract_clades

Given a backbone tree and a list with the taxa classifications,
decomposes the backbone tree into single clades, represented by lists of
taxon ids. The decomposition identifies monophyletic genera which will form individual clades.
If genera are paraphyletic in the tree (i.e. the the mrca of two exemplar species in one genus
is not the direct parent of the two exemplars), we traverse up the tree and include all species
from sister clades until the paraphyly is resolved.

=cut

sub extract_clades {
	my ($self, $tree, @records) = @_;
	my %genera;
    my $mt = 'Bio::SUPERSMART::Domain::MarkersAndTaxa';

    # first get all taxa that will be included in the clades
    my @valid_ranks = ("species", "subspecies", "varietas", "forma");
    my $level = $mt->get_highest_informative_level(@records);
    my @highest_taxa = $mt->get_distinct_taxa($level, @records);
    my %all_taxa = map { $_ => 1 } $mt->query_taxa_table(\@highest_taxa, \@valid_ranks, @records);

	# get exemplars for genera
	my %terminal_ids = map { $_ => 1 } map{$_->get_name} @{ $tree->get_terminals };

	foreach my $id ( keys(%terminal_ids) ) {
		my ($genus) = $mt->query_taxa_table($id, "genus", @records);
		if ( ! $genus ) {
			$log->warn("No genus found for taxon id  $id");
			next;
		}
		push @{$genera{$genus}->{'exemplars'}}, $id;
	}

	# get mrca and distance from root for all exemplars in each genus
	foreach my $genus ( keys %genera ) {
		my %exemplar_ids = map {$_=>1} @{$genera{$genus}->{'exemplars'}};
		my @nodes = grep { exists( $exemplar_ids{$_->get_name}) }  @{$tree->get_terminals};
		my $mrca = $tree->get_mrca( \@nodes );

		$genera{$genus}->{'mrca'} = $mrca;
		my $dist = $mrca->calc_nodes_to_root;
		$genera{$genus}->{'dist_from_root'} = $dist;
	}
	# sort genera by depth of exemplar mrca
	my @sorted_genera = sort { $genera{$a}->{'dist_from_root'} <=> $genera{$b}->{'dist_from_root'} } keys (%genera);

	# make sets of leaves that are in the subtree spanned from the respective mrca
	my @sets;
	for my $genus ( @sorted_genera ) {
		my $mrca = $genera{$genus}->{'mrca'};
		my $dist = $genera{$genus}->{'dist_from_root'};

		# replace mrca with parent node to include monotypic genera in sister clade
		#if ($mrca->is_terminal) {
		#	$mrca = @{$mrca->get_ancestors}[0];
		#}

		# get genera within the subtree:
		my @mrca_terminal_ids = map{$_->get_name} @{$mrca->get_terminals};
		my @mrca_genera =  $mt->query_taxa_table(\@mrca_terminal_ids, 'genus', @records);


		# get all species that are contained in the genera of the mrca
		my @s = $mt->query_taxa_table(\@mrca_genera, \@valid_ranks, @records);
		# only include taxa in a set if they are not already part of a clade spanned higher in the hierarchy
		my @remaining_terminals = keys(%all_taxa);

		my @clade = _intersect(\@remaining_terminals, \@s);
		delete @all_taxa{@clade};

		$log->info('Possible clade containing ' . scalar(@clade) . ' species');

		if ( my $support = $mrca->get_meta_object('fig:posterior') || $mrca->get_meta_object('fig:bootstrap') ) {
			$log->info("Support value of clade mrca: $support");
		}
		push @sets, \@clade;
	}

	# All sets with only one taxon represent monotypic genera. Since we have resolved
	#  paraphyly at this point, all monotypic genera nested in other clades are accounted for;
	#  therefore it is safe to remove all non-nested monotypic genera.

	# There is one special case, however: All genera are monotypic, so all species will be
	#  put together into one clade!
	my @num_elems = uniq map{scalar(@$_)} @sets;
	if ( scalar(@num_elems) == 1 and $num_elems[0] == 1 ) {

		# XXX assuming that all species do need to go into one clade, this presumably
		# means that the lists need to be flattened, i.e. inside the map {} block
		# the arrays need to be dereferenced. That wasn't the case, but now it is.
		@sets = [ map{ @$_ } @sets ];
	}
	else {
		# remove monotypic genera from clades and sets that have <3 species and thus cannot be resolved
		@sets = grep { scalar(@$_) > 2 } @sets;
	}
	$log->info("Extracted " . scalar(@sets) . " clades");
	return @sets;
}

=back

=cut

1;
