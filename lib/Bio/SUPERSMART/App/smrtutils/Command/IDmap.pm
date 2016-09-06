package Bio::SUPERSMART::App::smrtutils::Command::IDmap;

use strict;
use warnings;

use Bio::Phylo::IO qw(parse parse_tree unparse);

use Bio::SUPERSMART::Service::TreeService;
use Bio::SUPERSMART::Service::MarkersAndTaxaSelector;
use Bio::Phylo::Factory;
use Bio::Phylo::Util::CONSTANT ':objecttypes';

use Data::Dumper;

use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrtutils qw(-command);

# ABSTRACT: maps between taxon names and NCBI taxonomy identifiers

=head1 NAME

IDmap.pm - Maps a tree with NCBI taxonomy IDs to taxon names and vice versa

=head1 SYNOPSYS

smrt-utils 

=head1 DESCRIPTION

Maps between taxon names (as in the NCBI taxonomy and phylota) and NCBI taxonomy taxon IDs for a given tree.

=cut

sub options {    
	my ($self, $opt, $args) = @_;
	my $outformat_default = 'newick';
	my $outfile_default = 'tree-remapped.dnd';
	return (
		['treefile|t=s', "tree files", { arg => 'file' }],	
		['informat|i=s', "input format", { arg => 'format' }],			
		["outfile|o=s", "name of the output tree file (newick format) defaults to $outfile_default", { default=> $outfile_default, arg => "file"}],    	    
		['outformat|f=s', "file format of output tree, defaults to $outformat_default. Supported formats: newick, nexus, figtree (nexus)", { default => $outformat_default, arg => "format" }],
		['tnrs|s', "do taxonomic name resolution for taxon names in tree to yield a tree with SUPERSMART compatible names", {}],
	    );	
}

sub validate {
	my ($self, $opt, $args) = @_;			

	$self->usage_error('need tree file as argument') if not ($opt->treefile);
	$self->usage_error('tree file not found or empty') if not (-e $opt->treefile and -s $opt->treefile);
}

sub run {
	my ($self, $opt, $args) = @_;    
	my $logger = $self->logger;      	
	my $ts = Bio::SUPERSMART::Service::TreeService->new;
		
	my $outformat = $opt->outformat;	
	
	# parse tree(s)
	my %args  = $opt->informat ? ( '-format' => $opt->informat ) : ();
	$args{'-file'} = $opt->treefile;
	my @trees = $ts->read_tree( %args );

	if ( $opt->tnrs ) {
		# adjust to names in database
		$self->_resolve_names($_) for @trees;
	}
	else {
		# remap	
		$self->_remap($_) for @trees;
	}

	# write to file
	$ts->to_file( 
		'-file'   => $opt->outfile, 
		'-tree'   => \@trees, 
		'-format' => $outformat
		);

	$logger->info("DONE, tree written to " . $opt->outfile ); 
}

sub _resolve_names {
	my ( $self, $tree ) = @_;
	
	my $mts = Bio::SUPERSMART::Service::MarkersAndTaxaSelector->new;
	my $logger = $self->logger;      	
	
	# Traverse tree and check if names are in database or can be tnrs'ed.
	# Change name accordingly.
	# Keep names that could not be mapped.
	my @unmapped;
	my %all_names = map {$_=>1} map {$mts->decode_taxon_name($_->get_name)} @{$tree->get_terminals};
	$tree->visit(
        sub {
            my $node = shift;
			if ( $node->is_terminal ) {
				my $name =  $mts->decode_taxon_name($node->get_name);
				$logger->info("Trying to resolve name $name");
				my @dbnodes = $mts->get_nodes_for_names( $name );

				if ( ! scalar(@dbnodes) ) {
					$logger->warn("Could not map name $name. Will remove tip from tree.");
					push @unmapped, $name; 
				}
				else {
					# Change name of node
					my $newname = $dbnodes[0]->taxon_name;

					if ( $newname eq $name ) {
						$logger->info("Name $name exists in database.");
					}					
					else {
						if ( ! $all_names{$newname} ) {
							$logger->info("Changing tip name from $name to $newname.");
							$node->set_name($mts->encode_taxon_name($newname));
							delete $all_names{$name};
							$all_names{$newname} = 1;
						}
						else {
							# It can be the case that a name is mapped to some name that is already in the tree
							# In that case, the node must not be inserted.
							$logger->warn("Cannot change tip name $name to $newname because $newname is already in the tree. Removing tip $name.");
							push @unmapped, $name;
						}
					}
					$logger->warn("More than one node for $name found in database.") if scalar(@dbnodes) > 1;					
				}
			}
			$logger->debug(Dumper(\@unmapped));
		}
		);
	
	# Prune unmapped tips, if any
	if ( my $cnt = scalar(@unmapped) ) {
		@unmapped = map{$mts->encode_taxon_name($_)} @unmapped;
		$logger->info("Pruning $cnt unmapped tips from tree.");
		$tree = $tree->prune_tips(\@unmapped);
	}
	
	return $tree;
}

sub _remap {
	my ( $self, $tree ) = @_;
	
	my $ts = Bio::SUPERSMART::Service::TreeService->new;

	# determine if tree has ids or names and then remap	
	my @t = map { $_->get_name } @{ $tree->get_terminals };
	if ( scalar ( grep { /\D/ } @t) ) {				
		$self->logger->info("Remapping from taxon names to taxon identifiers");
		$ts->remap_to_ti( $tree );
	}	
	else {
		$self->logger->info("Remapping from taxon identifiers to taxon names");
		$ts->remap_to_name( $tree );
	}	
	return $tree;
}

1;
