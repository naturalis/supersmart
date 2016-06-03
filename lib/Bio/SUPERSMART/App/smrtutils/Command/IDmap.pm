package Bio::SUPERSMART::App::smrtutils::Command::IDmap;

use strict;
use warnings;

use Bio::Phylo::IO qw(parse parse_tree unparse);

use Bio::SUPERSMART::Service::TreeService;
use Bio::Phylo::Factory;
use Bio::Phylo::Util::CONSTANT ':objecttypes';

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
		["outfile|o=s", "name of the output tree file (newick format) defaults to $outfile_default", { default=> $outfile_default, arg => "file"}],    	    
		['outformat|f=s', "file format of output tree, defaults to $outformat_default. Supported formats: newick, nexus, figtree (nexus)", { default => $outformat_default, arg => "format" }],
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
	my @trees = $ts->read_tree( '-file'=>$opt->treefile );

	# remap	
	$self->_remap($_) for @trees;
	
	# write to file
	$ts->to_file( 
		'-file'   => $opt->outfile, 
		'-tree'   => \@trees, 
		'-format' => $outformat
		);

	$logger->info("DONE, tree written to " . $opt->outfile ); 
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
