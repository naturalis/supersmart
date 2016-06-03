package Bio::SUPERSMART::App::smrtutils::Command::FormatConv;

use strict;
use warnings;

use Bio::Phylo::IO qw(parse_tree);
use Bio::SUPERSMART::Service::TreeService; 

use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrtutils qw(-command);

# ABSTRACT: Convert tree file to other format

=head1 NAME

FormatConv.pm - Format converter for tree files

=head1 SYNOPSYS

smrt-utils formatconv [-t <file>] [-f <format>] [-o <file>] [-h] [-v] [-w <dir>] [-l <file>] [-y] 

=head1 DESCRIPTION

Takes treefile as argument and converts to user-specified format. Supported are NEXUS/figtree, NEXUS and newick

=cut

sub options {    
	my ($self, $opt, $args) = @_;
	return (
		['tree|t=s', "tree file", { arg => 'file' }],		
		['format|f=s', "file format to convert tree to", { arg => "format" }],
		['outfile|o=s', "filename of output tree", { arg => 'file' }],
		);	
}

sub validate {
	my ($self, $opt, $args) = @_;			
	$self->usage_error('need tree file as argument') if not ( $opt->tree );
	$self->usage_error('need format to convert to') if not ( $opt->format );	
    if ( $opt->format !~ /^(?:newick|nexus|figtree)$/i ) {
        $self->usage_error("only newick and nexus format are supported");
    }
}

sub run {
	my ($self, $opt, $args) = @_;    

	my $logger = $self->logger;      	
	my $ts = Bio::SUPERSMART::Service::TreeService->new; 	

	my @trees = $ts->read_tree( '-file'=>$opt->tree );
	$ts->to_file( '-tree'=>\@trees, '-format'=>$opt->format, '-file'=>$opt->outfile);
	
	$logger->info("DONE, tree written to " . $opt->outfile);
	
}

1;
