package Bio::SUPERSMART::App::smrt::Command::Taxize;

use strict;
use warnings;

use Bio::SUPERSMART::Config;
use Bio::SUPERSMART::Service::MarkersAndTaxaSelector;
use Bio::SUPERSMART::Domain::MarkersAndTaxa;

use Bio::SUPERSMART::App::SubCommand;
use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrt qw(-command);

# ABSTRACT: writes taxa table for given list of taxon names or root taxa

=head1 NAME

Taxize.pm - reconciles taxon names with information from NCBI taxonomy

=head1 SYNOPSYS

$ smrt taxize [-v ] [-w <dir>] -i <file> [-o <file>] [-e <rank>] 

=head1 DESCRIPTION

Performs taxonomic name resolution by mapping a provided list of input names onto the 
NCBI taxonomy. Optionally expands the higher input taxa to the specified lower taxon, 
e.g. to expand a named Order to its constituent species. Produces a table that lists 
all resolved (and expanded) taxa and their higher classification. Optionally, non-binomial 
species names can be filtered out. By default, taxize also filters out taxon names containing
the keywords 'unidentified' and 'environmental_sample'.


=cut

sub options {
	my ( $self, $opt, $args ) = @_;
	my $outfile_default = "species.tsv";
	
	my @ranks = reverse(Bio::SUPERSMART::Service::MarkersAndTaxaSelector->get_taxonomic_ranks);
	return (
		[
			"infile|i=s",
			"file with list of taxon names",
			{ arg => "file", galaxy_in => 1, galaxy_type => "data", galaxy_optional => 1 }
		],
		[
			"outfile|o=s",
			"name of the output file, defaults to '$outfile_default'",
			{ default => $outfile_default, arg => "file", galaxy_out => 1, galaxy_type => "data", galaxy_format => "tabular", galaxy_label => "taxa file"}
		],
		[
			"root_taxa|r=s",
			"one or multiple taxon names (seperated by commata) to be expanded to either 'Species' or to taxonomic rank "
			. "given in parameter 'expand_rank'. Taxon names containing spaces must be enclosed in quotes",
			{ arg => "taxon,taxon,...", galaxy_in => 1, galaxy_type => "text", galaxy_optional => 1 }
		],
		[
			"expand_rank|e=s",
			"rank to which root taxa are expanded. Possible values : " . join(', ', @ranks),
			{ default => 0, arg => "rank", galaxy_in => 1, galaxy_type => "text", galaxy_type => "select", galaxy_value => "species", galaxy_options => \@ranks, }
		],
		[
			"binomials_only|b",
			"filter taxon names that are not binomials",
			{ default => 0, galaxy_in => 1, galaxy_type => "boolean" }
		],
		[
		    "all_ranks|a",
		    "also write higher taxa than 'Species' to taxa file",
		    { default => 0, galaxy_in => 1, galaxy_type => "boolean" }
		]
	);
}

sub validate {
	my ( $self, $opt, $args ) = @_;

	# We only have to check the 'infile' argument.
	#  If the infile is absent or empty, abort
	my $file = $opt->infile;
	my $root = $opt->root_taxa;
	$self->usage_error("need infile and/or root_taxon argument") if not $file and not $root;
	if ( $file ) {
		$self->usage_error("file $file does not exist") unless ( -e $file );
		$self->usage_error("file $file is empty")       unless ( -s $file );
	}
}

sub run {
	my ( $self, $opt, $args ) = @_;

	# collect command-line arguments
	my $expand_rank = $opt->expand_rank;
	my $root_taxa   = $opt->root_taxa;

	# instantiate helper objects
	my $log    = $self->logger;
	my $config = Bio::SUPERSMART::Config->new;
	my $mts    = Bio::SUPERSMART::Service::MarkersAndTaxaSelector->new;
	my $mt     = Bio::SUPERSMART::Domain::MarkersAndTaxa->new;

	my @names;

	if ( my $infile = $opt->infile ) {
		@names = $mt->parse_names_file( $infile );
	}
	if ( $root_taxa ) {
		my @taxa = split(',', $root_taxa);
		push @names, @taxa;
		$log->info( scalar(@taxa) . " root taxa given in argument" );
	}	
	# expand root taxa if argument is provided
	if ( $expand_rank or $root_taxa ) {
		@names = $mts->expand_taxa( \@names, $expand_rank || "species" );
	}
	
	if ( ! scalar (@names)) {
		$log->warn("Could not find any taxa for names list or root taxon in database. No output produced.");
		return 1;
	}
	
	my @taxa_table = $mts->make_taxa_table( '-taxon_names' => \@names, '-binomial' => $opt->binomials_only, 
											'-all_ranks' => $opt->all_ranks, '-expand_rank' => lc($opt->expand_rank) );

	if ( ! scalar (@taxa_table) ) {
		my $msg = 'Taxa table does not contain any rows. ' .
			'Either taxon names could not be resolved or higher level taxa were not expanded. ' . 
			'Use --expand_rank to expand higher-level taxa or --all_ranks to include higher taxa in table.';			
		$log->fatal($msg);
		return 0;
	}

	$mts->write_taxa_file( '-file' => $opt->outfile, '-table' => \@taxa_table ); 


	$log->info("DONE, results written to " . $opt->outfile);
	return 1;
}

1;
