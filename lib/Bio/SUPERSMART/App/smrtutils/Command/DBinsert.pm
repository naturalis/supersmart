package Bio::SUPERSMART::App::smrtutils::Command::DBinsert;

use strict;
use warnings;

use File::Spec;

use Bio::Phylo::IO qw(parse unparse);
use Bio::Phylo::Util::CONSTANT ':objecttypes';

use Bio::SUPERSMART::Service::MarkersAndTaxaSelector;
use Bio::SUPERSMART::Domain::MarkersAndTaxa;

use base 'Bio::SUPERSMART::App::SubCommand';
use Bio::SUPERSMART::App::smrtutils qw(-command);

# ABSTRACT: inserts custom sequences and taxa into the database

=head1 NAME

DBinsert - insert user specified sequences and new taxa into database

=head1 SYNOPSYS

smrt-utils dbinsert [-h ] [-v ] [-w <dir>] [-l <file>] [-y ] [-a <file>] [-l <file>] [-t ] [-p ] [-d ] [-f ] [-g ] 

=head1 DESCRIPTION

This subcommand lets the user insert custom taxa and sequences into the local SUPERSMART database.
Sequences or alignments can be provided in any sequence format that is supported by Bio::Phylo, default 
is FASTA. IMPORTANT: Note that the inserted sequences are not added to any Phylota cluster, they are only added to
the 'seqs' table. 'smrt align' will therefore not find these sequences, if you want to include them into the analysis,
you will have to manually add the files to the 'smrt align' output file listing the alignments. A prefix for generated
accessions can be given, the default is 'SMRT'. 

If the 'generate_fasta' flag is set, alignment files compatible
with the 'smrt' commands are created. This includes the smrt-style fasta definition line 
(e.g. >gi|449785701|seed_gi|449785701|taxon|292712|mrca|292712/1-1839) and the filename beginning with the 
seed gi. Also note that the seed gi is artificial since no clustering as in phylota is performed.

If taxa present in an alignment to add are not present in the database, smrt-utils dbinsert attempts to add them
to the 'nodes' table in the database.

Also, custom taxa can be added if they are provided in a taxa file as produced by 'smrt taxize' or by 'smrt-utils replicate'.
Note that the taxon ids of the added taxa and its hierarchy must be given in the input taxa file, ids are not 
automatically generated. 

At the moment, this subcommand does not provide functionality to remove entries from the database.
Entries can be removed by hand by searching for instance for the prefixes and suffixes generated while inserting:

sqlite> delete from nodes_194 where common_name like "%SUPERSMART%";
sqlite> delete from seqs where acc like "SMRT%";

=cut

sub options {    
	my ($self, $opt, $args) = @_;
	my $format_default = 'fasta';
	my $prefix_default = 'SMRT';
	return (
		['alignment|a=s', "alignment file(s) to insert into database, multiple files should be separatet by commata", { arg => 'file' }],		
		['list|s=s', "list of alignment files to insert into database", { arg => 'file' } ],		
		['taxafile|t=s', "taxa file as produced by 'smrt-utils replicate'; if given, possible artificial taxa are inserted", {}],
		['prefix|p=s', "prefix for generated sequence accessions, defaults to $prefix_default", {default => $prefix_default}],
		['desc|d=s', "description for sequence(s)", {}],
		['format|f=s', "format of input alignemnt files, default: $format_default", { default => $format_default }],
		['generate_fasta|g', "generate FASTA files compatible with smrt commands. FASTA file names will have prefix <seed-gi> and the suffix '-smrt-inserted'. An alignment list containing the new fasta files is created under the name 'aligned-smrt-inserted.txt'. This option is enabled by default.", {default => 1} ]
	    );	
}

sub validate {
	my ($self, $opt, $args) = @_;			
	$self->usage_error('need either alignment, list or taxafile argument') if not ($opt->alignment or $opt->list or $opt->taxafile);
}

sub run {
	my ($self, $opt, $args) = @_;    
	
	my $logger = $self->logger;      	
		
	# insert potential artificial taxa, if given in taxa file
	$self->_insert_taxa($opt->taxafile, $opt->prefix) if $opt->taxafile;		

	my $aln_str = $opt->alignment;
	my $aln_list = $opt->list;

	# insert alignments, if given in argument
	if ( $aln_str || $aln_list ) {
		my @files;
		
		push @files, split(',', $aln_str) if $aln_str;
		
		if ( $aln_list ) {
			$logger->info("going to read alignment file list $aln_list");
			open my $fh, '<', $aln_list or die $!;
			my @al = <$fh>;
			close $fh;
			chomp @al; 
			push @files, @al;
		}
		$logger->info( 'Got file names of ' . scalar(@files) . ' files to proecess' );
	
		# fresh outfile for inserted alignments
		my $outfile = 'aligned-smrt-inserted.txt';
		unlink $outfile if ( -e $outfile );

		my @inserted;
		# retrieve matrix object(s) from file(s)
		for my $file ( @files ) {
			$logger->info ("Processing alignment file $file");
			# create matrix object from file
			my $project = parse(
				'-format'     => $opt->format,
				'-type'       => 'dna',
				'-file'       => $file,
				'-as_project' => 1,
			    );
			my ($matrix) = @{ $project->get_items(_MATRIX_) };
			
			# iterate over sequences and insert into database
			my @seqs = @{$matrix->get_entities};
			for my $seq ( @seqs ) {
				my $gi = $self->_insert_seq($seq, $matrix, $opt->prefix, $opt->desc );				
				push @inserted, $gi;
			}
			
			# write sequence with new defline to file if generate_fasta flag is given,
			# also write new file name into list of fasta files									
			if ( $opt->generate_fasta ) {
				# prepend artificial 'seed-gi' to filename.
				# This is needed by 'smrt orthologize', otherwise it won't find
				# the sequences! As the seed gi, we take the first gi that we can find
				my $seed_gi = $inserted[0];
				
				# make filename for output file
				my ($vol, $dirs, $newfile) = File::Spec->splitpath($file);
				#remove current file extension
				$newfile =~ s/\.[^\.]+$//;
				$newfile = "${seed_gi}-${newfile}-smrt-inserted.fa";
				
				$logger->info("Writing alignment to $newfile");
				unparse ( -phylo => $matrix, -file => $newfile, -format=>'fasta' );

				# print filename to outfile
				open my $fh, '>>', $outfile or die $!;
				print $fh $newfile . "\n";
				close $fh;
			}	
		}
		$logger->info("Inserted " . scalar(@inserted) . " sequences into database");   	
	}
	$logger->info("DONE");

}

sub _insert_seq {
	my ($self, $seq, $matrix, $prefix, $desc) = @_;
	
 	my $logger = $self->logger;
	my $mts = Bio::SUPERSMART::Service::MarkersAndTaxaSelector->new;
	
	my $id;
	my $name = $seq->get_name;
	if ( $name=~/^[0-9]+?/ ) {
		$id = $name;
	} 
	else {
		$logger->info("Descriptor of sequence $name does not look like taxon ID, trying to map id");
		$name =~ s/_/ /g;
		$logger->info("Trying to find $name");
		my ($node) = $mts->get_nodes_for_names($name);
		if ( ! $node ) {
			$logger->warn("Could not map taxon name $name");

			# If taxon not found, create a new entry  in the database
			my $ti_new = $self->_insert_taxon($name, $prefix);
			$node = $mts->find_node($ti_new);
		}				
		$id = $node->ti;
		$logger->info("Remapped $name to $id");
	}
	# gather fields required for 'seqs' table in phylota
	my $division = 'INV';
	my $acc_vers = 1;
	my $unaligned = $seq->get_unaligned_char;
	my $length = length($unaligned);
	my $gbrel = '000';
	my $def = $desc || "";
	
	$def .= " -- entry generated by SUPERSMART " . $Bio::SUPERSMART::VERSION . " --";
	
	my %seqids = $mts->generate_seqids($prefix);
	
	my $gi = $seqids{'gi'};
	my $acc = $seqids{'acc'};
	
	my $ti = $id;

	# set accession date to today
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
	my $acc_date = join("-", $year + 1900, $mon, $mday);
	
	$mts->insert_seq({gi=>$gi, ti=>$ti, acc=>$acc, acc_vers=>$acc_vers, length=>length($unaligned), division=>$division, acc_date=>$acc_date, gbrel=>$gbrel, def=>$def, seq=>$unaligned});
	
	# write FASTA file with definition line compatible with other smrt commands
	# since the newly inserted seqs have no cluster,  we set the seed gi to the gi of the sequence
	# and the mrca to the taxon id of the sequence. Set mrca to the taxon's ti.
	my $defline = "gi|$gi|seed_gi|$gi|taxon|$ti|mrca|$ti/1-" . length($unaligned);
	$seq->set_generic('fasta_def_line', $defline);		

	$logger->info("Inserted sequence with gi $gi and accession $acc into database");	
	
	# give back the newly created artificial GI
	return($gi);
}

sub _insert_taxon {
	my ($self, $taxon_name, $prefix) = @_;

	my $logger = $self->logger;
	my $mts = Bio::SUPERSMART::Service::MarkersAndTaxaSelector->new;
	
	# get genus ti
	my @names = split(/\s+/, $taxon_name);
	my $genus = $names[0];
	$logger->info("Genus for taxon $taxon_name : $genus");

	my $n = $mts->search_node({ taxon_name=>$genus })->first;
	
	# insert genus into database if not given
	if ( ! $n ) {
		$logger->warn("No genus found for taxon name $taxon_name. Inserting genus $genus into database.");
		my $ti = $mts->max_ti + 1;		
		$mts->insert_node({ti             => $ti, 
						   taxon_name     => $genus,
						   terminal_flag  => 0, 
						   rank_flag      => 1, 
						   rank           => 'genus', 
						   common_name    => "$prefix: $genus - entry generated by SUPERSMART "  . $Bio::SUPERSMART::VERSION . " -", 
						   ti_anc         => 0, 
					  });		
		$n = $mts->search_node({ taxon_name=>$genus })->first;
	}

	my $ti_genus = $n->ti;
	$logger->info("Taxon ID of genus : $ti_genus");

	# get new taxon id for inserted species
	my $ti_species = $mts->max_ti + 1;
	$logger->info("Taxon id for species $taxon_name to insert is $ti_species");
	
	# other fields 
	my $terminal_flag = 1;
	my $rank_flag = 1;
	my $rank = 'species';		
	# we will code in the common name that this entry was generated by smrt
	my $common_name = "$prefix: $taxon_name - entry generated by SUPERSMART "  . $Bio::SUPERSMART::VERSION . " -";
	
	$logger->info("Inserting node for taxon $taxon_name ($ti_species) into database ");
	$mts->insert_node({ti             => $ti_species, 
					   taxon_name     => $taxon_name,
					   terminal_flag  => $terminal_flag, 
					   rank_flag      => $rank_flag, 
					   rank           => $rank, 
					   common_name    => $common_name, 
					   ti_anc         => $ti_genus, 
					  });		
	return $ti_species;	
}

sub _insert_taxa {
	my ($self, $taxafile, $prefix) = @_;
	
	my $mt = Bio::SUPERSMART::Domain::MarkersAndTaxa->new;
	my $mts = Bio::SUPERSMART::Service::MarkersAndTaxaSelector->new;
	my $logger = $self->logger;
	
	my @records = $mt->parse_taxa_file($taxafile);
	my $counter = 0;
	# insert entries that are not yet in database
	for my $rec(@records) {
		# select lowest rank that has taxon ID and genus; note that we do not need to 
		#  look for subspecies etc., since 'smrt replicate' only
		#  produces artificial species!
		my $taxon_name = $rec->{'name'};
		my @ranks = reverse($mts->get_taxonomic_ranks);
		
		my @tis;
		for my $r ( @ranks ) {
			if ( $rec->{$r} =~ /[0-9]+/ ) {
				push @tis, $rec->{$r};			      			
			}
		}
		if ( not scalar(@tis) ) {
			$logger->warn("cannot process taxon $taxon_name, no taxon IDs found");
			next;
		}
		
		# set 
		my $ti = $tis[0];
		my $ti_anc = $tis[1] || 0;
		
		my $node = $mts->find_node($ti);
		if ( $node ) {
			$logger->info("Taxon $taxon_name ($ti) already in database, skipping.");
			next;
		}

		my $ti_genus = $rec->{'genus'};
		
		# gather other fields to enter
		my $terminal_flag = 1;
		my $rank_flag = 1;
		my $rank = 'species';		

		# we will code in the common name that this entry was generated by smrt		
		my $common_name = "$prefix: $taxon_name - entry generated by SUPERSMART "  . $Bio::SUPERSMART::VERSION . " -";	
		$logger->info("Inserting node for taxon $taxon_name ($ti) into database ");
		$mts->insert_node({ti             => $ti, 
				   taxon_name     => $taxon_name,
				   ti_anc         => $ti_anc, 
				   terminal_flag  => $terminal_flag, 
				   rank_flag      => $rank_flag, 
				   rank           => $rank, 
				   common_name    => $common_name, 
				   #ti_genus       => $ti_genus,
				  });		
		$counter++;
	}	
	$logger->info("Inserted $counter nodes into database");	
}

1;
