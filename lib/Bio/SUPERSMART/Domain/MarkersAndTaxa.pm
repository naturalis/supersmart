# this is an object oriented perl module

package Bio::SUPERSMART::Domain::MarkersAndTaxa;

use strict;
use warnings;

use Bio::SUPERSMART::Service::MarkersAndTaxaSelector;
use List::Util qw (min max sum);
use List::MoreUtils 'uniq';
use Data::Dumper;
use File::Spec;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::LocatableSeq;
use Bio::Phylo::Matrices::Datum;
use Bio::Phylo::Factory;
use Bio::SUPERSMART::Config;
use Bio::Phylo::Util::CONSTANT ':objecttypes';
use Bio::Phylo::Util::Exceptions 'throw';
use Bio::Phylo::Util::Logger;
use Bio::Phylo::IO 'parse_matrix';
use Moose;

has 'alnfiles'      => ( is => 'ro', 'isa' => 'ArrayRef' );
has 'alignments'    => ( is => 'ro', 'isa' => 'HashRef' );
has 'alns_for_taxa' => ( is => 'ro', 'isa' => 'HashRef' );
has 'taxa_for_alns' => ( is => 'ro', 'isa' => 'HashRef' );
has 'adjacencies'   => ( is => 'ro', 'isa' => 'HashRef' );
has 'candidates'    => ( is => 'ro', 'isa' => 'HashRef' );
has 'species'       => ( is => 'ro', 'isa' => 'ArrayRef' );
has 'min_cover'     => ( is => 'ro', 'isa' => 'Int' );
has 'logger'        => ( is => 'ro', 'isa' => 'Str' );

around BUILDARGS => sub {
	my $orig = shift;
	my $class = shift;

	if ( @_ == 2 && ! ref $_[0] ) {

		# here we prepare constructor args. first we parse the alignment file list.
		my $log  = Bio::Phylo::Util::Logger->new;
		my %args = ( 'logger' => $log );

		my $alndir = shift;
		$args{'alnfiles'} = [ $class->parse_aln_dir($alndir) ];
		$log->debug("read ".scalar(@{$args{'alnfiles'}})." alignments from $alndir");
		$args{'min_cover'} = shift;

		# then we parse the alignments
		$args{'alignments'} = {};
		my @taxa;
		for my $f ( @{ $args{'alnfiles'} } ) {
			my %fasta = $class->parse_fasta_file($f);
			$args{'alignments'}->{$f} = \%fasta;
			push @taxa, $class->get_taxa_from_fasta(%fasta);
			@taxa = uniq @taxa;
		}
		$args{'species'} = \@taxa;

		# then we index the alignments
		$class->_index_alignments(\%args);

		# now call the constructor
		return $class->$orig(%args);
	}
	else {
		return $class->$orig(@_);
	}
};

=head1 NAME

Bio::SUPERSMART::Domain::MarkersAndTaxa - Markers and Taxa

=head1 DESCRIPTION

This package defines a domain-specific object that creates and manages optimal sets of
markers and taxa to build supermatrices. It holds state, which is why it is not a
service class but a domain object. The following read-only properties (i.e. getters) are
populated during object creation provided the constructor is passed an alignment list file
and are therefore available as additional method calls:

=over

=item alnfiles

Returns an array ref of alignment file locations

=item alignments

Returns an array ref of parsed alignments

=item alns_for_taxa

Returns a hash ref whose keys are taxon IDs, values are alignment files for that taxon

=item taxa_for_alns

Returns a hash ref whose keys are alignment files, values are taxon IDs in that alignment

=item adjacencies

Returns a hash ref adjacency matrix

=item candidates

Returns a hash ref of candidate taxa

=item species

Returns an array ref of all potential species to be added to the supermatrix

=back

=head1 METHODS

=over

=item _index_alignments

=cut

sub _index_alignments {
    my ( $class, $args ) = @_;

    # instantiate helper objects
    my $config = Bio::SUPERSMART::Config->new;
    my $mts    = Bio::SUPERSMART::Service::MarkersAndTaxaSelector->new;
    my $log    = Bio::Phylo::Util::Logger->new;

    # make forward and reverse mappings between taxa and alignments
    my %taxa_for_alns;
    my %alns_for_taxa;
    for my $i ( 0 .. $#{ $args->{'alnfiles'} } ) {

        # dereference data
		my %fasta = %{ $args->{'alignments'}->{$args->{'alnfiles'}->[$i]} };

		# grep distinct taxa, store under alignment file name
		my @taxa = uniq $class->get_taxa_from_fasta(%fasta);
		$taxa_for_alns{$args->{'alnfiles'}->[$i]} = \@taxa;

		# store alignment file names for each taxon
		for my $t (@taxa) {
			$alns_for_taxa{$t} = [] if not $alns_for_taxa{$t};
			push @{ $alns_for_taxa{$t} }, $args->{'alnfiles'}->[$i];
		}
    }
    $args->{'taxa_for_alns'} = \%taxa_for_alns;
    $args->{'alns_for_taxa'} = \%alns_for_taxa;
    $log->debug("have ".scalar(keys(%taxa_for_alns))." alignments for ".scalar(keys(%alns_for_taxa))." taxa");

    # get adjacency matrix with taxa connected by markers
    my %adjacency_matrix = $mts->generate_marker_adjacency_matrix(
    	[ values %{ $args->{'alignments'} } ],
    	$args->{'species'},
    );
    $args->{'adjacencies'} = \%adjacency_matrix;

    # prune adjacency matrix: each taxon which has
    #  not sufficient coverage cannot possibly be an exemplar
    my @low_coverage_taxa;
    my $cover = $args->{'min_cover'};

    for my $taxon ( keys %adjacency_matrix ) {

    	# taxon has zero alignments or fewer than coverage
    	# XXX modify this to allow for user taxa
	if ( not $alns_for_taxa{$taxon} or scalar( @{ $alns_for_taxa{$taxon} } ) < $cover ) {
	    push @low_coverage_taxa, $taxon;
	}
    }
    $log->debug("have ".scalar(@low_coverage_taxa)."/".scalar(keys(%adjacency_matrix))." low coverage taxa");
    for my $t ( @low_coverage_taxa ) {

    	# remove forward occurrence in AM
	delete $adjacency_matrix{$t};
        for my $k ( keys %adjacency_matrix ) {

            # remove reverse occurrence
            delete $adjacency_matrix{$k}->{$t};
        }
    }
    $log->debug("adjacency matrix now has ".scalar(keys(%adjacency_matrix))." entries");

    # get all independent subsets of species that are connected by at least
    # one marker and select the largest subset as candidates for exemplars
    my $sets = $mts->get_connected_taxa_subsets( \%adjacency_matrix );
    if ( @$sets ) {
    	my %candidates = map { $_ => 1 } @{ ( sort { scalar(@$b) <=> scalar(@$a) } @$sets )[0] };
    	$args->{'candidates'} = \%candidates;
    }
    else {
    	$log->error("No connected species sets were constructed!");
    	$args->{'candidates'} = {};
    }
}

=item orthologize_cladedir

proforms orthologize on a clade directory with alignments present

=cut

sub orthologize_cladedir {
	my ( $self, %args ) = @_;

	my $dir        = $args{'dir'};
	my $outdir     = $args{'outdir'};
	my $maxdist    = $args{'maxdist'};

	my $alndir = File::Spec->catdir($dir, "alignments");

	my $sg = Bio::SUPERSMART::Service::SequenceGetter->new;

	# collect seed gis from alignment file names
	my @files =  glob("${alndir}/*.fa");
	my @gis =  grep { $_ ne '' } map {$1 if $_=~/\/([0-9]+)-clade[0-9]+\.fa/ } @files;

	$sg->merge_alignments( $maxdist, $dir, $alndir, $outdir, @gis );
}

=item write_clade_matrix

given a clade directory with merged alignments,
enriches with additional haplotypes (if 'enrich' argument given)
and writes a matrix of the merged clade alignments
in the specified format (phylip or nexml).

=cut

sub write_clade_matrix {
	my ($self, %args) = @_;

	my $markersfile = $args{'markersfile'};
	my $outfile = $args{'outfile'};
	my $enrich   = $args{'enrich'};
	my $min_markers = $self->min_cover;
	my $max_markers = $args{'max_markers'};
	my $outformat = $args{'format'};

	# initialize the container objects
	my $log = Bio::Phylo::Util::Logger->new;
	my $ns      = 'http://www.supersmart-project.org/terms#';
    my $factory = Bio::Phylo::Factory->new;
    my $mts     = Bio::SUPERSMART::Service::MarkersAndTaxaSelector->new;
    my $service = Bio::SUPERSMART::Service::TreeService->new;
	my $sg      = Bio::SUPERSMART::Service::SequenceGetter->new;

	my $project = $factory->create_project( '-namespaces' => { 'smrt' => $ns } );
	my $taxa    = $factory->create_taxa;
	$project->insert($taxa);

	my @matrices;

	# process all alignment files

	for my $file ( @{ $self->alnfiles } ) {

		# parse the file, enrich and degap it
		$log->info("Adding alignment $file");
		my $matrix = $self->parse_fasta_as_matrix(
			'-name' => $file,
			'-file' => $file,
			'-taxa' => $taxa,
			);
		$mts->enrich_matrix($matrix) if $enrich;
		$matrix = $mts->degap_matrix($matrix);
		push @matrices, $matrix if $matrix;

	}
	return undef if not @matrices;

	# pick CLADE_MAX_MARKERS biggest alignments
	@matrices = sort { $b->get_ntax <=> $a->get_ntax } @matrices;
	$log->info("Number of matrices : " . scalar(@matrices));
	if ( scalar(@matrices) > $max_markers ) {
		$log->info("Found more alignments for clade than CLADE_MAX_COVERAGE. Using the first $max_markers largest alignments.");
	}

	# keep track of taxon ids, the number of markers for each taxid,
	# because some markers might not be selected and taxa have to be removed
	my %markers_for_taxon;
	for my $i ( 0 .. $max_markers - 1 ) {
			if ( my $mat = $matrices[$i] ) {
				my @ids_for_mat;
				for ( @{ $mat->get_entities } ) {
					my $taxid = $_->get_name;
					$taxid =~ s/_.$//g;
					push @ids_for_mat, $taxid;
				}
				$markers_for_taxon{$_}++ for uniq (@ids_for_mat);
				$project->insert($mat);
			}
	}

	# remove a taxon from matrix if it has none or less markers than given in CLADE_TAXON_MIN_MARKERS
	# also remove all rows in the matrices where this taxon appears
	my ($tax) = @{ $project->get_items(_TAXA_) } ;
	for my $t ( @ {$tax->get_entities} ) {
		my $taxname = $t->get_name;
		my $marker_cnt = $markers_for_taxon{$taxname};
		if ( ! $marker_cnt || $marker_cnt < $min_markers ) {
			$log->info("Removing taxon $taxname, not enough markers" );
			$tax->delete($t);

			# remove rows from matrix containing the taxon
			for my $mat ( @matrices ) {
				for my $row ( @{$mat->get_entities} ) {
					if ( my $taxid = $row->get_name =~ /$taxname/ ) {
						$log->info("Removing row for taxon " . $taxname . " from matrix");
							$mat->delete($row);
					}
				}
			}
		}
	}

	# write table listing all marker accessions for taxa
	my @marker_table = $mts->get_marker_table( @{ $project->get_items(_MATRIX_) } );
	$mts->write_marker_table( $markersfile, \@marker_table );

	# write the merged nexml
	if ( lc $outformat eq 'nexml' ) {
		$log->info("Going to write file $outfile");
		open my $outfh, '>', $outfile or die $!;
		print $outfh $project->to_xml( '-compact' => 1 );
	}

	# write supermatrix phylip
	elsif ( lc $outformat eq 'phylip' ) {
		my @matrices = @{ $project->get_items(_MATRIX_) };
		my ($taxa) = @{ $project->get_items(_TAXA_) };
		$log->info("Going to write file $outfile");
		$service->make_phylip_from_matrix($taxa, $outfile, @matrices);
	}
	return $outfile;
}

=item pick_exemplars

Given a taxa table and a comma-separated list of user taxa and the amount of
exemplar species per genus, returns a list of exemplar species.

=cut

sub pick_exemplars {
	my ( $self, $taxafile, $usertaxa, $cnt_per_genus ) = @_;

	$cnt_per_genus = 9**9**9 if $cnt_per_genus == -1;

	# instantiate helper objects
	my $log = $self->logger;
	my $mts = Bio::SUPERSMART::Service::MarkersAndTaxaSelector->new;

	# read and pre-process taxa
	my @records = $self->parse_taxa_file($taxafile);
	my %user_taxa = $self->parse_user_taxa( $usertaxa, @records );
	my %species_for_genus = $self->get_species_for_genera( @records );

	# iterate over genera
	my %candidates       = %{ $self->candidates };
	my %adjacency_matrix = %{ $self->adjacencies };
	my @alignments       = @{ $self->alnfiles };
	my @exemplars;
	for my $genus ( keys %species_for_genus ) {
		# lookup genus name
		my $gname = '';
		if (my $n = $mts->find_node($genus) ) {
			$gname = $n->taxon_name;
		}
		$log->info("Looking for exemplars in genus $gname ($genus)");

		# select taxa that are in this genus and are in the candidate exemplars
		my %genus_taxa       = map { $_ => 1 } @{ $species_for_genus{$genus} };
		my %genus_candidates = map { $_ => 1 } grep { $candidates{$_} } keys %genus_taxa;
		$log->debug("Initial exemplar candidates : " . Dumper(\%genus_candidates));

		# now only keep the ones that have connection within the genus
		for my $can ( sort keys %genus_candidates ) {
			my %adj = %{ $adjacency_matrix{$can} };
			my @connected_within = grep { $genus_taxa{$_} and $adj{$_} > 0 } keys %adj;
			delete $genus_candidates{$can} unless @connected_within;
		}

		# at this point, we might have lost exemplars in monotypic genera and exemplars
		# that do not share markers with other taxa in the genus; if the latter one
		# is the case, we can only add one single taxon to the exemplars, since two taxa
		# without marker overlap would most likely cause a paraphyly in the backbone
		if ( scalar keys(%genus_candidates) < 2 ) {
			if ( my @valid_taxa = grep { exists $candidates{$_} } keys %genus_taxa ) {
				push @exemplars, $valid_taxa[0];
				$log->info("Added taxon $valid_taxa[0] as (monotypic) exemplar");
			}
			else {
				$log->warn("No exemplars found for genus $gname ($genus) ");
			}
		}

		# if at this point we have two candidates, these are our exemplars
		elsif ( scalar keys(%genus_candidates) <=  $cnt_per_genus ) {
			push @exemplars, keys %genus_candidates;
			$log->info( "Added taxa ".join(',',keys %genus_candidates)." as exemplars" );
		}

		# if we still have more than two candidates, take the one which are furthest
		# apart within this genus with respect to the available sequences
		elsif ( scalar keys(%genus_candidates) > $cnt_per_genus ) {
			$log->info("Found more than $cnt_per_genus candidates, choosing the most distant ones");
			my %distance;
			for my $aln ( sort @alignments ) {
				if ( my $d = $self->calc_aln_distances( $aln, [ sort keys %genus_candidates ] ) ) {
					my %dist = %{$d};

					# pick the most distal pair, weight it by number of pairs minus one
					my ($farthest) = sort { $dist{$b} <=> $dist{$a} } sort( keys %dist );
					$distance{$farthest} += scalar( keys(%dist) ) - 1;
				}
			}
			if ( scalar keys %distance ) {
				my  @specs  = uniq map { split( /\|/, $_ ) } sort { $distance{$b} <=> $distance{$a} } sort keys %distance;
				my $idx = min($cnt_per_genus, scalar(@specs)) - 1;
				my @added = @specs[ 0..$idx];
				push @exemplars, @added;
				$log->info("Added taxa "  . join(',', @added) .  " as exemplars");
			}
			else {

				# if most distant pair cannot be found, add them all
				# (not sure if this case can ever happen...)
				push @exemplars, keys %genus_candidates;
				$log->info( "Could not retrieve distances, added all taxa "
					  . join( ',', keys %genus_candidates )
					  . " as exemplars" );
			}
		}
	}
	return @exemplars;
}

=item select_high_coverage_markers

Choose the n markes (alignemnts) that cover the most species.
If exemplar species are not part of at least BACKBONE_MIN_COVERAGE of
these alignments, they get discarded

=cut

sub select_high_coverage_markers {
    my ( $self, $nmarkers, @exemplars ) = @_;
    # instantiate helper objects
    my $log = $self->logger;
    my $config = Bio::SUPERSMART::Config->new;

    # dereference precomputed data
    my %alns_for_taxa = %{ $self->alns_for_taxa };
    my %taxa_for_alns = %{ $self->taxa_for_alns };

    # reduce aft and tfa hashes to only include exemplars
    my %ex = map { $_ => 1 } @exemplars; # quick lookup
    for my $tax ( keys %alns_for_taxa ) {
        delete $alns_for_taxa{$tax} unless $ex{$tax};
    }
    for my $aln ( keys %taxa_for_alns ) {
        my @taxa = grep { $ex{$_} } @{ $taxa_for_alns{$aln} };
        if ( scalar @taxa ) {
            $taxa_for_alns{$aln} = \@taxa;
        }
        else {
            delete $taxa_for_alns{$aln};
        }
    }

	# one alignment must have at least two exemplar species, 
	# remove alignments from tfa and aft which have less
	for my $aln ( keys %taxa_for_alns ) {
		my @taxa = @{ $taxa_for_alns{$aln} };
		if ( scalar (@taxa) < 2 ) {
			delete $taxa_for_alns{$aln};
			# delete alignments from aft
			for my $tax ( @taxa ) {
				my @alns_red = grep { ! ($_ eq $aln) } @{ $alns_for_taxa{$tax} };
				$alns_for_taxa{$tax} = \@alns_red;								
			}			
		}
	}
	
	# count the number of taxa that are in each alignment, sort the alignments by taxon coverage
	# and select the first $nmarkers alignments
	my %num_taxa_for_alns = map { $_=>scalar(@{$taxa_for_alns{$_}}) } keys(%taxa_for_alns);
	
	my @sorted_alns = sort { $num_taxa_for_alns{$b} <=> $num_taxa_for_alns{$a} } keys %num_taxa_for_alns;
	my $cutoff = min($nmarkers, scalar(@sorted_alns));
	my @high_coverage_alns = @sorted_alns[0..$cutoff-1];
	
	# if an exemplar does not have at least BACKBONE_MIN_COVERAGE alns in 
	# that are part of the high coverage alignments, it is discarded
	my @pruned_exemplars;
	for my $ex ( @exemplars ) {

		my %alns = map { $_=>1 } @{$alns_for_taxa{$ex}};
		my $coverage = scalar( grep { exists($alns{$_}) } @high_coverage_alns );
		if ( $coverage < $self->min_cover ) {
			$log->warn("Discarding exemplar $ex due to low coverage in the selected $nmarkers alignments");
		}
		else {
			push @pruned_exemplars, $ex;
		}		
	}
	
	return \@pruned_exemplars, \@high_coverage_alns;	
}

=item optimize_packing_order

Given an array of exemplars, return an array reference of optimally sorted exemplars
and an array ref of optimally sorted alignments

=cut

sub optimize_packing_order {
    my ( $self, @exemplars ) = @_;

    # instantiate helper objects
    my $log = $self->logger;
    my $config = Bio::SUPERSMART::Config->new;

    # dereference precomputed data
    my %alns_for_taxa = %{ $self->alns_for_taxa };
    my %taxa_for_alns = %{ $self->taxa_for_alns };

    # reduce aft and tfa hashes to only include exemplars
    my %ex = map { $_ => 1 } @exemplars; # quick lookup
    for my $tax ( keys %alns_for_taxa ) {
        delete $alns_for_taxa{$tax} unless $ex{$tax};
    }
    for my $aln ( keys %taxa_for_alns ) {
        my @taxa = grep { $ex{$_} } @{ $taxa_for_alns{$aln} };
        if ( scalar @taxa ) {
            $taxa_for_alns{$aln} = \@taxa;
        }
        else {
            delete $taxa_for_alns{$aln};
        }
    }

	# one alignment must have at least two exemplar species, 
	#  remove alignments from tfa and aft which have less
	for my $aln ( keys %taxa_for_alns ) {
		my @taxa = @{ $taxa_for_alns{$aln} };
		if ( scalar (@taxa) < 2 ) {
			delete $taxa_for_alns{$aln};
			# delete alignments from aft
			for my $tax ( @taxa ) {
				my @alns_red = grep { ! ($_ eq $aln) } @{ $alns_for_taxa{$tax} };
				$alns_for_taxa{$tax} = \@alns_red;								
			}			
		}
	}
	
    # sort alignments for each exemplar by decreasing taxon coverage
    for my $taxon ( @exemplars ) {
        if ( my $alns = $alns_for_taxa{$taxon} ) {
            my @sorted = sort {
                  scalar( @{ $taxa_for_alns{$b} } ) <=>
                  scalar( @{ $taxa_for_alns{$a} } )
            } @{ $alns };
            $alns_for_taxa{$taxon} = \@sorted;
        }
    }

    # sort the exemplar taxa by increasing occurrence in alignments, so rarely sequenced
    # species are treated preferentially by including their most speciose alignments first
    my @sorted_exemplars = sort {
        scalar( @{ $alns_for_taxa{$a} } ) <=>
        scalar( @{ $alns_for_taxa{$b} } )
    } grep { $alns_for_taxa{$_} } @exemplars;

    # now collect the alignments (just as many to give all taxa asufficient coverage!)
    # starting with the least well-represented taxa
    my ( %aln, %seen );
  TAXON: for my $taxon ( @sorted_exemplars ) {
        $log->info("Checking alignment coverage for taxon $taxon");

        # take all its not-yet-seen alignments...
        my @alns = grep { !$aln{$_} } @{ $alns_for_taxa{$taxon} };
        $seen{$taxon} = 0 if not defined $seen{$taxon};

        # add alignments until seen enough. XXX we may have a
        # BACKBONE_MAX_COVERAGE so that all our exemplars have at least
        # a minimum which we, optionally, add to. Alternatively, we might
        # implement some sort of iterative expand() function that computes,
        # given the selected exemplars and alignments, which taxa are most
        # distant from each other in marker sharing, and try to select
        # markers that shorten that distance.
        my $max = $config->BACKBONE_MAX_COVERAGE or $self->min_cover;
      ALN: while ( $seen{$taxon} < $max ) {

            # most speciose alignments first: we sorted aft already
            my $aln = shift @alns;
            if ( not $aln or not -e $aln ) {

                # might happen, since we selected exemplars for which there
                # are sufficient minimum alignments, but not necessarily
                # BACKBONE_MAX_COVERAGE!
                $log->info("No more alignments available for $taxon, now have: ".$seen{$taxon});
                next TAXON;
            }
            $aln{$aln}++;

            # increment coverage count for all taxa in this alignment
            $seen{$_}++ for @{ $taxa_for_alns{$aln} };
			$log->debug("Have now " . $seen{$taxon}  . " alignments for taxon $taxon");
			last ALN if not @alns;
        }
    }
	my @sorted_alignments = sort keys %aln;
	
	# Calculate density of Supersmatrix with respect to marker occurrence
	my $nonzero = sum( values(%seen) );
	my $dens = $nonzero / ( scalar(@sorted_exemplars) * scalar(@sorted_alignments) );
	$log->info("Density of supermatrix : $dens");

    return \@sorted_exemplars, \@sorted_alignments;
}

=item dedup

Removes duplicate sequences (by GI)

=cut

sub dedup {
    my ( $class, %fasta ) = @_;
    my %seen;
    my %result;
    for my $defline ( keys %fasta ) {
        if ( $defline =~ /gi\|(\d+)/ ) {
            my $gi = $1;
            $result{$defline} = $fasta{$defline} unless $seen{$gi}++;
        }
    }
    return %result;
}


=item delete_empty_columns

Removes gap-only columns

=cut

sub delete_empty_columns {
    my ( $self, $fasta ) = @_;
    my $log = $self->logger;
    my %allseqs = %$fasta;

    # Delete columns that only consist of gaps
    my $nchar = length $allseqs{ ( keys(%allseqs) )[0] };
    my $vals  = values %allseqs;
    my %ind;
    for my $v ( values %allseqs ) {
        # check if column only consists of gap characters
        $ind{ pos($v) - 1 }++ while $v =~ /[-Nn\?]/g;
    }
    my @to_remove = sort { $b <=> $a } grep { $ind{$_} == $vals } keys %ind;
    # remove columns at specified indices
    for my $v ( values %allseqs ) {
        substr( $v, $_, 1, "" ) for @to_remove;
	}
    my $removed = $nchar - length $allseqs{ ( keys(%allseqs) )[0] };
    $log->info("Removed $removed gap-only columns from matrix");
	return %allseqs;
}

=item keep_taxa

Retains the records for the provided taxon identifiers that occur in a hash such as is
produced by parse_fasta_string

=cut

sub keep_taxa {
    my ( $class, $taxa, $fasta ) = @_;
    my %result;
    my %taxa = map { $_ => 1 } @{ $taxa };
    for my $defline ( keys %{ $fasta } ) {
        if ( $defline =~ /taxon\|(\d+)/ ) {
            my $taxon = $1;
            $result{$defline} = $fasta->{$defline} if $taxa{$taxon};
        }
    }
    return %result;
}

=item query_taxa_table

given a taxon id or a reference to a list of taxon identifiers, and a single taxonomic
rank or a reference to a list of ranks, returns all entries in the taxon table matching
the rank(s) that are associated with the input taxon ID(s).

=cut

sub query_taxa_table {
    my ($self, $tids, $ranks, @records) = @_;

    # input can either be scalars or array refs
    my %tids = ref($tids) ? map {$_=>1} @${tids} : ($tids=>1);
    my %ranks = ref($ranks) ? map {$_=>1} @${ranks} : ($ranks=>1);

    my @matching = grep { my %h = %{$_}; grep{ exists($tids{$_}) } values%h } @records;

    # XXX why would this be necessary?
    no warnings 'uninitialized';
    my @ids = uniq map { @{$_}{keys(%ranks)} } @matching;
    # remove NA values from result
    my @result = grep { defined $_ and /[0-9]+/ } @ids;
    return @result;
}

=item parse_names_file

Reads a file listing taxon names and returns them in an array

=cut

sub parse_names_file {
    my ( $class, $file ) = @_;
    my $log  = Bio::Phylo::Util::Logger->new;
    $log->info("going to read taxon names from file $file");
	my @names;

	# read names from file or STDIN, clean line breaks
	open my $fh, '<', $file or die $!;
	while(<$fh>) {

		# strip line breaks and leading/trailing whitespace
		chomp;
		s/^\s*//;
		s/\s*$//;
		push @names, $_ if /\S/;
	}
	close $fh;
	$log->info( "Read " . scalar(@names) . " taxon names from $file" );

	return @names;
}

=item parse_aln_file

Reads the flat list of alignments, returns an array of validated file names
(i.e. not empty strings, must be valid locations)

TODO: Deprecated. We are switching to use separate directories instead of listing the
	  alignment fasta file in a separate file.

=cut

sub parse_aln_file {
    my ( $class, $file ) = @_;
    my $log  = Bio::Phylo::Util::Logger->new;
    $log->debug("going to read alignment list $file");
	$log->error("file $file does not exist") if not -e $file;
    my @alignments;
    open my $fh, '<', $file or die $!;
    while (<$fh>) {
        chomp;
        if ( /\S/ ) {
        	if ( -e $_ ) {
        		push @alignments, $_;
        	}
        	else {
        		$log->warn("Entry $_ not found!");
        	}
        }
        else {
        	$log->debug("skipping blank line");
        }
    }
    return @alignments;
}


=item extract_seed_gis

Given an array of file names of the form 'XXXXX-XXXXX-XX-subtree.fa',
extracts the seed_gi from each file, which is by convention the number 
before the first '-'. Returns an array of seed gis

=cut

sub extract_seed_gis { 
	my ( $self, @filenames ) = @_;
	
	my $log    = Bio::Phylo::Util::Logger->new;
	$log->debug("Going to extract seed GIs from filenames");

	my @gis;
	for my $name ( @filenames ) {
		if ( $name =~ /(\d+)-.+\.fa/ ) {
			my $gi = $1;
			push @gis, $gi;
		}
	}
	return @gis;
}

=item parse_aln_dir

Reads a directory of alignments and returns an array with the file locations.
Note that alignments need the ending ".fa" in order to be recognised.
The directory can be a zip file which will then be unzipped.

=cut

sub parse_aln_dir {
	my ( $self, $dirname ) = @_;
	my $log    = Bio::Phylo::Util::Logger->new;
    $log->debug("going to read fasta files in alignment directory '$dirname'" );
	$log->error("directory $dirname does not exist") if not -d $dirname;

	# extract filenames
	my @files = glob( "$dirname/*.fa" );

	return @files;
}

=item parse_user_taxa

Parses the list of user-provided taxa. Arguments: comma-separated list of user taxa,
array of taxon records from species.tsv

=cut

sub parse_user_taxa {
    my ( $self, $include_taxa, @records ) = @_;
    my $log = $self->logger;
    my $mts = Bio::SUPERSMART::Service::MarkersAndTaxaSelector->new;
    my @ranks = qw(species subspecies varietas forma);
    my %user_taxa;
    if ($include_taxa) {
    	$log->info("Going to parse list of user taxa");
        my @taxa = split /,/, $include_taxa;
        my @ids = map {
            my @nodes = $mts->search_node( { taxon_name => $_ } )->all;
            $nodes[0]->ti
        } @taxa;
        my @species = $self->query_taxa_table( \@ids, \@ranks, @records );
        %user_taxa = map { $_ => 1 } @species;
        $log->info("Identified ".scalar(keys(%user_taxa))." species to retain");
    }
    return %user_taxa;
}

=item parse_fasta_file

Reads the provided FASTA file, returns a hash where keys are FASTA definition lines,
values are sequences.

=cut

sub parse_fasta_file {
    my ( $class, $file ) = @_;
    open my $fh, '<', $file or die $!;
    my $string = do { local $/; <$fh> };
    return $class->parse_fasta_string($string);
}

=item parse_fasta_string

Reads the provided FASTA string, returns a hash where keys are FASTA definition lines,
values are sequences.

=cut

sub parse_fasta_string {
    my ( $class, $string ) = @_;
    my @lines = split /\n/, $string;
    my %fasta;
    my $current;
    for my $line ( @lines ) {
        chomp $line;
        if ( $line =~ /^>(.+)/ ) {
            $current = $1;
            if ( exists $fasta{$current} ) {
                $fasta{$current} = '';
            }
        }
        else {
            $fasta{$current} .= $line;
        }
    }
    return %fasta;
}

=item parse_fasta_as_matrix

Reads FASTA file as Bio::Phylo::Matrices::Matrix. Arguments:

 -taxa => taxa block
 -file => input file
 -name => (optional) matrix name

=cut

sub parse_fasta_as_matrix {
    my ( $class, %args ) = @_;
    my $file = $args{'-file'} or throw 'BadArgs' => "Need -file argument";
    my $taxa = $args{'-taxa'} or throw 'BadArgs' => "Need -taxa argument";
    my $factory = Bio::Phylo::Factory->new;

    # read fasta data
    my $matrix = parse_matrix(
        '-type'       => 'dna',
        '-format'     => 'fasta',
        '-file'       => $file,
        '-as_project' => 1,
    );
    $matrix->set_name($args{'-name'}) if $args{'-name'};

    # create taxon lookup
    my %t = map { $_->get_name => $_ } @{ $taxa->get_entities };

    # cleanup rows
    $matrix->visit(sub{
        my $row = shift;
        my $name = $row->get_name;
        my %fields = split /\|/, $name;
        $fields{$_} =~ s/^(\d+).*$/$1/ for keys %fields;
        my $binomial = $fields{'taxon'};

        # create new taxon object if none seen yet
        if ( not $t{$binomial} ) {
            $t{$binomial} = $factory->create_taxon( '-name' => $binomial );
            $t{$binomial}->set_meta_object( 'smrt:tid' => $fields{'taxon'} );
            $taxa->insert($t{$binomial});
        }
        $row->set_taxon( $t{$binomial} );
        $row->set_name( $binomial );
        $row->set_meta_object( 'smrt:gi'      => $fields{'gi'} );
        $row->set_meta_object( 'smrt:mrca'    => $fields{'mrca'} );
        $row->set_meta_object( 'smrt:seed_gi' => $fields{'seed_gi'} );
    });
    $matrix->set_taxa($args{'-taxa'});
    return $matrix;
}

=item parse_taxa_file

Reads the provided taxa file, returns a list of records where each consists of a
hash where keys are the headers in the taxa table.

B<Note>: this method, and C<parse_taxa_string> are completely agnostic about the
column headers so the can equally be used to read other tab-separated files, such
as the backbone markers table, for example.

=cut

sub parse_taxa_file {
    my ( $class, $file ) = @_;
    open my $fh, '<', $file or die $!;
    my $string = do { local $/; <$fh> };
    return $class->parse_taxa_string($string);
}

=item parse_taxa_string

Reads the provided taxa string, returns a list of records where each consists of a
hash where keys are the headers in the taxa table.

=cut

sub parse_taxa_string {
    my ( $class, $string ) = @_;
    my @header;
    my @result;
    LINE: for my $line ( split /\n/, $string ) {
        chomp $line;

        # skip comments and blank lines
        next LINE if $line =~ /^\s*#/;
        next LINE if $line =~ /^\s*$/;

        # store header
        if ( not @header ) {
            @header = split /\t/, $line;
            next LINE;
        }
        my %record = ( 'keys' => [ @header ] );

        # store fields
        my @fields = split /\t/, $line;
        for my $i ( 0 .. $#header ) {
            $record{$header[$i]} = $fields[$i];
        }

        # store record
        push @result, \%record;
    }
    return @result;
}

=item calc_mean_distance

Calculates the average pairwise distance within the alignment.

=cut

sub calc_mean_distance {
    my ( $class, $fastastr ) = @_;
    my $stream = Bio::AlignIO->new(
    	'-string' => $fastastr,
    	'-format' => "fasta",
    );
    my $aln = $stream->next_aln;

	# return 0 if only one sequence provided
	if ( $aln->num_sequences == 1 ) {
		return 0;
	}
	else {
		return ( 1 - $aln->average_percentage_identity / 100 );
	}
}

=item calc_aln_distances

Given an alignment file (fasta) and a reference to an array of taxon ids,
calculates the distance between the taxa with respect to the sequences
given (if taxa share markers). Returned is a hashref with keys being a combination
of taxa, separated by '|', the values being the molecular distance between these taxa,
normalized by sequence length

=cut

sub calc_aln_distances {
    my ( $self, $aln, $tax ) = @_;
    my @taxa  = @$tax;
    my $ts    = Bio::SUPERSMART::Service::TreeService->new;
    my $dat   = 'Bio::Phylo::Matrices::Datum';
    my $log   = $self->logger;
    my %fasta = $self->parse_fasta_file($aln);
    my %dist;

    # aggregate sequence objects by species within this genus
    my %sequences_for_species;
    for my $taxon (@taxa) {
        my $name = $ts->find_node($taxon)->taxon_name;
        if ( my @raw = $self->get_sequences_for_taxon( $taxon, %fasta ) ) {
            $log->debug("$aln contained ".scalar(@raw)." sequences for $taxon ($name)" );
            my @seq = map { $dat->new( '-type' => 'dna', '-char' => $_ ) } @raw;
            $sequences_for_species{$taxon} = \@seq;
        }
    }

    # check if we've seen enough sequenced species
    if ( scalar( keys(%sequences_for_species) ) < 3 ) {
        $log->debug("fewer than 3 species in $aln, no distance calculated");
        return;
    }

    # calculate the distances between give taxa, take the
    # average if species have multiple sequences
    for my $i ( 0 .. ( $#taxa - 1 ) ) {
        for my $j ( ( $i + 1 ) .. $#taxa ) {
            my $sp1 = $taxa[$i];
            my $sp2 = $taxa[$j];
            if ( $sequences_for_species{$sp1} and $sequences_for_species{$sp2} ) {
                $log->debug("going to compute average distance between $sp1 and $sp2");
                my @seqs1 = @{ $sequences_for_species{$sp1} };
                my @seqs2 = @{ $sequences_for_species{$sp2} };
                my $dist;
                for my $seq1 (@seqs1) {
                    for my $seq2 (@seqs2) {
                        $dist += $seq1->calc_distance($seq2);
                    }
                }
                $dist /= ( scalar(@seqs1) * scalar(@seqs2) );
                my $key = join '|', sort { $a <=> $b } ( $sp1, $sp2 );
                $dist{$key} = $dist;
            }
        }
    }
    return \%dist;
}

=item get_alignment_subset

given an alignment as a hash (as for instance produced by parse_fasta_file),
extracts sequences that match a specified keyword in the fasta definition line.
The keyword is specified by the second argument, e.g. {'taxon'=>[10,20,30]}
to extract all entries for taxa 10, 20 and 30, respectively.

=cut

sub get_alignment_subset {
    my ($self, $aln, $arg) = @_;

    # determine for what type (gi, taxon, ..., given in fasta header) to chose
    my ($type) = keys %$arg;
    my @ids = @{$arg->{$type}};

    my %result;
    for my $k ( keys %$aln ) {
        foreach my $id ( @ids ) {
            if ( $k =~ m/$type\|$id/ ) {
                $result{$k} = $aln->{$k};
            }
        }
    }
    return %result;
}

=item get_supermatrix_numtaxa

Given the file location of a supermatrix, returns the number of taxa

=cut

sub get_supermatrix_numtaxa {
    my ($self, $supermatrix) = @_;
    # parse number of sites in alignment from first line in supermatrix
    open my $file, '<', $supermatrix or die $!;
    my $firstline = <$file>;
    $firstline =~ s/^\s+//;
    my $numsites = (split(/\s/, $firstline))[0];
    close $file;
    return $numsites
}

=item get_supermatrix_numsites

Given the file location of a supermatrix, returns the number of sites

=cut

sub get_supermatrix_numsites {
    my ($self, $supermatrix) = @_;
    # parse number of sites in alignment from first line in supermatrix
    open my $file, '<', $supermatrix or die $!;
    my $firstline = <$file>;
    $firstline =~ s/^\s+//;
    my $numsites = (split(/\s/, $firstline))[1];
    close $file;
    return $numsites
}

=item get_supermatrix_taxa

Given the file location of a supermatrix, returns the taxa in the matrix

=cut

sub get_supermatrix_taxa {
    my ( $self, $supermatrix ) = @_;
    my ( $ntax, @taxa );
    open my $fh, '<', $supermatrix or die $!;
    LINE: while(<$fh>) {
        chomp;
        if ( not $ntax and /^\s*(\d+)\s+\d+\s*$/ ) {
            $ntax = $1;
            next LINE;
        }
        if ( @taxa < $ntax and /^\s*(\S+)/ ) {
            my $taxon = $1;
            push @taxa, $taxon;
        }
        last LINE if @taxa == $ntax;
    }
    return @taxa;
}


=item get_distinct_taxa

Returns the distinct IDs for the provided taxonomic level.

=cut

sub get_distinct_taxa {
    my ( $class, $taxon, @records ) = @_;
    my %ids = map { $_ => 1 } map { $_->{$taxon} } @records;
    # filter for e.g. 'NA' values
    my @taxa = grep /[0-9]+/,  keys %ids;
    return @taxa;
}

=item get_species_for_genera

Given a set of taxa table records, returns a hash where each key is a
genus and the values are species in that genus

=cut

sub get_species_for_genera {
    my ( $self, @records ) = @_;
    my $log = $self->logger;
    my %species_for_genus;
    my @ranks = qw(species subspecies varietas forma);
    for my $genus ( $self->get_distinct_taxa( 'genus' => @records ) ) {

        # extract the distinct species (and lower) for the focal genus
        my @species = $self->query_taxa_table( $genus, \@ranks, @records );
        $species_for_genus{$genus} = \@species;
    }
    $log->info( "Read " . scalar(keys(%species_for_genus)) . " genera" );
    return %species_for_genus;
}

=item get_highest_informative_level

Determines, for a given species table, which taxonomic level
(column in the table) is the highest one that can give information about
distinct taxa in this table. This will be the column which does not contain NA values and
which has at least two different taxon identifiers. This level forms a good basis to check
for paraphyly in the branches spawning from there.

=cut

sub get_highest_informative_level{
    my ( $class, @records) = @_;
    my $mts = Bio::SUPERSMART::Service::MarkersAndTaxaSelector->new;
    my $result;

    # get all possible taxonomic ranks, ordered from highest to lowest
    my @all_ranks = $mts->get_taxonomic_ranks;
    my %ranks_in_table = map { $_=>1 } keys %{$records[0]};

    # iterate over all ranks; the first informative taxon level is the level that
    # does not contain NA values and that has distinct taxon IDs
    foreach my $rank (@all_ranks){

        # skip rank if not in our taxa table
        next if not exists $ranks_in_table{$rank};

        # extract column for taxonomic level
        my @column = map {my %h=%$_; my $entry=$h{$rank}; $entry=~s/\s//g; $entry; } @records;

        # omit columns that contain 'NA' values
        next if grep {$_ eq 'NA'} @column;

        if (scalar (uniq @column) > 1){
            $result = $rank;
            last;
        }
    }
    return $result;
}

=item get_root_taxon_level

Returns the level of the 'root' taxon. This is the lowest taxonomic rank which
is the same for all entries in the species table.

=cut

sub get_root_taxon_level {
    my ( $class, @records) = @_;
    my $mts = Bio::SUPERSMART::Service::MarkersAndTaxaSelector->new;
    my $result;

    # get all possible taxonomic ranks, ordered from lowest to highest
    my @all_ranks = reverse $mts->get_taxonomic_ranks;
    my %ranks_in_table = map { $_=>1 } keys %{$records[0]};

    # iterate over all ranks; the root taxon level is the first rank for which all
    #  entries in the taxa table are the same!
    for my $rank ( @all_ranks ){

        # skip rank if not in our taxa table
        next if not exists $ranks_in_table{$rank};

        # extract column for taxonomic level
        my @column = map {my %h=%$_; my $entry=$h{$rank}; $entry=~s/\s//g; $entry; } @records;

        # omit columns that contain 'NA' values
        next if grep {$_ eq 'NA'} @column;
        if (scalar (uniq @column) == 1){
            $result = $rank;
            last;
        }
    }
    return $result;
}

=item get_taxa_from_fasta

Extracts taxon identifiers from a hash such as is produced by parse_fasta_string

=cut

sub get_taxa_from_fasta {
    my ( $class, %fasta ) = @_;
    my @taxa;
    for my $key ( keys %fasta ) {
        if ( $key =~ /taxon\|(\d+)/ ) {
            my $taxon = $1;
            push @taxa, $taxon;
        }
    }
    return @taxa;
}

=item get_gis_from_fasta

Extracts sequence identifiers from a hash such as is produced by parse_fasta_string.
Returns a hash ( taxon => [ gi1, gi2 ] )

=cut

sub get_gis_from_fasta {
    my ( $class, %fasta ) = @_;
    my %gis;
    for my $key ( keys %fasta ) {
        if ( $key =~ /taxon\|(\d+)/ ) {
            my $taxon = $1;
            $gis{$taxon} = [] if not $gis{$taxon};
            if ( $key =~ /(>|\|)gi\|(\d+)/ ) {
                my $gi = $1;
                push @{ $gis{$taxon} }, $gi;
            }
        }
    }
    return %gis;

}

=item get_sequences_for_taxon

Extracts a list of sequence strings for the provided taxon from a hash such as is
produced by parse_fasta_string

=cut

sub get_sequences_for_taxon {
    my ( $class, $taxon, %fasta ) = @_;
    my @sequences;
    for my $key ( sort keys %fasta ) {
        if ( $key =~ /taxon\|$taxon[^\d]/ ) {
            push @sequences, $fasta{$key};
        }
    }
    return @sequences;
}

=item get_nchar

Returns the number of characters in the fasta hash.

=cut

sub get_nchar {
    my ( $class, %fasta ) = @_;
    my ($seq) = values %fasta;
    return length $seq;
}

=item to_fasta_string

Returns the provided fasta hash as a concatenated string.

=cut

sub to_fasta_string {
    my ( $class, %fasta ) = @_;
    return join "\n", map { ">$_\n$fasta{$_}" } keys %fasta;
}

=item write_supermatrix

Given a set of alignment files and exemplar taxa, a file format (default phylip),
and file locations for marker table and supermatrix, writes out a supermatrix to file

=cut

sub write_supermatrix {
    my ( $self, %args ) = @_;

    # instantiate helper objects
    my $logger = $self->logger;
    my $mts    = Bio::SUPERSMART::Service::MarkersAndTaxaSelector->new;
    my $config = Bio::SUPERSMART::Config->new;

    # dereference argument data structures
    my @exemplars  = @{ $args{'exemplars'} };
    my @alignments = @{ $args{'alignments'} };
	# if 'taxon_names' given, print taxon names instead of IDs
	my $taxon_names = $args{'taxon_names'};

    # make hash with concatenated sequences per exemplar
    my %allseqs = map { $_ => "" } @exemplars;

    # iterate over alignments
    my @marker_table;
    for my $aln ( sort { $a cmp $b } @alignments ) {
    	my $fasta = $self->alignments->{$aln};
		my $nchar = length((values %$fasta)[0]);
		my %marker;

		# iterate over taxa
		for my $taxon ( sort { $a <=> $b } @exemplars ) {
			# iterate over sequences
			my @best;
			for my $defline ( keys %$fasta ) {
				my %def = grep { /\S/ } split /\|/, $defline;
				if ( $def{'taxon'} =~ /^$taxon\/?/ ) {
					my $seq = $fasta->{$defline};
					my $missing = ( $seq =~ tr/?/?/ );
					push @best, [ $missing, $def{'gi'}, $seq ];
				}
			}
			# pick best sequence
			if ( @best ) {
				my ( $best ) = sort { $a->[0] <=> $b->[0] } @best;
				$marker{$taxon} = [ $best->[1] ]; # store GI
				$allseqs{$taxon} .= $best->[2]; # grow matrix
			}
			else {
				$allseqs{$taxon} .= '?' x $nchar;
			}
		}
		push @marker_table, \%marker;
    }

    # write table listing all marker accessions for taxa
    $mts->write_marker_table( $args{'markersfile'}, \@marker_table, $args{'exemplars'} );

    # prune gap only columns
	%allseqs = $self->delete_empty_columns(\%allseqs);

    # Write supermatrix to file
    my $aln = Bio::SimpleAlign->new();
	map {
		my $id = $taxon_names ? $mts->find_node($_)->taxon_name : $_;
		$id =~ s/\s/_/g;
        $aln->add_seq(Bio::LocatableSeq->new(
	    '-id'   => $id,
	    'seq'   => $allseqs{$_},
	    'start' => 1
        ));
    } keys %allseqs;
    $aln->sort_alphabetically;

    my $filename = $args{'outfile'};
	open my $fh, '>', $filename or die $!;

	my $format = $args{'format'} || 'phylip';
	$format = 'nexus' if lc $format eq 'mrbayes';
	my $idlength = 10;
	if ( $taxon_names ) {
		$idlength = max ( map { length $_->id } $aln->each_seq );
	}
	my %output_args = (
		-format => $format,
		-fh => $fh,
		-idlength => $idlength,
		-show_symbols => 0,
		-show_endblock => 0,
		);

    my $stream = Bio::AlignIO->new( %output_args );

	# add run block for MrBayes if requested
	my $runblock;
	if ( lc $args{'format'} eq 'mrbayes' ) {
		$aln->missing_char('?');
		my $numruns = $config->EXABAYES_NUMRUNS;
		my $numgens = $config->EXABAYES_NUMGENS;
		my $sfreq = 1000;
		my $file = "mrbayes-out.nex";
		$runblock = "begin mrbayes;\n" .
			"set autoclose=yes nowarn=yes;\n" .
			"lset nst=6 rates=gamma;\n" .
			"mcmc nruns=$numruns ngen=$numgens samplefreq=$sfreq file=$file;\n" .
			"end;\n";
	}
	$stream->write_aln($aln);
	print $fh $runblock if $runblock;
}

=back

=cut

1;
