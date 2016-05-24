package Bio::SUPERSMART;

=head1 NAME

Bio::SUPERSMART - self-updating platform for estimation of rates of speciation, migration, 
and relationships of taxa

=head1 SYNOPSIS

 # get help message and quit
 $ smrt
 $ smrt -h
 $ smrt [subcommand] -h
 

=head1 DESCRIPTION

The SUPERSMART pipeline is operated by the command line tool C<smrt>. Consult its help
messages to learn how to use the pipeline. The package whose documentation you are reading 
now is only of interest to developers.

=cut

use version; our $VERSION = version->declare('v0.1.35');

1;
