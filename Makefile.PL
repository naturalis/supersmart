#!/usr/bin/perl
use strict;
use warnings;
use lib 'lib';
use ExtUtils::MakeMaker;

WriteMakefile(
	'NAME'          => 'Bio::SUPERSMART',
	'ABSTRACT_FROM' => 'lib/Bio/SUPERSMART.pm',
	'VERSION_FROM'  => 'lib/Bio/SUPERSMART.pm',
	'AUTHOR'        => 'http://www.supersmart-project.org/#people',
	'LICENSE'       => 'mit',	
	'PREREQ_PM' => {
		'DBIx::Class'           => 0,
		'DBI'                   => 0,
		'Template'              => 0,
		'Config::Tiny'          => 0,
		'Moose'                 => 0,
		'XML::Twig'             => 0,
		'HTML::Parser'          => 0,
		'JSON'                  => 0,
		'Math::Random'          => 0,
		'App::Cmd'              => 0,
		'String::RewritePrefix' => 0,
		'IO::String'            => 0,
		'Statistics::R'         => 0,
		'Parallel::ForkManager' => 0,
		'List::MoreUtils'       => 0,
		# 'Parallel::MPI::Simple' => 0,
		
		# XXX would be good if there was some way to define a dependency on
		# git repositories such that `cpanm --installdeps .` would be able
		# to pick it up. Doesn't work here, but maybe in the extended
		# meta merge stuff?
		# 'BioPerl'               => 0,
		# 'BioPerl-Run'           => 0,
		# 'Bio::Phylo'            => 0, 
	}
);
