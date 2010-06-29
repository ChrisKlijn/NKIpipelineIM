#!/usr/bin/perl
#===============================================================================
#
#         FILE:  db2fastq.pl
#
#        USAGE:  ./db2fastq.pl  
#
#  DESCRIPTION:  
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  C Klijn, c.klijn@nki.nl
#      COMPANY:  NKI
#      VERSION:  0.1
#      CREATED:  06/29/2010
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;
use DBI;

# Run parameters
# 1 - database name (default insert.db)

my $dbf = shift || "insert.db";

# Connect to database

my $dbh = DBI->connect("dbi:SQLite:dbname=$dbf","","", { RaiseError => 1, AutoCommit => 0});

# Design and run query
my $q1 = "select * from read";
my $h1 = $dbh->prepare($q1);
$h1->execute();

while (my @row = $h1->fetchrow_array) {
	my ($readID, $runbc, $seq) = @row[1,2,4];
	print ">", join(":", $readID, $runbc) , "\n", $seq, "\n";
}

