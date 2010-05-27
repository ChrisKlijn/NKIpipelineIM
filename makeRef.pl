#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  makeRef.pl
#
#        USAGE:  ./makeRef.pl  
#
#  DESCRIPTION:  
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  D Sie (DLS), d.sie@nki.nl
#      COMPANY:  NKI/CMF
#      VERSION:  1.0
#      CREATED:  05/14/2009 02:02:35 PM
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;
use Bio::SeqIO;
use Bio::SearchIO;
use DBI;

$| = 1;

my $dbh = DBI->connect("dbi:SQLite:dbname=insert.db","","", { RaiseError => 1, AutoCommit => 0});
my $base = shift || "out";
my $run_name = shift || $base;

my $q1 = "select qStart, qEnd, sequence, ris.seq_id, iName, iStart, iEnd, iOri from ris left join mapping using (mappingID) left join read using (readID) where run_name = ?";
my $h1 = $dbh->prepare($q1);

$h1->execute($run_name);
open REF, ">$base.ref.fa";
while (my @res = $h1->fetchrow_array) {
	my $l = $res[1] - $res[0];

	# Take only the mapped part of the sequence as reference
	########################################################
	my $seq = substr($res[2], $res[0], $l);
	print REF ">", $res[3], "-", join(":", @res[4 .. 7]), "\n$seq\n";
}
close REF;

open QUE, ">$base.que.fa";

my $q2 = "select read.seq_id, sequence from read left join ris using (readID) where risID is null and run_name = ?";
my $h2 = $dbh->prepare($q2);
$h2->execute($run_name);
while (my ($n, $s) = $h2->fetchrow_array) {
	print QUE ">$n\n$s\n";
}

close QUE;
