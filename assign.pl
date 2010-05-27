#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  assign.pl
#
#        USAGE:  ./assign.pl  
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
#      CREATED:  05/15/2009 01:00:54 PM
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;
use DBI;

$| = 1;
my $psl = shift || "out1.1.psl";
my $maxR = shift || 100;
my $run_name = shift;
my $dbh = DBI->connect("dbi:SQLite:dbname=insert.db","","", { RaiseError => 1, AutoCommit => 0});


#######################
# Fetch reads from DB #
#######################

my $q2 = "select seq_id, readID, barcode from read where run_name = ?";
my $h2 = $dbh->prepare($q2);
my %dbread = ();
$h2->execute($run_name);
while (my @row = $h2->fetchrow_array) {
	$dbread{$row[0]}{readID} = $row[1];
	$dbread{$row[0]}{barcode} = $row[2];
}
$h2->finish;

###################
# Read PSL output #
###################

my $q1 = "insert into ris (mappingID, readID, seq_id, iName, iStart, iEnd, iOri, flag) values (?,?,?,?,?,?,?,?)";
my $h1 = $dbh->prepare($q1);

open PSL, "<$psl";
while (<PSL>) {last if /-----------/} # Skip header

my %stats = ();
my $read = undef;
my $bc = undef;
while (<PSL>) {
	my @row = split(/\t/, $_);
	$read ||= $row[9];
	$bc = $dbread{$read}{barcode};

	# Process mappings per read
	###########################
	if ($read ne $row[9]) {
		next unless ($read);
		if (not exists($stats{chromosomes})) {
			$read = $row[9];
			%stats = ();	
			next;
		}

		my $skip = 0;
		my @cs = keys(%{$stats{chromosomes}});
		my @os = keys(%{$stats{ori}});
		my @ss = sort { $a <=> $b } @{$stats{starts} || []};
		my $srange = $ss[-1] - $ss[0];
		my @se = sort { $a <=> $b } @{$stats{ends} || []};
		my $erange = $se[-1] - $se[0];


		# Skip if more than 1 chromosome is mapped
		##########################################
		$skip++ if ($#cs != 0); 
		
		# Skip on orientation mismatch
		##############################
		$skip++ if ($#os != 0);

		# Select start/end range
		########################
		my $range = ($os[0] eq "+")?$srange:$erange;
		
		# Skip if target insertions differ too much
		###########################################
		$skip++ if ($range > $maxR);

		
		if ($skip == 0) {
			print join(":", $dbread{$read}{readID}, $read, $cs[0], $ss[0], $se[-1], $os[0], $srange, $erange), "\t";
			my $rv = $h1->execute(1, $dbread{$read}{readID}, $read, $cs[0], $ss[0], $se[-1], $os[0], 0);
			print $rv, " inserted\n";

		}
		else {
			print join(":", $read, $cs[0], $ss[0], $se[-1], $os[0], $srange, $erange, @cs, @os), " skipped\n";	
			
		}

		# Reset for next read
		#####################
		$read = $row[9];
		%stats = ();	
	}

	# Collect info/stats for read
	#############################
	my ($idc, $s, $e, $o) = split(/:/, $row[13]);
	my ($id, $c) = split(/-/, $idc);
	
	next if ($bc ne $dbread{$id}{barcode});

	$stats{chromosomes}{$c} = 0; 
	$stats{ori}{$o} = 0; 
	push @{$stats{starts}}, $s; 
	push @{$stats{ends}}, $e; 

}
close PSL;
$h1->finish;
$dbh->commit;
$dbh->disconnect;

