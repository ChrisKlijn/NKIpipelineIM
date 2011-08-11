#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  process.pl
#
#        USAGE:  ./process.pl 
#
#  DESCRIPTION:  
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:   (), <>
#      COMPANY:  
#      VERSION:  1.0
#      CREATED:  09/11/2008 05:55:08 PM CEST
#     REVISION:  ---
#===============================================================================
use strict;
use warnings;
use Bio::SeqIO;
use Bio::SearchIO;
use DBI;

sub filter {
	# Define filter here, return a zero to skip
	###########################################

=head
	index definitions:

	$hits is based on PLS output:

	$hits->[0] => match 
	$hits->[1] => mismatch 
	$hits->[2] => repmatch 
	$hits->[3] => N's 
	$hits->[4] => Q gap count 
	$hits->[5] => Q gap bases
	$hits->[6] => T gap count
	$hits->[7] => T gap bases
	$hits->[8] => strand
	$hits->[9] => Q name
	$hits->[10] => Q size 
	$hits->[11] => Q start
	$hits->[12] => Q end
	$hits->[13] => T name 
	$hits->[14] => T size
	$hits->[15] => T start
	$hits->[16] => T end 
	$hits->[17] => block count 
	$hits->[18] => block sizes
	$hits->[19] => Q starts
	$hits->[20] => T starts

	$vector is based on VULGAR output:

	$vector->[0] => vulgar:
	$vector->[1] => Q name (readID)
	$vector->[2] => Q start
	$vector->[3] => Q end
	$vector->[4] => Q ori
	$vector->[5] => T name (vector)
	$vector->[6] => T start 
	$vector->[7] => T end
	$vector->[8] => T ori
	$vector->[9] => Raw score
	$vector->[10 .. $#$vector] => Variable number of triplets
=cut

	my $hits = shift;
	my $vecs = shift;

	###############
	# > 95% match #
	###############
	
	return(0) if (($hits->[0] / ($hits->[0] + $hits->[1])) < 0.95);

	#################
	# Vector lineup #
	#################

	my $vector = $vecs->{'TRNSP5LTR'} || undef;

	return(0) unless $vector;
	return(0) if $vector->[6] != 0; # must start at begining
	return(0) if $vector->[7] < 23; # must be longer than 23 nt long
	return(0) if (abs($hits->[11] - $vector->[3]) > 5); # distance from hit must be shorter than 5
	return(1);
}

$| = 1;

my $dbh = DBI->connect("dbi:SQLite:dbname=insert.db","","", { RaiseError => 1, AutoCommit => 0});
my $fa = shift;
my $base = shift;
my $psl = $base.".psl";
my $exo = $base.".exo";
my $bcexo = $base.".bc.exo";
my $report = $base.".report.txt";
my $prebc = shift || "";
my $run_name = shift || $prebc;
my $mkdb = shift || 1;
my $minCnt = shift || 150;

if ($mkdb == 1) {
my $dbdef = [qq/
	
	create table read (
		readID integer primary key autoincrement,
		seq_id varchar(64),
		run_name varchar(64),
		barcode varchar(8),
		sequence text not null
	);
/,qq/
	create table mapping (
		mappingID integer primary key autoincrement,
		readID int(11) not null,
		match int(11) not null,
		misMatch int(11) not null,
		repMatch int(11) not null,
		N int(11) not null,
		qGapCount int(11) not null,
		qGapBases int(11) not null,
		tGapCount int(11) not null,
		tGapBases int(11) not null,
		strand varchar(4) not null,
		qName varchar(64) not null,
		qSize int(11) not null,
		qStart int(11) not null,
		qEnd int(11) not null,
		tName varchar(64) not null,
		tSize int(11) not null,
		tStart int(11) not null,
		tEnd int(11) not null,
		blockCount int(11) not null,
		blockSizes varchar(64) not null,
		qStarts varchar(128) not null,
		tStarts varchar(128) not null
	);
/,qq/
	create table vulgar (
		vulgarID integer primary key autoincrement,
		readID int(11) not null,
		qName varchar(64) not null,
		qStart int(11) not null,
		qEnd int(11) not null,
		qOri varchar(4) not null,
		tName varchar(64) not null,
		tStart int(11) not null,
		tEnd int(11) not null,
		tOri varchar(4) not null,
		score int(11) not null,
		trip varchar(128)
	);
/,qq/
	create table ris (
		risID integer primary key autoincrement,
		mappingID int(11) not null,
		readID int(11) not null,
		seq_id varchar(64),
		iName varchar(64),
		iStart int(11),
		iEnd int(11),
		iOri int(11),
		flag int(11)
	);

/];

$dbh->do($_) foreach @$dbdef;
$dbh->commit;
}

#################
# Read barcodes #
#################
print STDERR "Loading barcodes\n";
my %r2bc = ();
open BCEXO, "<$bcexo";
while (<BCEXO>) {
	next unless (/^vulgar:/);
	chomp;
	my @row = split(/\s+/, $_);
	my $i = 0;
	my $readID = $row[1];
	my $bc = $prebc . $row[5];

	# If one read has mulitiple barcodes, select highest score
	##########################################################
	if (exists($r2bc{$readID})) {
		next if ($row[9] < $r2bc{$readID}->[9]);
	}
	$r2bc{$readID} = [@row];
}
close BCEXO;

my %bc = ();


##############
# Load reads #
##############

my $in = Bio::SeqIO->new(-file=>$fa, -format => "fasta");

my $cnt = 0;
my $q1 = "insert into read (seq_id, barcode, sequence, run_name) values (?, ?, ?, ?)";
my $h1 = $dbh->prepare($q1);

print STDERR "Loading reads\n";
while (my $seq = $in->next_seq()) { 
	my $s = $seq->seq();
	my $id = $seq->primary_id();
	my $code = "NA";
	$code = $prebc . $code;

	if (exists($r2bc{$id})) {
		$code = $prebc . $r2bc{$id}->[5];
	}

	$bc{$code}{'seq'}{$seq->primary_id} = $s;
	$h1->execute($id, $code, $s, $run_name);
	print STDERR "\r$cnt" if ((++$cnt % 1000) == 0);
}
print STDERR "\rTotal reads: $cnt\n";
$dbh->commit;

my %read = @{$dbh->selectcol_arrayref("select readID, seq_id from read", {Columns => [2,1]} )};

print STDERR "Total barcodes: " . scalar(keys(%bc)) . "\n";

my $total = 0;
my %valid = ();

########################################
# Discard reads for defective barcodes #
########################################

foreach my $k (keys(%bc)) { 
	my @seqs = keys(%{$bc{$k}{'seq'}});
	if (scalar(@seqs) >= $minCnt) {
		$valid{$_} = $k foreach @seqs;
	}
	else {
		my $n = scalar(@seqs);
		$total += $n;
		print STDERR "$k with $n reads discarded\n";
		delete($bc{$k});
	}
}

print STDERR "$total reads discarded due to under representation of barcode`\n";


#################
# Read mappings #
#################

my %mappings = ();
open PSL, "<$psl";
while (<PSL>) {
	last if (/----/);
}

my $q2 = "insert into mapping ( readID,
match, misMatch, repMatch, N, qGapCount, qGapBases, tGapCount, tGapBases, strand, qName, 
qSize, qStart, qEnd, tName, tSize, tStart, tEnd, blockCount, blockSizes, qStarts, tStarts
) values (". join(",", ("?") x 22) .")";
my $h2 = $dbh->prepare($q2);

$cnt = 1;
print STDERR "Loading mappings\n";
while (<PSL>) {
	chomp;
	my @row = split(/\t/, $_);
	my $i = 0;
	$h2->execute($read{$row[9]}, @row);
	my $mID = $dbh->func('last_insert_rowid');
	push @{$mappings{$row[9]}}, [@row, $mID];
	next if (not exists($valid{$row[9]}));
	$bc{$valid{$row[9]}}{'hits'}{$row[9]} = 0;
	print STDERR "\r$cnt" if ((++$cnt % 1000) == 0);
}
close PSL;
$dbh->commit;
print STDERR "\rTotal mappings: $cnt\n";


################
# Read vectors #
################

my %vec = ();
my $q3 = "insert into vulgar ( readID,	
		qName,
		qStart,
		qEnd,
		qOri,
		tName,
		tStart,
		tEnd,
		tOri,
		score,
		trip
) values (?,?,?,?,?,?,?,?,?,?,?)";
my $h3 = $dbh->prepare($q3);

$cnt = 0;
open EXO, "<$exo";
while (<EXO>) {
	next unless (/^vulgar:/);
	chomp;
	my @row = split(/\s+/, $_);
	$h3->execute($read{$row[1]}, @row[1 .. 9], join(" ", @row[10 .. $#row]));
	print STDERR "\r$cnt" if ((++$cnt % 1000) == 0);
	next if (not exists($valid{$row[1]}));
	$vec{$row[1]}{$row[5]} = [@row];
}
$dbh->commit;
print STDERR "\rTotal vecs: $cnt\n";


####################
# Process all data #
####################

open REP, ">$report";

my $q4 = "insert into ris (mappingID, readID, seq_id, iName, iStart, iEnd, iOri, flag) values (?,?,?,?,?,?,?,?)";
my $h4 = $dbh->prepare($q4);

print STDERR "BC\treads\tmappings\tinserts\tuniqueInserts\n";
foreach my $code (keys(%bc)) {

	print STDERR $code . "\t" . scalar(keys(%{$bc{$code}{'seq'}})) . "\t" . scalar(keys(%{$bc{$code}{'hits'}})) . "\t";
	
	my $seq = $bc{$code}{'seq'};
	my @insert = ();
	my %site = ();
	foreach my $seq_id (keys(%$seq)) {
		my $hits = $mappings{$seq_id} || undef;
		unless ($hits) {
			next; # no hit, no game
		}

		my @putative = ();
		my $putativeScore = undef;

		my %vecs = %{$vec{$seq_id} || {}};

		my $sequence = $seq->{$seq_id};
		my $vecSeq = $sequence;

		##############
		# Print read #
		##############
		
		print REP ">$seq_id\n$sequence\n";
		my %vhsh = ();
		foreach my $v (values(%vecs)) {
			next unless ($v);
			push @{$vhsh{$v->[5]}}, $v;

			##########################
			# Print vector alignment #
			##########################

			my $name = "[--" . $v->[5];
			my $nl = length($name);
			my $ml = $v->[3] - $v->[2];

			if ($ml > ($nl + 2)) {
				$name .= ("-") x ($ml - $nl - 1) . "]";
			}
			else {
				print STDERR "$v\n$name\t$ml";
				substr($name, $ml - 1) = "]";
			}
			substr($vecSeq,$v->[2], $ml, $name);
			
		}
		print REP $vecSeq . "\n";	



		foreach my $h (sort {
			$b->[0] <=> $a->[0] ||
			$a->[1] <=> $b->[1] 
		} @$hits) {
			
			###################################################
			# Break if score is decreasing from putativeScore #
			###################################################

			last if ($putativeScore && $putativeScore != $h->[0]);
			
			###########################
			# Print genomic alignment #
			###########################

			my $al = $sequence;
			my $m = $h;
			my $name = "[--" . join(":", @{$m}[13,15,16,0,1]);
			my $nl = length($name);
			my $ml = $m->[12] - $m->[11];

			if ($ml > ($nl + 2)) {
				$name .= ("-") x ($ml - $nl - 1) . "]";
			}
			else {
				substr($name, $nl - 1) = "]";
			}
			substr($al,$m->[11], $ml, $name);
			print REP $al, "\n";
			next if (&filter($h, \%vecs) == 0);

			my $flag = 0;	
			push @putative, [$h->[-1], $read{$h->[9]}, @$h[9, 13, 15, 16, 8], $flag];
			$putativeScore = $h->[0];
			print STDERR "$h\n";
		}
		
		if (@putative == 1) {
			$h4->execute(@{$putative[0]});	
			push @insert, @putative;
			if ($putative[0][6] eq "+") {
				$site{join(":", @{$putative[0]}[3,4,6,7])}++;
			}
			else {
				$site{join(":", @{$putative[0]}[3,5,6,7])}++;
			}
			print REP "Elected: ", join(":", @{$putative[0]}), "\n\n" ;
		}
		else {
			print REP "Putative inserts: ", scalar(@putative), "\n\n";
		}
	}
	print STDERR scalar(@insert) . "\t"; 
	print STDERR scalar(keys(%site)) . "\n";
}
$dbh->commit;
