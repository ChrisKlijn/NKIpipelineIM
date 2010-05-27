#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  export.pl
#
#        USAGE:  ./export.pl  
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
#      CREATED:  05/20/2009 04:18:07 PM
#     REVISION:  ---
#===============================================================================


#BEGIN {
	#my $ROOT = "/net/ensembl-api/ensembl-49";
	#unshift @INC, $ROOT . "/modules";
	#unshift @INC, $ROOT . "/ensembl/modules";
	#unshift(@INC, $ROOT . "/bioperl-live");
#}

use strict;
use warnings;
use DBI;
#use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp expand);
#use Bio::EnsEMBL::Registry;
#Bio::EnsEMBL::Registry->load_registry_from_db(
		#-host   =>      "legion.nki.nl",
		#-port   =>      3306,
		#-user   =>      "ensro",
		#-pass   =>       "ensro",
		#-verbose => 0,
		#);

=head
	Index to column mapping:

	0:risID
	1:mappingID
	2:readID
	3:seq_id
	4:iName
	5:iStart
	6:iEnd
	7:iOri
	8:flag
	9:seq_id
	10:run_name
	11:barcode
	12:sequence
	13:length2vec
	14:position
	15:flag

=cut

#######################################
# Make your sorting modification here #
#######################################

# On which columns to sort:
###########################

sub sortfun {
	$a->[11] cmp $b->[11] || 	#bc
	$a->[4] cmp $b->[4] ||		#chromosome
	$a->[7] cmp $b->[7] ||		#orientation
	$a->[13] <=> $b->[13] ||	#length
	$a->[15] <=> $b->[15] ||	#flag 
	$a->[14] <=> $b->[14]		#position 
}

# On which columns to cluster:
##############################

my @grpCols = (
	11,		#bc
	4,		#chromosome
	7,		#orientation
	13,		#length
	15,		#flag
#	14,		#position overrides cw
);

#############################################
# No modifications beyond this point please #
#############################################

my $run_name = shift || "out1";
my $dbf = shift || "insert.db";
my $squat = shift || 0;
my $window = shift || 10000;
my $cw = shift || 50; # clustering window
my $spl = shift || "T7";

$| = 1;

#my $SA = Bio::EnsEMBL::Registry->get_adaptor("mus_musculus", "core", "Slice");

my %genes = ();
my $dbh = DBI->connect("dbi:SQLite:dbname=$dbf","","", { RaiseError => 1, AutoCommit => 0});

my $q1 = "select * from ris left join read using (readID) where run_name = ?";
my $h1 = $dbh->prepare($q1);
$h1->execute($run_name);

my %splink = ();
my $q2 = "select * from vulgar left join read using (readID) where run_name = ? and tName = ?";
$splink{$_->[2]} = [@$_] foreach @{$dbh->selectall_arrayref($q2, {}, $run_name, $spl)};

my %vec = ();
my $q3 = "select * from vulgar left join read using (readID) where run_name = ? and (tName = ? or tName = ?)";
$vec{$_->[2]} = [@$_] foreach @{$dbh->selectall_arrayref($q3, {}, $run_name, "SB", "MMTV7")};

print STDERR "Nr.spl: " . scalar(keys(%splink)), "\n";
print STDERR "Nr.vec: " . scalar(keys(%vec)), "\n";

my @data = ();
while (my @row = $h1->fetchrow_array) {
	my ($readID, $c, $s, $e, $o, $bc, $seq) = @row[3,4,5,6,7,11,12];
	my $l = 0;
	my $extrp = 0;
	my $splinkFlag = 0;

#	toggle splink exclusion
#	next unless (exists($splink{$readID}));
#	Determine length:
#	- with splink, extrapolated junction
#	- without splink, blast alignment length
############################################

	# Determine length
	##################

	if (exists($splink{$readID})) {
		$extrp = $splink{$readID}->[7] - 1; # Extrapolate the beginning of the splink
		substr($seq, $splink{$readID}->[3], length($seq), "");
		$splinkFlag = 1;
		if (exists($vec{$readID})) {
			$l = length($seq) - $vec{$readID}->[4];
		}
		else {
			$splinkFlag += 2;
			#die "Oh noo! $readID";	
		}
		
	}
	else {
		# tEnd - tStart
		$l = $e - $s;
	}
	
	# To squat or not to squat
	##########################

	if ($squat == 1) {
		my $nseq = substr($seq, $s);
		(my $squated = $nseq) =~ tr/[AGCT]/[AGCT]/s;
		$l = length($squated) - $extrp;
	}
	else {
		$l -= $extrp;	
	}

	# Detemine ris position
	#######################
	my $p = ($o eq "+")?$s:$e;

	push @data, [@row, $l, $p, $splinkFlag];
}

my $total = 0;


# Sorting and clustering here
#####################################
my @sdata = sort sortfun @data;

my $skipStr = join(":", @{$sdata[0]}[@grpCols]);
my @cache = ();
my $clust = -1;
foreach my $o (0 .. $#sdata) {
	my @row = @{$sdata[$o]};
	
	if ($skipStr ne join(":", @row[@grpCols])) {	
		$clust++;
		$skipStr = join(":", @row[@grpCols]);
		push @{$cache[$clust]}, $sdata[$o];
		next;
	}
	
	my $p = $row[14];
	my $np = $sdata[$o - 1][14];

	if (abs($np - $p) > $cw) {
		$clust++; 
		$skipStr = join(":", @row[@grpCols]);
		push @{$cache[$clust]}, $sdata[$o];
		next;
	}
	push @{$cache[$clust]}, $sdata[$o];
}

print join("\t", qw/chr start end ori length index depth flag length2vec rReadID/), "\n";
foreach my $c (@cache) {
	
	# determine most occuring start end combination 
	my %hsh = ();
	foreach my $ref (@{$c}) {
		my @row = @$ref;
		my $se = join(":", @row[5,6]);
		push @{$hsh{$se}}, $row[3];
	}

	my $se = (sort { scalar(@{$hsh{$b}}) <=> scalar(@{$hsh{$a}}) } keys(%hsh))[0];
	my ($s, $e) = split(/\:/, $se);
	my $readID = (sort @{$hsh{$se}})[0];
	
	my $depth = scalar(@$c);

	my @prow = @{$c->[0]};
	# my ($pre, $post) = &proxGenes($prow[4], $s, $e,$window);
	# print join("\t", $prow[4], $s, $e, $prow[7], $e - $s, $prow[11], $depth, @prow[15, 13], $readID, join("|", @$pre), join("|", @$post)) , "\n";
	print join("\t", $prow[4], $s, $e, $prow[7], $e - $s, $prow[11], $depth, @prow[15, 13], $readID) , "\n";
}

print STDERR scalar(@cache), "\n";

=head
sub proxGenes {
	my $c = shift;
	my $s = shift;
	my $e = shift;
	my $window = shift;
	
	my $slicePre = $SA->fetch_by_region("chromosome", $c, $s - $window, $s);
	my $slicePost = $SA->fetch_by_region("chromosome", $c, $e, $e + $window);

	my @pre = ();
	my $seqPre = $slicePre->seq;
	my $seqPost = $slicePost->seq;
	reverse_comp(\$seqPre);

	my %rePre = ();
	my %rePost = ();
	my $p = 0;

	my $genesPre = $slicePre->get_all_Genes();
	foreach my $g (@$genesPre) {
		my $id = $g->external_name || $g->display_id;
		push @pre, $id;
	}
	my @post= ();
	my $genesPost = $slicePost->get_all_Genes();
	foreach my $g (@$genesPost) {
		my $id = $g->external_name || $g->display_id;
		push @post, $id;
	}
	return (\@pre, \@post);
}
=cut
