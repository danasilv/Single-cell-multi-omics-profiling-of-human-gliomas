use strict;
use IO::Handle; STDOUT->autoflush;

my $trimExtraR1PreAdapter = 3; #trim 3bp (2 filled in bp + "A") before adapter location
my $trimExtraR2PreAdapter = 1; #trim 1bp ("A") before adapter location
my $postTrimMinReadLength = 36;
my $adapterMinMatchPct = .8;
my $avQualityCutoff = 48;

my @barcodes  = qw(
ACAACC
ACTCAC
AGGATG
ATCGAC
CAAGAG
CATGAC
CCTTCG
CGGTAG
CTATTG
CTCAGC
GCATTC
GTGAGG
GTTGAG
TATCTC
TCTCTG
TGACAG
TGCTGC
ACAGAC
AGAAGG
ATCAAG
CCATAG
GAAGTC
GGTAAC
TGTAGG
);
my %validBarcodes = ();
foreach my $barcode (@barcodes)
{
	$validBarcodes{$barcode} = 1;
}

my $infile1 = shift;
my $infile2 = shift;
my $outfile = "$infile1.counts";
my $hiOutfile = "$infile1.hiCounts";
my $summOutfile = "$infile1.summ";

my $usage = "usage: splitFastqPair.pl {file1} {file2}\n";

die $usage if $infile1 eq "";
die $usage if $infile2 eq "";
die "file1 '$infile1' eq file2 '$infile2'\n$usage" if $infile1 eq $infile2;

open IN1, $infile1 or die "$!\n$usage";
open IN2, $infile2 or die "$!\n$usage";

my %barcodeCounts = ();
my %printedCounts = ();
my $total = 0;
my $validBarcodeCount = 0;
my %files1 = ();
my %files2 = ();
my $matchSecondOnly = 0;

my $adapterR1Count = 0;#read 1 count in which adapter was found and excised
my $adapterR2Count = 0;
my $adapterDimerCount = 0;
my $lowQtrim1Count = 0;#read 1 count in which low-quality bases were trimmed
my $lowQtrim2Count = 0;
my $printedReads = 0; #printed after trimming exceeds min length

my $idLine1;
my $fastqLine1;
my $strandLine1;
my $qual1;

my $idLine2;
my $fastqLine2;
my $strandLine2;
my $qual2;
print "processed 0";
while ($idLine1 = <IN1>)
{
	$total++;
	$fastqLine1 = <IN1>;
	chomp $fastqLine1;
	$strandLine1 = <IN1>;
	$qual1 = <IN1>;
	chomp $qual1;
	
	$idLine2 = <IN2>;
	$fastqLine2 = <IN2>;
	chomp $fastqLine2;
	$strandLine2 = <IN2>;
	$qual2 = <IN2>;
	chomp $qual2;

	if ($total %100000 == 0)
	{
		print "\rprocessed $total";
		my $id1 = (split " ", $idLine1)[0];
		my $id2 = (split " ", $idLine2)[0];
		die "ERROR: IDs are off\nid1: $id1\nid2: $id2\n" unless $id1 eq $id2; 
	}

	my $barcode = substr($fastqLine1, 0,6);
	my $barcode2 = substr($fastqLine2,0,6);

	my $letter = substr($fastqLine1,6,1);
	my $letter2 = substr($fastqLine2,6,1);
	#print $fastqline;
	#print "letter: $letter barcode: $barcode\n";
	if ($letter eq "T" && $validBarcodes{$barcode})
	{
		printBarcode($barcode);
		$barcodeCounts{$barcode}++;	
	}
	elsif ($letter2 eq "T" && $validBarcodes{$barcode2})
	{
		$matchSecondOnly++;
		printBarcode($barcode2);
		$barcodeCounts{$barcode2}++;
	}
	else
	{
		$barcodeCounts{$barcode}++;
	}
}

sub printBarcode
{
	my $printBarcode = shift;
	my $trimmedFastq1 = substr($fastqLine1,7);
	my $trimmedQual1 = substr($qual1,7);
	my $trimmedFastq2 = substr($fastqLine2,9);
	my $trimmedQual2 = substr($qual2,9);

	my $exactMatchNum = length($printBarcode) + 0;

	#trim adapter from read 1
	my $revcompBarcode = reverse($printBarcode);
	$revcompBarcode =~ tr/ACGTacgt/TGCAtgca/;
	my $adapter1 = $revcompBarcode."AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
	my $adapter1Len = length($adapter1);

	my $match1 = substr($adapter1,0,$exactMatchNum);
	my $match1Loc = rindex($trimmedFastq1,$match1);
	while ($match1Loc > -1) #frequent double adapters.
	{
		my $trimmed1Len = length($trimmedFastq1);
		my @read1Arr = split("", $trimmedFastq1);
		my @adapter1Arr = split "", $adapter1;
		my $read1MatchCount = $exactMatchNum;
		my $read1MismatchCount = 0;
		my $read1SlackVar = 0; #for times when the sequencer stutters and a bp is read twice
		#print "matching for adapter1 $match1\n";
		#print "lastmatchloc: $match1Loc trimmedLen: $trimmedLen\n";
		for (my $adapter1Iter = $exactMatchNum; $adapter1Iter < $trimmed1Len; $adapter1Iter++)
		{
			my $readIter = $adapter1Iter + $match1Loc + $read1SlackVar;
			last unless ($readIter < $trimmed1Len);
			last unless ($adapter1Iter < $adapter1Len);
			#print "testing readArr[$readIter]: $readArr[$readIter] vs adapter1Arr[$adapter1Iter]:$adapter1Arr[$adapter1Iter]\n";
			if ($read1Arr[$readIter] eq $adapter1Arr[$adapter1Iter])
			{
				$read1MatchCount++;
			}
			else
			{
				if ($read1Arr[$readIter + 1] eq $adapter1Arr[$adapter1Iter])
				{
					$read1SlackVar++;
					$read1MatchCount++;
				}
				else
				{	
					$read1MismatchCount++;
				}
			}
		}

		#make sure match pct is above threshold and no more than 2bp have been stuttered
		my $match1Pct = ($read1MismatchCount > 0)?$read1MatchCount/$read1MismatchCount:1;
		if ($match1Pct > $adapterMinMatchPct && $read1SlackVar < 3)
		{
			#make sure read is long enough
			my $newFastq1Len = $match1Loc - $trimExtraR1PreAdapter;
			if ($newFastq1Len < 0)
			{
				$newFastq1Len = 0;
				$adapterDimerCount++;
			}
			$trimmedFastq1 = substr($trimmedFastq1,0,$newFastq1Len);
			$trimmedQual1 = substr($trimmedQual1,0,$newFastq1Len);
			$adapterR1Count++;

			#update match loc for cases of double adapters
			$match1Loc = rindex($trimmedFastq1,$match1);
		}
		else
		{
			$match1Loc = -1;
		}
	}

	#trim adapter from read 2
	my $revcompBarcode = reverse($printBarcode);
	$revcompBarcode =~ tr/ACGTacgt/TGCAtgca/;
	my $adapter2 = $revcompBarcode."AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA";
	my $adapter2Len = length($adapter2);

	my $match2 = substr($adapter2,0,$exactMatchNum);
	my $match2Loc = rindex($trimmedFastq2,$match2);
	while ($match2Loc > -1)
	{
		my $trimmed2Len = length($trimmedFastq2);
		my @read2Arr = split("", $trimmedFastq2);
		my @adapter2Arr = split "", $adapter2;
		my $read2MatchCount = $exactMatchNum;
		my $read2MismatchCount = 0;
		my $read2SlackVar = 0; #for times when the sequencer stutters and a bp is read twice
		#print "matching for adapter2 $match2\n";
		#print "lastmatchloc: $match2Loc trimmed2Len: $trimmed2Len\n";
		for (my $adapter2Iter = $exactMatchNum; $adapter2Iter < $trimmed2Len; $adapter2Iter++)
		{
			my $readIter = $adapter2Iter + $match2Loc + $read2SlackVar;
			last unless ($readIter < $trimmed2Len);
			last unless ($adapter2Iter < $adapter2Len);
			#print "testing readArr[$readIter]: $readArr[$readIter] vs adapter2Arr[$adapter2Iter]:$adapter2Arr[$adapter2Iter]\n";
			if ($read2Arr[$readIter] eq $adapter2Arr[$adapter2Iter])
			{
				$read2MatchCount++;
			}
			else
			{
				if ($read2Arr[$readIter + 1] eq $adapter2Arr[$adapter2Iter])
				{
					$read2SlackVar++;
					$read2MatchCount++;
				}
				else
				{	
					$read2MismatchCount++;
				}
			}
		}

		my $match2Pct = ($read2MismatchCount > 0)?$read2MatchCount/$read2MismatchCount:1;
		if ($match2Pct > $adapterMinMatchPct && $read2SlackVar < 3)
		{
			my $newFastq2Len = $match2Loc - $trimExtraR2PreAdapter;
			$newFastq2Len = 0 if ($newFastq2Len < 0);
			$trimmedFastq2 = substr($trimmedFastq2,0,$newFastq2Len);
			$trimmedQual2 = substr($trimmedQual2,0,$newFastq2Len);
			$adapterR2Count++;

			#update match loc for cases of double adapters
			$match2Loc = rindex($trimmedFastq2,$match2);
		}
		#if last instance of adapter didn't match, exit loop
		else
		{
			$match2Loc = -1;
		}
	}


	#trim low-quality bases from 3' end of read 1
	#4bp sliding window from 3' end until quality is greater than 15
	my @qual1Arr = split("", $trimmedQual1);
	my $lastGoodLoc = 0;
	my $trimmedQual1Len = length($trimmedQual1);
	for (my $i = $trimmedQual1Len - 1; $i >= 3; $i--)
	{
		my $av = (ord($qual1Arr[$i]) + ord($qual1Arr[$i-1]) + ord($qual1Arr[$i-2]) + ord($qual1Arr[$i-3])) / 4;
		# Low quality cutoff = 15 + 33 = 48
		if ($av > $avQualityCutoff)
		{
			$lastGoodLoc = $i + 1;
			last;
		}
	}
#	print "lastGoodLoc: $lastGoodLoc (len $trimmedQual1Len)\n";

	if ($lastGoodLoc < $trimmedQual1Len)
	{
		$lowQtrim1Count++;
		$trimmedFastq1 = substr($trimmedFastq1,0,$lastGoodLoc);
		$trimmedQual1 = substr($trimmedQual1,0,$lastGoodLoc);
	}
#	print "Post$idLine1\nfastq: $trimmedFastq1\nqual:  $trimmedQual1\n";

	#trim low-quality bases from 3' end of read 2
	#4bp sliding window from 3' end until quality is greater than 15
	my @qual2Arr = split("", $trimmedQual2);
	my $lastGoodLoc = 0;
	my $trimmedQual2Len = length($trimmedQual2);
	for (my $i = $trimmedQual2Len - 1; $i >= 3; $i--)
	{
		my $av = (ord($qual2Arr[$i]) + ord($qual2Arr[$i-1]) + ord($qual2Arr[$i-2]) + ord($qual2Arr[$i-3])) / 4;
		# Low quality cutoff = 15 + 33 = 48
		if ($av > $avQualityCutoff)
		{
			$lastGoodLoc = $i + 1;
			last;
		}
	}

	if ($lastGoodLoc < $trimmedQual2Len)
	{
		$lowQtrim2Count++;
		$trimmedFastq2 = substr($trimmedFastq2,0,$lastGoodLoc);
		$trimmedQual2 = substr($trimmedQual2,0,$lastGoodLoc);
	}

	#print if lengths are long enough
	if (length($trimmedFastq1) > $postTrimMinReadLength && length($trimmedFastq2) > $postTrimMinReadLength)
	{
		if (!exists($files1{$printBarcode}))
		{
			my $fh1;
			open $fh1, ">$infile1.$printBarcode.fastq" or die $!;
			$files1{$printBarcode} = $fh1;
			my $fh2;
			open $fh2, ">$infile2.$printBarcode.fastq" or die $!;
			$files2{$printBarcode} = $fh2;
		}
		my $fh1 = $files1{$printBarcode};
		print $fh1 "$idLine1$trimmedFastq1\n$strandLine1$trimmedQual1\n";

		my $fh2 = $files2{$printBarcode};
		print $fh2 "$idLine2$trimmedFastq2\n$strandLine2$trimmedQual2\n";

		$printedCounts{$printBarcode}++;
		$printedReads++;
	}

	$validBarcodeCount++;

}


open ALL, ">$outfile";
open HI, ">$hiOutfile";
my $allCatCount = 0;
my $hiCatCount = 0;
my $allCount = 0;
my $hiCount = 0;
foreach my $key (sort keys %barcodeCounts)
{
	my $count = $barcodeCounts{$key};
	$allCatCount++;
	$allCount+= $count;
#	print "$key :$count\n";
	print ALL "$key: $count\n";
	if ($barcodeCounts{$key} > 1000)
	{
		$hiCatCount++;
		$hiCount += $count;
		print HI "$key: $count\n";
	}
}

open INFO, ">$summOutfile";
print INFO "validBarcodes\ttotalReads\n";
print INFO "$validBarcodeCount\t$total\n";
print INFO "\nTotal reads read: $total\n";
print INFO "Reads with any barcode: $allCount ($allCatCount unique barcodes)\n";
my $validBarcodeArrCount = @barcodes;
print INFO "Reads with valid barcodes: $validBarcodeCount ($validBarcodeArrCount)\n";
print "\nprinted \n$allCount/$total ($allCatCount unique barcodes) hi(>1000) barcodes to $outfile and \n$hiCount/$total ($hiCatCount unique barcodes) to $hiOutfile\n";
my $validPct = $validBarcodeCount/$total;
print "got $validBarcodeCount/$total valid barcodes ($validPct)\n";
print "$matchSecondOnly reads matched second barcode only\n";
print INFO "Valid Percent: $validBarcodeCount/$total = $validPct\n";
print INFO "Reads with barcode in second read only: $matchSecondOnly\n";
print INFO "--ADAPTER EXCISING--\n";
print INFO "AdapterMinMatchPct: $adapterMinMatchPct\n";
print INFO "Reads with run-over adapter in read 1: $adapterR1Count\n";
print INFO "Reads with run-over adapter in read 2: $adapterR2Count\n";
print INFO "Adapter Dimer Count: $adapterDimerCount\n";
print INFO "--3' QUALITY TRIMMING--\n";
print INFO "avQualityCutoff: $avQualityCutoff (4bp sliding window average)\n";
print INFO "Reads trimmed for low 3' quality in read1: $lowQtrim1Count\n";
print INFO "Reads trimmed for low 3' quality in read2: $lowQtrim2Count\n";
print INFO "--TOTAL PRINTED READS--\n";
print INFO "Total printed reads: $printedReads (longer than $postTrimMinReadLength"."bp after adapter excising and trimming)\n";
print INFO "Barcode\tprinted reads\n";
foreach my $barcode (@barcodes)
{
	my $count = $printedCounts{$barcode} || 0;
	print INFO "$barcode\t$count\n";
}

