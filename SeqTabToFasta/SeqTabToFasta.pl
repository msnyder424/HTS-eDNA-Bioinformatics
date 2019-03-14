#usr/bin/perl

#Script to make DADA2 sequence tables into fasta files
#M. R. Snyder 2017

use Cwd;
use File::Basename;

my $Out = "Dada2OTUSummary.txt";
open (OUTSUM, ">", $Out) || die "Can't open $Out: $!\n";
print OUTSUM "Sample\tN no chim reads\tN no chim OTUs\n";

my $Infile = "SeqTab.txt";

open (IN, $Infile) || die "Can't open $Infile: $!\n"; 

my $OTUsDir = basename(getcwd)."Dada2OTUs";
mkdir $OTUsDir unless -d $OTUsDir;

my $c=0;
my @Seqs = ();
while (<IN>){
	chomp ($_);
	my $line = $_;
	my @col = split (/\t/, $line);
	my @Reads = ();
	if ($c == 0){
		for $x (0..$#col){
			my $seq = substr($col[$x], 1, -1);
			push (@Seqs, $seq);
		}
	}
	else {
		my $sample = substr($col[0], 1, -1);
		for $x (1..$#col){
			push (@Reads, $col[$x]);
			#print "$col[$x]\n";
		}
		my $d=1;
		my $Outfile = $OTUsDir."/"."$sample"."OTUs.fasta";
		open (OUT, ">", $Outfile) || die "Can't open $Outfile: $!\n";
		my $totalcount = 0;
	    for $x (0..$#Seqs){
	    	my $CurrCount = $Reads[$x];
	    	$totalcount = $totalcount + $CurrCount;
	    }
	    for $x (0..$#Seqs){
	    	my $CurrLength = length ($Seqs[$x]);
	    	my $CurrCount = $Reads[$x];
	    	if ($CurrCount > 0){
				print OUT ">$d|$Reads[$x]\n$Seqs[$x]\n";
				$d++;
			}
		}
		my $NOTUs = $d-1;
		print OUTSUM "$sample\t$totalcount\t$NOTUs\n";
		close OUT;
	}
	$c++;
}
$nsamples = $c-1;
print "$#Seqs unique sequences in $nsamples samples processed.\n";
close IN;
close OUTSUM;

