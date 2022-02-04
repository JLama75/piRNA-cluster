#!/usr/bin/perl
use strict;
use warnings;

die unless @ARGV == 2;
my ($file, $cluster_loci) = @ARGV;

#$file variable stores file containing piRNA reads uniquely mapped to reference sequence.
#$cluster_loci stores file containing chromosomal positions of annotated piRNA clusters

#multiple invocation arguments <>
chomp $file;
my $temp = $file;
$temp =~ /^21/;
my ($strain) = $temp =~ /^(.*)\.piRNA/;	# extracts and stores the file name (eg, $strain = 21147-R1 or 21147-R2 and so on.)
	

#Following script counts the total reads mapped to each chromosomal position.
# piRNA reads with common sequence and particular strand bias (+/-) are counted and stored in hash named reads. 
# Also storing the unique chromosomal co-ordinates and strand information in hash: f1
# eg. f1{'chr2L'}{'1'} {+} = n.  
my %reads;
my %f1;

open IN, '<', $file or die;
while (my $line = <IN>) {
	chomp $line;
    my ($chrom, $start, $seq, $strand) = (split /\t/, $line)[1, 2, 3, 4]; #here $chrom is the chromosomal arm; $start --> chromosomal co-ordinate; $seq --> nucleotide sequence; $straind --> +/-
	if (exists ${$reads{$seq}}{$strand}) {
		${$reads{$seq}}{$strand}++;	
		${${$f1{$chrom}}{$start}}{$strand}++;
	}
	else {
		${$reads{$seq}}{$strand} =1;
		${${$f1{$chrom}}{$start}}{$strand} = 1;
		
		
	}
}
close IN;


#---------------------------------------------------------------------------------------------------------------------

#opening output file with name from $strain + "...". (eg, 21147-R1.piRNA_cluster_expression.txt)
my $output = $strain."_ReadCountCoverage.txt"; 
open OUT, '>', $output or die;


#Now reading the piRNA cluster annotation file and extracting the piRNA cluster boundary information
#extracting and listing out only the piRNA read counts mapping to chromosomal positions within piRNA cluster boundary.

open CLS, '<', $cluster_loci or die;

while (my $line = <CLS>) {
    chomp $line;
	my ($chrom, $start, $end) = (split /\t/, $line) [1,2,3];
	my $key = $chrom.":".$start."..".$end;
	#using loop to print out read count from piRNA cluster boundary (eg- from "start = 1" to "stop =7020ch" cluster boundaries)
	foreach my $i ($start..$end) {
		if (exists ${${$f1{$chrom}}{$i}}{'+'}) {
			print OUT "$key\t$i\t\+\t${${$f1{$chrom}}{$i}}{'+'}\n";
		}
		else {
			print OUT "$key\t$i\t\+\t0\n";
		}

		if (exists ${${$f1{$chrom}}{$i}}{'-'}) {
            print OUT "$key\t$i\t\-\t${${$f1{$chrom}}{$i}}{'-'}\n";
        }
        else {
            print OUT "$key\t$i\t\-\t0\n";
        }
	}
}
close CLS;
close OUT;

print "The end!\n\n";











