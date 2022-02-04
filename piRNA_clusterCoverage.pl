#!/usr/bin/perl
use strict;
use warnings;

die unless @ARGV == 1;
my ($file) = @ARGV;

#$file variable stores file containing piRNA read counts for each chromosomal location within the piRNA cluster boundary.

chomp $file;
my $temp = $file;
$temp =~ /^21/; 


my ($strain) = $temp =~ /^(.*)_cluster/;	
	
# Following script will sum up the piRNA read counts of each chromosomal location within the piRNA cluster boundary.

my %piRNA_cluster;
my %TotalCount;
open IN, '<', $file or die;
while (my $line = <IN>) {
	chomp $line;
	my ($cluster, $i, $strand, $coverage) = (split /\t/, $line) [0,1,2,3];
	
	${$piRNA_cluster{$cluster}}{$strand} += $coverage;
	$TotalCount{$cluster} += $coverage;
	
	}
close IN;

my $output = $strain."_piRNAcluster_strandBaisedDensity.txt";
open OUT, '>', $output or die;

foreach my $key (sort keys %piRNA_cluster){
	foreach my $strand (keys %{$piRNA_cluster{$key}}) {
		print OUT "$key\t$strand\t${$piRNA_cluster{$key}}{$strand}\n";
	}
	

}

close OUT;




my $ot = $strain."_piRNAcluster_TotalDensity.txt";
open OUT, '>', $ot or die;

foreach my $cluster (sort keys %TotalCount){
	
	print OUT "$cluster\t$TotalCount{$cluster}\n";

}
close OUT;