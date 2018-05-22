#!usr/bin/perl

use warnings;
use strict;


# set the number of candidate based on the expression level (emPAI)
my $cand = 1000;

# load the candidate
my $ge; my %can; my $count = 0;
open (can1, "emPAI_gene.csv") or die $!;
while (<can1>)     {
	next if /^name/;
	last if $count > $cand;
	chomp;
	my @a = split /,|\t/;
	$can{$a[0]} = $a[2];
	$count++;      }
close can1;



# load reference sequences
my %cds; 
open (hand1, "CDS_pep.fa") or die $!;
while (<hand1>) {
    $_=~ s/\s+$//;
    if (/^>/)       {
      $ge=$_;
      $ge=~ s/^>//;
      $ge =~ s/^.+_//;
      next;}
    $cds{$ge} .= $_; }

close hand1;

my %ref_pep;

# create the reference peptide by trypsin cleavage. Any sequences that is identical to this reference will be discarded. 
foreach my $n (keys %cds)  {
	my @fragment = split /(?<=[RK])(?!P)/, $cds{$n};
	foreach my $i (@fragment)  {
		$ref_pep{$i} = $n if length($i) > 4 or length($i) < 40;  }  
    }

# get the candidate sequences


mkdir "mass-lib";

open (out1, ">mass-lib/mass_mut.$cand.cds.fa");
open (out2, ">mass-lib/mass_mut.$cand.cds.fragment.fa");
my %pep; 
my @aa = qw(A C D E F G H I K L M N P Q R S T V W Y);
   print join(',', @aa)."\n";
foreach my $ncu (keys %can)     {

	if (exists $cds{$ncu} )  {
		my $si = length($cds{$ncu}) - 2;  
		# replacing every aa from 2nd position (not M) to end. * mean stop, so not include)
		for my $i (1 .. $si) {
			my $amino_acid = substr($cds{$ncu}, $i, 1);
			foreach my $j (@aa)  {
				if ($j ne $amino_acid) {
					my $seq = $cds{$ncu};
					substr($seq, $i, 1, $j);
					#print "$seq\n";
					my $name = "$ncu\_$i"."$amino_acid-$j"; 
					print out1 ">$name\n$seq\n";
					# cut rules, after R or K and not followed by P.
					my @fragment = split /(?<=[RK])(?!P)/, $seq;
					
					foreach my $f (@fragment)   {
						next if length($f) < 4 or length($f) > 40;
						next if exists $ref_pep{$f}; 
									
						$pep{$f} .= "$name\t";  }  
					} } }  }  }


foreach my $p (keys %pep)  {
	my @a = split /\t/, $pep{$p};
	next if scalar @a > 1; 
	print out2 ">$pep{$p}\n$p\n"; }


foreach my $i (keys %ref_pep) {
	print out2 ">$ref_pep{$i}\n$i\n"; }

close out1; close out2;

	
						

