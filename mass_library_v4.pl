#!usr/bin/perl

use warnings;
use strict;
use Bio::Util::codonUsage qw(translate cai);

# set the number of candidate based on the expression level (emPAI)
my $cand = 2000;
my $shift = "-2";

# load the candidate
my $ge; my %can; my $count = 0;
open (can1, "emPAI_gene.csv") or die $!;
while (<can1>)     {
	next if /^name/;
	last if $count > $cand;
	chomp;
	my @a = split /,/;
	$can{$a[0]} = $a[2];
	$count++;      }
close can1;



# load reference sequences
my %gene; my %cds; 
open (hand1, "v12T0_ORF_100.fa") or die $!;
while (<hand1>) {
    $_=~ s/\s+$//;
    if (/^>/)       {
      $ge=$_;
      $ge=~ s/^>//;
      next;}
    my $seq = $_;
    my $orf = substr($seq,100,(length($seq)-200));
    
    my ($dna, $cds) = trans_pep($orf);
    #print "loading $ge sequence\n$orf\n$cds\n";
    #print "$seq\n$cds\n";
    $gene{$ge}=$seq;
    $cds{$ge}=$cds; }

close hand1;

# remove the items from %gene and %cds that are not in %can
foreach my $id (keys %gene)  {
	delete $gene{$id} if !exists $can{$id}; }
foreach my $id (keys %cds)  {
	delete $cds{$id} if !exists $can{$id};  }


open (out1, ">shift-ratio-all/mass-lib/mass_$shift.overlapped.$cand.genes.txt");
open (out2, ">shift-ratio-all/mass-lib/mass_$shift.peptide.$cand.genes.fa");
print out1 "peptide\tfrom\tmapped_gene\n"; 
my %pep; 
   
foreach my $ncu (sort keys %can)     {
	   print "processing $ncu\n";
	   next if !exists $gene{$ncu};
       my $len = length($gene{$ncu}) if exists $gene{$ncu};
       #step 1: retrieve DNA sequence from upstream
       for (my $i=103; $i < $len - 100; $i+=3)  {
          my $upstream = substr($gene{$ncu}, 100, $i-100) if exists $gene{$ncu};
          my ($updna, $uppep) = trans_pep($upstream);
          next if $uppep =~ /[RK]$/;
          $uppep =~ s/^.*[RK]//g;  # cleave the upstream sequence 
          $updna = substr($updna,-(length($uppep)*3));

       #step 2: retrieve DNA sequence from downstream
          my $shift_start;
          if ($shift eq "+1")   {
             $shift_start = $i + 1; }
          elsif ($shift eq "-1") {
             $shift_start = $i - 1; }
          else {
			 $shift_start = $i - 2; }

          my $downstream = substr($gene{$ncu}, $shift_start, int((length($gene{$ncu})-$shift_start)/3)*3) if exists $gene{$ncu}; 
          my ($downdna, $downpep) = trans_pep($downstream);
          
          if ($downpep =~ /^(.*?[RK])/)    {   # perl defaultly is greedy search, i.e. give the longest match. to get shortest one, append a ? to * so that force it to get the shortest match
              $downpep = $1; }  # cleave the downstream sequence (trypsin)
          next if $downpep eq "*";   # not include the case that shift directly cause stop.
          $downdna = substr($downdna,0, (length($downpep)*3));
          

          my $pep = $uppep.$downpep;
          
          next if length($pep) < 7;
          
          my $id = "$ncu\_$shift\_pos$i";   # the ID contain the gene, shift and position information. 
          
          #step 3: authenticate if protein sequence is equal to in frame sequence of any other proteins in genome
          my $set = "on";
          my $pattern = $pep;
          if ($pattern =~ /(.+)\*$/)  {
			  $pattern = $1.'\*';
			  #print "$pattern\n"; 
			   }
			  
          Search: foreach my $j (keys %cds)   {
             if ($cds{$j} =~ /($pattern)/)  {
                 print out1 "$pattern\t$id\t$j\n";
                 $set = "off";
                 last Search; } }
          next if $set eq "off";

          
          $pep{$pep} .= "$id\t";
              }
                                  }


# report results
# first report the full length in frame
foreach my $id (keys %cds)   {
	print out2 ">$id\n$cds{$id}\n";  }

# the shift peptide	       
foreach my $seq (keys %pep)   {
          my @a = split /\t/, $pep{$seq};
          my $r = scalar(@a)-1;
          my $id =$a[0]."r$r";
          print out2 ">$id\n$seq\n";   }


         

close out1; close out2;     
              

       
sub trans_pep     {

    my $seq = shift;

          my $pep; my $dna;
          TRANSLATE: for (my $j=0; $j < length($seq); $j+=3)     {
              my $cod = substr($seq, $j, 3);
              $dna .= $cod;
              my $aa = translate($cod);
              $pep .= $aa ;   
              last TRANSLATE if $aa eq '*';                   }
              #print "$dna\n$pep\n";
    return ($dna, $pep);       }

              

       


