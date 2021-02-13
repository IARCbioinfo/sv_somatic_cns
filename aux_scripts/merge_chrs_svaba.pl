###############################################################################
# Author: Alex Di Genova
# Laboratory: IARC/SCG/RCG
# Copyright (c)
# year: 2021
###############################################################################
use Data::Dumper;
use Getopt::Std;
use strict;

sub usage {
   print "$0 usage : -s sampleID -p <prefix>\n";
   print "Error in use\n";
   exit 1;
}

my %opts = ();
getopts( "s:p:", \%opts );
if ( !defined $opts{s} or !defined $opts{p}){
   usage;
}

my @chrom=("chr1","chr2","chr3", "chr4",
           "chr5","chr6", "chr7", "chr8",
           "chr9", "chr10", "chr11", "chr12",
           "chr13","chr14","chr15", "chr16",
           "chr17", "chr18", "chr19", "chr20",
           "chr21", "chr22","chrX","chrY");


#s_run_chr1.alignments.txt.gz  #s_run_chr1.svaba.germline.indel.vcf	      s_run_chr1.svaba.unfiltered.germline.sv.vcf
#s_run_chr1.bps.txt.gz	      #s_run_chr1.svaba.germline.sv.vcf		      s_run_chr1.svaba.unfiltered.somatic.indel.vcf
#s_run_chr1.contigs.bam	      #s_run_chr1.svaba.somatic.indel.vcf	      s_run_chr1.svaba.unfiltered.somatic.sv.vcf
#s_run_chr1.discordant.txt.gz  #s_run_chr1.svaba.somatic.sv.vcf
#s_run_chr1.log		      #s_run_chr1.svaba.unfiltered.germline.indel.vcf

my @svaba_e=("svaba.germline.indel.vcf","svaba.germline.sv.vcf",
              "svaba.somatic.indel.vcf","svaba.somatic.sv.vcf",
              "svaba.unfiltered.germline.indel.vcf","svaba.unfiltered.germline.sv.vcf",
              "svaba.unfiltered.somatic.indel.vcf","svaba.unfiltered.somatic.sv.vcf",
              "alignments.txt.gz","bps.txt.gz",
              "log");
              #"discordant.txt.gz");
#we merge all the chromosomes
foreach my $e(@svaba_e){
    my ($merge)=check_file($opts{p},$e,@chrom);
    #we merge the vcf files
    if($e =~m/\.vcf/){
        merge_vcf($opts{s},$e,@$merge);
    }elsif($e=~m/bps.txt.gz/){
        merge_bps($opts{s},$e,@$merge);
    }elsif($e=~m/alignments.txt.gz/){
        merge_aln($opts{s},$e,@$merge);
    }elsif($e=~/log/){
       merge_log($opts{s},$e,@$merge);
    }else{
      #nothing to do with files
    }
}

#merge the log files
sub merge_log{
  my ($sample, $e,@files)=@_;
  my $out=join(".",$sample, $e);
  my $cmd=join(" ","cat ",@files," > ",$out);
  #print $cmd."\n";
  system($cmd);
  if ($? == -1) {
    print "failed to execute: $!\n";
    exit 1;
  }
}


#we merge the aln file
sub merge_aln{
  my ($sample,$e,@files)=@_;
  my $out=join(".",$sample, $e);
  my $cmd=join(" ","gzip -dc ",@files," | gzip > ",$out);
  #print $cmd."\n";
  system($cmd);
  if ($? == -1) {
    print "failed to execute: $!\n";
    exit 1;
  }
}


#function that merge bps files and create the genome
#file at sample level
sub merge_bps{
    my ($sample, $e,@files)=@_;
    #print Dumper($sample,$e,@files);
    my $out=join(".",$sample, $e);
    $out=~s/.gz//;
    open(OUT,">".$out) or die "Cannot create $out file\n";
    for(my $i=0;$i<scalar(@files);$i++){
        open(IN,"gzip -dc $files[$i] | ") or die "cannot open $files[$i]\n";
        while (my $line=<IN>) {
          # we print the header of the file
          if($i==0 and $line=~m/mapq1/){
            print OUT $line;
          }
          next if($line=~m/mapq1/);
          #we print the regular lines
          print OUT $line;
        }
        #we close the CHR vcf file
        close(IN);
    }
    #we close the out file
    close(OUT);
    system("gzip $out");
    if ($? == -1) {
      print "failed to execute: $!\n";
      exit 1;
    }
}

#Function that merge VCF files
#and create the genome file at sample level
sub merge_vcf{
    my ($sample, $e,@files)=@_;
    #print Dumper($sample,$e,@files);
    my $out=join(".",$sample, $e);
    open(OUT,">".$out) or die "cannot create $out file\n";
    for(my $i=0;$i<scalar(@files);$i++){
        open(IN,$files[$i]) or die "cannot open $files[$i]\n";
        while (my $line=<IN>) {
          # we print the header of the file
          if($i==0 and $line=~m/^#/){
            print OUT $line;
          }
          next if($line=~m/^#/);
          #we print the regular lines
          print OUT $line;
        }
        #we close the CHR vcf file
        close(IN);
    }
    #we close the out file
    close(OUT);
}


#Good
sub check_file{
    my ($p,$e,@chrom)=@_;
    my $files=();
    foreach my $c (@chrom){
      my $f=join(".",$p."_".$c,$e);
      if(-e $f){
        push(@{$files},$f);
      }else{
        print "$f do not exist\n";
        exit 1;
        #push(@{$files},$f);
      }
    }
    return $files;
}
