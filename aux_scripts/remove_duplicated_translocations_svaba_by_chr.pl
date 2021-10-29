
###############################################################################
# Author: Alex Di Genova 
# Laboratory: IARC/SCG/RCG
# Copyright (c)
# year: 2020
###############################################################################
use Data::Dumper;
use Getopt::Std;
use strict;

sub usage {
   print "$0 usage : -a  -b  -c\n";
   print "Error in use\n";
   exit 1;
}

my %opts = ();
getopts( "a:b:c:", \%opts );
if ( !defined $opts{a}  ) {
   usage;
}

open(FILE, $opts{a}) or die "cannot open file\n";
my $hash=();
while(my $line=<FILE>){
	chomp $line;
	if($line =~m/^#/){
		print $line."\n";
	}else{
		my @d=split("\t",$line);
		$hash->{join("_",$d[0],$d[1])}->{count}++;
		$hash->{join("_",$d[0],$d[1])}->{id}=$d[2];
	}
}
close(FILE);
#second pass to print only one selected pair 
#print Dumper($hash);
my $sids=();

$sids->{$hash->{$_}->{id}}=1  foreach (keys %{$hash});

open(FILE, $opts{a}) or die "cannot open file\n";
my @svs=();
while(my $line=<FILE>){
	chomp $line;
	if($line =~m/^#/){
		#print $line."\n";
	}else{
		my @d=split("\t",$line);
		my ($id)=split(":",$d[2]);
		if(defined $sids->{$id.":1"} or defined $sids->{$id.":2"}){
		#we save in an array to sort by chromosome and position the output	
		my $tmp=();
		$tmp->{chr}=$d[0];
		$tmp->{pos}=$d[1];
		$tmp->{line}=$line;
		$tmp->{chr}=~s/chr//g;
		if($tmp->{chr} eq "X"){
			$tmp->{chr}=23;
		}
		if($tmp->{chr} eq "Y"){
			$tmp->{chr}=24;
		}
		if($tmp->{chr} eq "M"){
			$tmp->{chr}=25;
		}
		 #print $line."\n";
		push(@svs,$tmp);
		}
	}
}

foreach my $sv (sort{$a->{chr} <=>$b->{chr} || $a->{pos} <=> $b->{pos}} @svs){
	print $sv->{line}."\n";
}
