#!/usr/bin/perl -w
use strict;
use warnings;
my $file = shift;
open (IN,"<$file");
while (<IN>) {
	my @genos = split;
	my $loci = shift(@genos);
	my @alles;
	my %alleFre;
	my $samNo = 0;
	foreach my $geno (@genos) {
		my @alle;
		if ($geno ne "-") {
			$samNo++;
			my @genoinfo = split(/;/,$geno);
			@alle = split(/,/,$genoinfo[0]);
			if ($loci =~ /ab/) {
				foreach my $x (@alle) {
					push (@alles,$x);
				}
			}else{
				push (@alles,$alle[0]);
			}
		}
		
	}
	foreach (@alles){
		$alleFre{$_}++;
	}
	print "$loci\n";
	foreach my $alle (sort{ $a <=> $b } keys %alleFre){
			my $cnt = $alleFre{$alle};
			$alleFre{$alle} = $alleFre{$alle}/@alles;
			print "$alle\t$alleFre{$alle}\t$cnt\n";
	}
	my @a = values %alleFre;
	my $gd = &GD_cal(\@a,$samNo);
	print "$loci\tGD\t$gd\n";
}
close IN;

sub PE_cal()
{
	my ($alle) = @_;
	my $pe1 = 0;
	my $pe2 = 0;
	for (my $i = 0; $i < @$alle; $i++) {
		my $pi = $$alle[$i];
		$pe1 += $pi*(1-$pi)*(1-$pi);
		for (my $j = $i+1; $j < @$alle; $j++) {
			my $pj = $$alle[$j];
			$pe2 += $pi*$pi*$pj*$pj*(3*$pi+3*$pj-4);
		}
	}
	my $pe = $pe1 + $pe2/2;
}

sub GD_cal()
{
	my ($alle,$n) = @_;
	my $pi2 = 0;
	for (my $i = 0; $i < @$alle; $i++) {
		$pi2 += $$alle[$i]*$$alle[$i];
	}
	my $gd = $n*(1-$pi2)/($n-1);
}


