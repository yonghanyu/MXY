#!/usr/bin/perl -w
use strict;
die "Fatal Error:perl Call_STRParam.pl population.txt para.xls" if (@ARGV != 2);
my $file = shift;
my $out = shift;
open (IN,"<$file");
open (OUT,">$out");
print OUT "Loci\tcount\tMP\tDP\tPIC\tHE\tPE\n";
while (<IN>) {
	chomp;
	unless (/^#/) {
		my @line = split/\t/;
		my %genoFre;
		my %alleFre;
		my $genoNo = 0;
		for (my $i = 1; $i < @line ; $i++) {
			if ($line[$i] ne '-') {
				$genoNo++;
				my @genos = split(/;/,$line[$i]);
				$genoFre{$genos[0]}++;
				my @al = split(/,/,$genos[0]);
				foreach my $x (@al) {
					$alleFre{$x}++;
				}
			}
			
		}
		foreach my $x (sort keys %genoFre) {
			$genoFre{$x} /= $genoNo;
		}
		foreach my $x (sort keys %alleFre ) {
			$alleFre{$x} /= (2*$genoNo);
		}
		&printAlleFre($line[0],\%alleFre);
		my @geno = values %genoFre;
		my @alle;
		foreach my $x (sort{$a <=> $b} keys %alleFre) {
		 	push (@alle,$alleFre{$x});
		}
		my $dp = &DP_cal(\@geno);
		my $mp = 1 - $dp;
		my $pic = &PIC_cal(\@alle);
		my $he = &HE_cal(\@alle);
		my $pe = &PE_cal(\@alle);
		print OUT "$line[0]\t$genoNo\t$mp\t$dp\t$pic\t$he\t$pe\n";
	}	
}
close IN;

sub DP_cal()
{
	my ($geno) = @_;
	my $pi2 = 0;
	for (my $i = 0; $i < @$geno; $i++) {
		my $pi = $$geno[$i];
		$pi2 += $pi*$pi;
	}
	my $dp = 1 - $pi2;
}

sub PIC_cal()
{
	my ($alle) = @_;
	my $pi2 = 0;
	for (my $i = 0; $i < @$alle; $i++) {
		my $pi = $$alle[$i];
		$pi2 += $pi*$pi;
	}
	my $n = @$alle;
	my $pq2 = 0;
	for (my $i = 0; $i < @$alle-1; $i++) {
		$pq2 += 2*$$alle[$i]*$$alle[$i+1]*$$alle[$i]*$$alle[$i+1];
	}
	my $pic = 1 - $pi2 - $pq2;
}

sub HE_cal()
{
	my ($alle) = @_;
	my $pi2 = 0;
	for (my $i = 0; $i < @$alle; $i++) {
		my $pi = $$alle[$i];
		$pi2 += $pi*$pi;
	}
	my $n = @$alle;
	my $he = $n/($n-1)*(1 - $pi2);
}

sub PE_cal()
{
	my ($alle) = @_;
	my $pe1 = 0;
	my $pe2 = 0;
	for (my $i = 0; $i < @$alle; $i++) {
		my $pi = $$alle[$i];
		$pe1 += $pi*(1-$pi)*(1-$pi);
		for (my $j = 0; $j < @$alle; $j++) {
			if ($j != $i) {
				my $pj = $$alle[$j];
				$pe2 += $pi*$pi*$pj*$pj*(3*$pi+3*$pj-4);
			}
		}
	}
	my $pe = $pe1 + $pe2/2;
}

sub printAlleFre()
{
	my ($loci,$alleFre) = @_;
	print "$loci\tallele\tfreq\n";
	foreach my $x (sort{$a <=> $b} keys %$alleFre) {
		print "$loci\t$x\t$$alleFre{$x}\n";
	}
	print "\n";
}
