#!/usr/bin/perl
use strict;
my $FileIn = $ARGV[0] or die "Usage: hyades2crash.pl HYADESFILE.dat";
my $FileOut = $FileIn; $FileOut =~ s/.dat$/.out/;

print "$FileIn --> $FileOut\n";

open (FILEIN, $FileIn) or die "Could not open FileIn=$FileIn\n";
# Ignore header
<FILEIN>;
<FILEIN>;
my @Numbers = split(' ', join('', <FILEIN>));
close FILEIN;

print "nNumbers = ",scalar @Numbers,"\n";

my $n1 = int shift @Numbers;
my $n2 = int shift @Numbers;

print "Number of coordinates: $n1, $n2\n";

my @Coord1 = splice(@Numbers,0,$n1);
my @Coord2 = splice(@Numbers,0,$n2);
my @Var1   = splice(@Numbers,0,$n1*$n2);
my @Var2   = splice(@Numbers,0,$n1*$n2);

print "Remaining numbers = ",scalar @Numbers,"\n";

$Coord1[0] = $Coord1[1]**2/$Coord1[2] if $Coord1[0] == 0; 
$Coord2[0] = $Coord2[1]**2/$Coord2[2] if $Coord2[0] == 0; 

print "Coord1[0..2]=@Coord1[0..2]\n";
print "Coord2[0..2]=@Coord2[0..2]\n";

open (FILEOUT, ">$FileOut") or die "Could not open FileOut=$FileOut\n";
print FILEOUT "Hyades $FileIn, units: [kg/m3] [K] [Pa] [J/kg]\n";
print FILEOUT "0 0.0 -2 1 2\n";
print FILEOUT "$n1 $n2\n";
print FILEOUT "0.0\n";
print FILEOUT "logRho logTe p Eint version\n";
my $i;
my $j;
for $j (0..$n2-1){
    for $i (0..$n1-1){
	my $Coord1 = log($Coord1[$i]*1000)/log(10.);     # g/cm3 -> kg/m3
	my $Coord2 = log($Coord2[$j]*1.1604e7)/log(10.); # keV -> K
	my $Var1 = shift @Var1; $Var1 =  0.1*$Var1;      # dyne/cm2 -> Pa;
	my $Var2 = shift @Var2; $Var2 = 1e-4*$Var2;      # erg/g    -> J/kg;
	printf FILEOUT "%E %E %E %E\n", $Coord1, $Coord2, $Var1, $Var2;
    }
}
close FILEOUT;
