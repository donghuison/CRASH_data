#!/usr/bin/perl
use strict;
my $FileIn = $ARGV[0] or die "Usage: crash2hyades.pl CRASHFILE.out";

open (FILEIN, $FileIn) or die "Could not open FileIn=$FileIn\n";

# Ignore header
<FILEIN>;
$_ = <FILEIN>; 
/(\d+)$/;
my $nVar      = $1;
my $nMaterial = $nVar/2;
$_ = <FILEIN>;
/^\s*(\d+)\s+(\d+)/;
my $n1 = $1;
my $n2 = $2;
<FILEIN>;
$_ = <FILEIN>;
my @NameVar = split;
my $i2;
my $i1;
my $i;
my $iMaterial;
my @Coord1;
my @Coord2;
my @p;
my @e;
for $i2 (0..$n2-2){
    for $i1 (0..$n1-1){
	my @Numbers = split(' ',<FILEIN>);

	# The first two columns are the 10 based log of coordinates
	$Coord1[$i1] = 10.0**$Numbers[0]/1000.0;         # kg/m3 -> g/cm3
	$Coord2[$i2] = 10.0**$Numbers[1]/1.1604e7;       # K     -> keV

	for $iMaterial (0..$nMaterial-1){
	    # The first two columns are followed by pXe EintXe pBe EintBe ...
	    $p[$iMaterial][$i] = 10.*$Numbers[2+2*$iMaterial];  # Pa->dyne/cm2
	    $e[$iMaterial][$i] = 1e4*$Numbers[2+2*$iMaterial+1];# J/kg->erg/g
	}
	$i++;
    }
}
close FILEIN;

print "nVar=$nVar nMaterial=$nMaterial n1=$n1 n2=$n2 NameVar=@NameVar\n";

my $iMaterial;

for $iMaterial (0..$nMaterial-1){

    $NameVar[2+2*$iMaterial] =~ /(\w\w)$/;
    my $Material = $1;
    my $FileOut = "eos_hyades_$Material.dat";

    open (FILEOUT, ">$FileOut") or die "Could not open FileOut=$FileOut\n";

    print FILEOUT "$Material CRASH EOS table for p(rho,T) and rho*e(rho,T)\n";
    print FILEOUT "some numbers ", 2 + $n1 + $n2 + 2*$n1*$n2, "\n";
    printf FILEOUT "%e %e ",$n1, $n2;
    my $Count = 2;
    for $i1 (0..$n1-1){
	printf FILEOUT "%e ", $Coord1[$i1];
	if($Count++ == 4){printf FILEOUT "\n"; $Count=0;}
    }
    for $i2 (0..$n2-1){
	printf FILEOUT "%e ", $Coord2[$i2];
	if($Count++ == 4){printf FILEOUT "\n"; $Count=0;}
    }
    for $i (0..$n1*$n2-1){
	printf FILEOUT "%e ", $p[$iMaterial][$i];
	if($Count++ == 4){printf FILEOUT "\n"; $Count=0;}
    }
    for $i (0..$n1*$n2-1){
	printf FILEOUT "%e ", $e[$iMaterial][$i];
	if($Count++ == 4){printf FILEOUT "\n"; $Count=0;}
    }
    printf FILEOUT "\n" if $Count;

    close FILEOUT;

    print "Finished $FileIn --> $FileOut\n";

}
