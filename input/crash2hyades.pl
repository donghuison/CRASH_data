#!/usr/bin/perl
use strict;
my $FileIn = $ARGV[0] or die "Usage: crash2hyades.pl CRASHFILE.out";

# Header line for Hyades output files
my %Header = ( "Xe" => 
	       "    371     5.40000000E+01 1.31300000E+02 3.05512951E+00", 
	       "Be" => 
	       "     81     4.00000000E+00 9.01200000E+00 1.84500000E+00",
	       "Pl" =>
	       "     32     3.50000000E+00 6.51000000E+00 1.04400000E+00",
	       "Au" =>
	       "     51     7.90000000E+01 1.96967000E+02 1.93000000E+01");

# Read CRASH data file
open (FILEIN, $FileIn) or die "Could not open FileIn=$FileIn\n";
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
for $i2 (0..$n2-1){
    for $i1 (0..$n1-1){
	# Read line
	my @Numbers = split(' ',<FILEIN>);

	# The first two columns are the 10 based log of coordinates. Read once.
	$Coord1[$i1] = 10.0**$Numbers[0]/1000.0   unless $i2;  # kg/m3 -> g/cm3
	$Coord2[$i2] = 10.0**$Numbers[1]/1.1604e7 unless $i1;  # K     -> keV

	for $iMaterial (0..$nMaterial-1){
	    # The first two columns are followed by pXe EintXe pBe EintBe ...
	    $p[$iMaterial][$i] = 10.*$Numbers[2+2*$iMaterial]; # Pa->dyne/cm2
	    $e[$iMaterial][$i] = 1e4*$Numbers[3+2*$iMaterial]; # J/kg->erg/g
	}
	$i++;
    }
}
close FILEIN;

print "nVar=$nVar nMaterial=$nMaterial n1=$n1 n2=$n2 NameVar=@NameVar\n";

my $iMaterial;

# Write separate Hyades files 1 by 1
for $iMaterial (0..$nMaterial-1){

    $NameVar[2+2*$iMaterial] =~ /(\w\w)$/;
    my $Material = $1;
    my $FileOut = "eos_crash_$Material.dat";

    open (FILEOUT, ">$FileOut") or die "Could not open FileOut=$FileOut\n";

    print FILEOUT "$Material CRASH EOS table for p(rho,T) and rho*e(rho,T)\n";
    printf FILEOUT "%s%8i\n", $Header{$Material}, 2 + $n1 + $n2 + 2*$n1*$n2;
    printf FILEOUT "%15.8e%15.8e",$n1, $n2;
    my $Count = 2;
    for $i1 (0..$n1-1){
	printf FILEOUT "%15.8e", $Coord1[$i1];
	if($Count++ == 4){printf FILEOUT "\n"; $Count=0;}
    }
    for $i2 (0..$n2-1){
	printf FILEOUT "%15.8e", $Coord2[$i2];
	if($Count++ == 4){printf FILEOUT "\n"; $Count=0;}
    }
    for $i (0..$n1*$n2-1){
	printf FILEOUT "%15.8e", $p[$iMaterial][$i];
	if($Count++ == 4){printf FILEOUT "\n"; $Count=0;}
    }
    for $i (0..$n1*$n2-1){
	printf FILEOUT "%15.8e", $e[$iMaterial][$i];
	if($Count++ == 4){printf FILEOUT "\n"; $Count=0;}
    }
    printf FILEOUT "\n" if $Count;

    close FILEOUT;

    print "Finished $FileIn --> $FileOut\n";

}
