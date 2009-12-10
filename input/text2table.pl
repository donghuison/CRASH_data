#!/usr/bin/perl -n
# Convert from text to table format and from CGS to SI units
# Usage: text2table.pl Beryllium_table.txt > Beryllium_table.out

print
"Opacities from $ARGV: [kg/m3] [K] [m2/kg]
0 0.0 2 1 1
31 50
0.0
logrho logT Opacity none
" if($. == 1);

$temp = $1 if /Temperature\s*=\s*(\S+)/;
if(s/g\/cc -/$temp/){
    @a = split(' ',$_);
    printf "%12.5e %12.5e %12.5e\n",log($a[0]*1000),log($a[1]),$a[2]*0.1;
}
