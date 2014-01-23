#!/usr/bin/perl
die "I NEED TWO ARGUMENTS:
pathto/FinalParams.out pathto/mhm_parameters.nml
" if $#ARGV != 1;
print "Reads Finalparams.out and mhm_parameters.nml,
merges new values into mhm_parameters.nml.opt
shows analysis and ranges.\n";
print "".("-"x113)."\n";
printf "%-38s | %9s << %17s -> %22s << %9s |%s\n", "NAME", "LOWER", "OLD VALUE", "NEW VALUE", "UPPER", " WARNING";
print "".("-"x113)."\n";

my $in = shift; my $out = shift;
my @line = ();
open my $IN, $in;
while(<$IN>) {
  s/^\s+//;
  my @a = ();
  push @a, split(/\s+/, $_);
  push @line, \@a; }
close $IN;

my $i=0;
open my $OUT, ">$out.opt";
open my $READ, $out;
while (<$READ>) {
  s/^\s+//;
  my @b = ();
  push @b, split(/\s+/, $_);
  if (/^\!/ or $#b != 6 or $i >= $#{$line[0]} ) { print $OUT $_; next; }
  if (/^(\S+)\s+=\s+(\S+),\s+(\S+),\s+(\S+),(.+)$/) {
    $i++;
    my $edge = " 1% close to edge!" if abs(${$line[1]}[$i]-$1) < abs(0.01*($3-$1)) or abs(${$line[1]}[$i]-$3) < abs(0.01*($3-$1));
    #print "! WARNING: p$i names differ: ".${$line[0]}[$i]." $1\n" if ${$line[0]}[$i] ne $1;
    printf "%-38s | %9s << %17s -> %22s << %9s |%s\n", $1, $2, $4, ${$line[1]}[$i], $3, $edge;
    printf $OUT "%-38s = %9s, %9s, %22s,%s\n", $1, $2, $3, ${$line[1]}[$i], $5; }}
close $OUT;
close $READ;
