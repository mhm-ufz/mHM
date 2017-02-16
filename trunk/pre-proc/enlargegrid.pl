#!/usr/bin/perl
# SHOW CURRENT COLS, ROWS:  enlargegrid.pl FILE.asc
# ENLARGE TO COLS, ROWS:    enlargegrid.pl 50 40
# Adds columns and rows to all .asc files in the current directory
$|=1;
my $NODATA = '-9999';
if ($ARGV[0]!~/^\d+$/) {
  my $datastart = 0; my %H = ();

  open(my $F, shift) or die $!;
  while (<$F>) {
    if ( /^NODATA/ ) { $datastart++; next; }
    if ($datastart>0) { my @col = split /\s+/, $_; $H{$#col+1}++; }}
  close $F;

  foreach (keys %H) { print sprintf "%4s rows with %4s cols.\n", $H{$_}, $_; }
}
else {
  my ($x, $y) = @ARGV; my $datastart = 0; my ($xdiff, $ydiff) = -1;
  {
    local $^I = ".bak"; # backup files
    local @ARGV = glob("*.asc");
    print "$_\n" foreach @ARGV;
    while (<>) {
      if (/^ncols/) { $datastart=0;$xdiff=-1;$ydiff=-1; print STDOUT "$ARGV\n"; }
      if ($datastart) { s/[\r\n]+//; print $_.(" $NODATA" x $xdiff)."\n";  }
      else {
        $xdiff = $x-$2 if s/^(ncols\s+)(\d+)/$1$x/;
        $ydiff = $y-$2 if s/^(nrows\s+)(\d+)/$1$y/;
        print;
        if ( /^NODATA_value\s+(\S+)/ ) { $NODATA = $1 if $1; $datastart++;
          print ''.($NODATA.(" $NODATA" x ($x-1))."\n")x$ydiff; }
} } } }
                                                                        