#!/usr/bin/env perl
#
# Generate Makefile object dependency file %.d from files ending on .f90, .f, and .for
#
# Usage: make.d.pl SRCFILE [OBJPATH] [SRCPATH] [INCPATH]
#
# LICENSE
#    This file is part of the UFZ makefile project.
#
#    The UFZ makefile project is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    The UFZ makefile project is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public License
#    along with the UFZ makefile project. If not, see <http://www.gnu.org/licenses/>.
#
#    Copyright 2010-2012 - Uwe Schulzweida, Matthias Cuntz
#    Modified from Uwe Schulzweida's createMakefiles.pl for Echam4
#
$PROG=`basename "$0"`;
#
# Create Makefile
#
# Dependency listings
if (@ARGV < 1) {
  print "Error make.d.pl: source file must be given.\n";
  print "    make.d.pl SRCFILE [OBJPATH] [SRCPATH]\n";
  exit;
}
$srcfile = "$ARGV[0]";
if (@ARGV >= 2) {
  $objpath = "$ARGV[1]/";
}
else {
  $objpath = "";
}
if (@ARGV >= 3) {
  $srcpath = "$ARGV[2]/";
}
else {
  $srcpath=`dirname "$srcfile"`;
  #($srcpath = $srcfile) =~ s/.*\///; # does not work yet -> look for dirname
}
if ("$srcpath" eq "") {
  $srcpath = "./";
}
if (@ARGV >= 4) {
  $incpath = "$ARGV[3]/";
}
else {
  $incpath = "./";
}
($srcfile1 = $srcfile) =~ s/.*\///;     # rm leading path
($dfile1 = $srcfile1)  =~ s/\.f90$/.d/; # .f90 -> .d
($dfile1 = $dfile1)    =~ s/\.f$/.d/;   # .f   -> .d
($dfile1 = $dfile1)    =~ s/\.for$/.d/; # .for -> .d
$dfile = "$objpath$dfile1";             # .d in objpath
($ofile1 = $srcfile1) =~ s/\.f90$/.o/; # .f90 -> .o
($ofile1 = $ofile1)   =~ s/\.f$/.o/;   # .f   -> .o
($ofile1 = $ofile1)   =~ s/\.for$/.o/; # .for -> .o
$ofile = "$objpath$ofile1";            # .o in objpath
#print " SF:$srcfile SP:$srcpath OP:$objpath IP:$incpath S1:$srcfile1 D1:$dfile1 D:$dfile O1:$ofile1 O:$ofile\n";
open(MAKEFILE, "> $dfile");
&MakeDependsf90($ARGV[2]);
#
exit;
#
# &PrintWords(current output column, extra tab?, word list); --- print words nicely
#
sub PrintWords {
   local($columns) = 78 - shift(@_);
   local($extratab) = shift(@_);
   local($wordlength);
   #
   print MAKEFILE @_[0];
   $columns -= length(shift(@_));
   foreach $word (@_) {
      $wordlength = length($word);
      if ($wordlength + 1 < $columns) {
         print MAKEFILE " $word";
         $columns -= $wordlength + 1;
         }
      else {
         #
         # Continue onto a new line
         #
         if ($extratab) {
            print MAKEFILE " \\\n\t\t$word";
            $columns = 62 - $wordlength;
            }
         else {
            print MAKEFILE " \\\n\t$word";
            $columns = 70 - $wordlength;
            }
         }
      }
   }

#
# &toLower(string); --- convert string into lower case
#
sub toLower {
   local($string) = @_[0];
   $string =~ tr/A-Z/a-z/;
   $string;
   }

#
# &uniq(sorted word list); --- remove adjacent duplicate words
#
sub uniq {
   local(@words);
   foreach $word (@_) {
      if ($word ne $words[$#words]) {
         push(@words, $word);
         }
      }
   @words;
   }

#
# &MakeDependsf90(f90 compiler); --- FORTRAN 90 dependency maker
#
sub MakeDependsf90 {
   local(@dependencies);
   local(%filename);
   local(@incs);
   local(@modules);
   #
   # Associate each module with the name of the file that contains it
   #
   chdir($srcpath);
   foreach $file (<*.f90>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         /^\s*module\s+([^\s!]+)/i &&
            ($filename{&toLower($1)} = $file) =~ s/\.f90$/.o/;
         }
      }
   foreach $file (<*.f>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         /^\s*module\s+([^\s!]+)/i &&
            ($filename{&toLower($1)} = $file) =~ s/\.f$/.o/;
         }
      }
   foreach $file (<*.for>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         /^\s*module\s+([^\s!]+)/i &&
            ($filename{&toLower($1)} = $file) =~ s/\.for$/.o/;
         }
      }
   #
   # Print the dependencies of the input file, with include's and modules
   #
   open(FILE, $srcfile1);
   while (<FILE>) {
     /^\s*#*include\s+["\']([^"\']+)["\']/i && push(@incs, $1);
     /^\s*use\s+([^\s,!]+)/i && push(@modules, &toLower($1));
   }

   if (defined @incs) {
     #print "  Incs\n";
     $jn = $#incs;
     while ($jn >= 0) {
       if (-e "${incpath}$incs[$jn]") {
	 #print "  $incfile: ${incpath}$incs[$jn]\n";
	 $incs[$jn] = "${incpath}$incs[$jn]";
       }
       else {
	 print "  $file: include file $incs[$jn] not in ${incpath}\n";
	 pop(@incs);
       }
       $jn--;
     }
   }
   if (defined @incs && $#incs < 0) {undef @incs;}

   if (defined @modules) {
     #print "  Modules\n";
     $jn = $#modules;
     while ($jn >= 0) {
       if ("$ofile1" eq "$filename{$modules[$jn]}") {
	 #print "  $objfile: $modules[$jn]\n";
	 pop(@modules);
       }
       $jn--;
     }
   }
   if (defined @modules && $#modules < 0) {undef @modules;}

   if (defined @incs || defined @modules) {
     print MAKEFILE "$ofile $dfile : $srcfile ";
     undef @dependencies;
     foreach $module (@modules) {
       if ("$filename{$module}" ne "") {
	 push(@dependencies, "$objpath$filename{$module}");
       }
     }
     @dependencies = &uniq(sort(@dependencies));
     &PrintWords(length($objfile) + 2, 0,
		 @dependencies, &uniq(sort(@incs)));
     print MAKEFILE "\n";
     undef @incs;
     undef @modules;
   }
   else {
     print MAKEFILE "$ofile $dfile : $srcfile\n";
   }
 }
