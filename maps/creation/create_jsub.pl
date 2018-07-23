#!/usr/bin/perl -w

use strict;
use warnings;

my @filenames;

@filenames = `ls -1 /volatile/clas/clase2/ahrnjic/jobs/simulation/*`;

open (OUTFILE, '>', "JSUB") or die "Can't open test output file\n";

foreach my $filename (@filenames){
    print OUTFILE "jsub $filename";
}

system("chmod","744","JSUB");
