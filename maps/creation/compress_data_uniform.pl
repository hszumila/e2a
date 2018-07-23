#!/usr/bin/perl -w

use strict;
use warnings;

my @runs;

for(my $i=101;$i <=150; $i++){
    push @runs,$i;
}

#print "@runs \n";

foreach my $runnumber (@runs){
    
    print STDOUT "$runnumber\n";
    
    open (OUTFILE, '>', "/volatile/clas/clase2/ahrnjic/jobs/simulation/f_01_${runnumber}") or die "Can't open output file\n";
    
    print OUTFILE "PROJECT: e2a\n";
    print OUTFILE "JOBNAME: f_01_$runnumber\n";
    print OUTFILE "OS: centos7\n";
    print OUTFILE "DISK_SPACE: 15 GB\n";
    print OUTFILE "MEMORY: 2000 MB\n";
    print OUTFILE "TRACK: simulation\n\n";
    
    print OUTFILE "COMMAND: /u/home/ahrnjic/barak_stuff/sim.csh $runnumber\n\n";
    
    print OUTFILE "INPUT_FILES: /u/home/ahrnjic/barak_stuff/convert_uniform\n\n";
    
    print OUTFILE "OTHER_FILES:\n";
    print OUTFILE "/u/home/ahrnjic/barak_stuff/ffread.in\n";
    print OUTFILE "/u/home/ahrnjic/barak_stuff/e2_pass2.tcl\n\n";
  
    print OUTFILE "OUTPUT_DATA: gsimtest.txt\n";
    print OUTFILE "OUTPUT_TEMPLATE: /volatile/clas/clase2/ahrnjic/outfiles/simulation/info/gsiminfo_$runnumber.txt\n\n";

    print OUTFILE "OUTPUT_DATA: gpp_test.txt\n";
    print OUTFILE "OUTPUT_TEMPLATE: /volatile/clas/clase2/ahrnjic/outfiles/simulation/info/gppinfo_$runnumber.txt\n\n";

    print OUTFILE "OUTPUT_DATA: recsis_test.txt\n";
    print OUTFILE "OUTPUT_TEMPLATE: /volatile/clas/clase2/ahrnjic/outfiles/simulation/info/recsisinfo_$runnumber.txt\n\n";

    print OUTFILE "OUTPUT_DATA: look.txt\n";
    print OUTFILE "OUTPUT_TEMPLATE: /volatile/clas/clase2/ahrnjic/outfiles/simulation/info/bosdump_$runnumber.txt\n\n";

    print OUTFILE "OUTPUT_DATA: outfile.root\n";
    print OUTFILE "OUTPUT_TEMPLATE: /volatile/clas/clase2/ahrnjic/outfiles/simulation/f_01_$runnumber.root\n";
    
    close(OUTFILE);
}

