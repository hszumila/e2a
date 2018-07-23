#!/bin/csh -f

#If running locally for testing
#unlink InputFile
#rm -f histfile gsimtest.txt gsimtest.bos
#rm -f gpp.*
#rm -f anamonhist logfile recsis_test.txt
#rm -f outfile1 outfile.root look.txt

#Step 1:Generator
if ($#argv != 1) then
./convert_uniform -l -b -q
else
./convert_uniform -l -b -q $1
endif

/u/home/baraks/bin/Linux64RHEL7/txt2part -m -omctk_uniform.evt < mctk_uniform.txt

#Step 2:GSIM
#source /u/home/baraks/packages/environment.csh
setenv CLAS_PARMS /group/clas/parms
setenv RECSIS_RUNTIME /group/clas/clsrc/recsis/runtime

setenv CLAS_CALDB_RUNINDEX calib.RunIndex

echo "Running..."
echo $CLAS_CALDB_RUNINDEX

/u/home/baraks/bin/Linux64RHEL7/gsim_bat -ffread ffread.in -mcin mctk_uniform.evt -bosout gsimtest.bos >& gsimtest.txt

#Step 3:GPP
setenv CLAS_CALDB_RUNINDEX calib_user.RunIndexe2a  # checked calib_user.RunIndexe2 also, but gave seg. fault
echo $CLAS_CALDB_RUNINDEX

/u/home/baraks/bin/Linux64RHEL7/gpp  -P0x1f -Y -ogpp_test.bos -a1.5 -b1.5 -c1.5 -f0.1 -R18016 gsimtest.bos >& gpp_test.txt
#/u/home/baraks/bin/Linux64RHEL7/gpp -P0x0 -ogpp_test.bos -R18016 gsimtest.bos >& gpp_test.txt #efficiency is set to 1

#Step 4:Recsis
setenv CLAS_CALDB_RUNINDEX calib.RunIndex
echo $CLAS_CALDB_RUNINDEX

ln -s gpp_test.bos InputFile

/volatile/halla/gmp12/baraks/e2_compile/bin/Linux64RHEL7/user_ana -t ./e2_pass2.tcl >& recsis_test.txt

/u/home/baraks/bin/Linux64RHEL7/bosdump outfile1 >& look.txt

#Step 5:ClasTool
#source /u/home/baraks/.cshrc
$CLASTOOL/bin/Linux64RHEL7/WriteRootDst -GSIM -o outfile.root outfile1


