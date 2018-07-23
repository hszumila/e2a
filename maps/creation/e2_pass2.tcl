### 10/18/09: Turn GEM package ON, Solenoid current 534A;
### 12/08/09: Add TPCC bosbank for RTPC calibration purposes
### 01/21/10: Add/change tcl flags to take into acount the change on building the TBID bank.
### 09/21/12: Add new tcl flag to fill the new RTPC ntuple & Add RTPC bank on the list of outbanknames
### source /group/clas/builds/PRODUCTION/packages/tcl/recsis_proc.tcl;
### 12/15/13: Ifarm change from 32 to 64 bit, try different tcl link

source /group/clas/builds/test3/src/clas6-trunk/reconstruction/recsis/recsis_proc.tcl;
#
# define packages
turnoff ALL;
global_section off;
turnon seb trk cc tof egn user;
set lec1_do -1;
#
#set trigger_particle_s 211;
#
inputfile           InputFile;
setc chist_filename histfile;
setc log_file_name  logfile;
#setc helfile        hel_table.txt;
#
# Allows to keep track of the run ext number.
#
#set runfile 1;
#
#
# Remove either TRPB (layer 1 info or ECPC EC calibration info)
#
setc outbanknames(1) "LCPBHEADRTPCSCRCRUNCTBIDCL01TRPBTBTRCC01CCPBSHHBSHPBDCPBECPBEVNTFBPMHEVTICHBICPBPARTSCPBTBERECHBGCPBBMPRVERTMVRTTGBIEC1REC1 MCTKMCVX";
outputfile outfile1 PROC1 2047;
# For CC calib test, add CCRC bank in addition to CC01 which is already included in (off)
#setc outbanknames(1) "HEADRTPCSCRCHELCRUNCTBIDCL01TRPBTBTRCCRCCC01CCPBSHHBSHPBDCPBECPBEVNTFBPMHEVTICHBICPBPARTSCPBTBERECHBGCPBBMPRVERTMVRTTGBIHLS ";
#outputfile outfile1 PROC1 2047;
#
#setc prlink_file_name  "prlink_dcvs_centered_tgt10cm.bos";
setc prlink_file_name  "prlink_90_75.bos";
setc bfield_file_name  "bgrid_T67to33.fpk";
#
setc poltarget_field_file_name "bgrid_s.fpk";
#
#set torus_current       750;
set torus_current       2250;
#set torus_current       -2250;
#set torus_current      1500;
#set poltarget_current   450000;
set mini_torus_current  6000;
#set TargetPos(3)       -64.0;
#set TargetMagPos(3)    -64.0;
#
set ntswitch 1;
#
set touch_id 0;

# Franz's tcl variables
set trk_maxiter       8;
set trk_minhits(1)    2;
set trk_lrambfit_chi2 50.;
set trk_tbtfit_chi2   70.;
set trk_prfit_chi2    70.;
set trk_statistics    3 ;
#
# New tracking:
set trk_fitRegion   7;
#
#for simulations only
set dc_xvst_choice  0;
#
set lseb_nt_do      -1;
set lall_nt_do      -1;
set lseb_hist       -1;
set lseb_h_do       -1;
set lmon_hist       -1;
#
# SC & CC time cut
#
set CUT_T_SCCC 300.;
#
# New tcl flag to get RF time from RFT bank
#
set lseb_get_rf     -1;
# New set of a sampling fraction after the EC Gains cal. with a bit
# higher gain factors for the energy deposited on the Outer Layers
# set sampl_frac 0.3; It didn't work! Stepan suggested to change it 
# on evnt_set_def.F routine
#
# CC flags for CC calib test
#
#set lcc_do        -1;
#set lcc_nt_do     -1;
#set lccr_nt_do    -1;
#
#
#
# PID and TBID flags
#
set lpid_make_trks   0;
set ltbid_nost_do    0;
set ltbid_do        -1;    
#
# EC flags
#
set inner_surf 3.;
set outer_surf 17.;
#
#
# Add EC reconstruction banks for Bayram
#
#set legn_do	-1;
#set lechb_nt_do -1;
#
# Turn MYSQL ON
#
#set lmysql          -1;
#set nmysql          -1;
#
# Turn MYSQL OFF
#
#set lmysql          0;
#set nmysql          0;
#
# helicity information
#
#set lhelcor        -1;
#
# tell FPACK not to stop if it thinks you are running out of time
#
fpack "timestop -9999999999"
#
#
# do not send events to event display
set lscat $false;
set ldisplay_all $false;
go 2000000;
#go 10000;
#
#
#
setc rec_prompt "CLASCHEF_recsis> ";
exit_pend;
