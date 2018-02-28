// ======================================================
// CODE TO DETERMINE PARAMETERS NEEDED FOR PROTON PID
// CUTS IN THE CASE OF SIMULATIONS
// ======================================================
#include <fstream>
#include <iostream>
#include <cmath>
#include <sstream>
#include <TStyle.h>
#include"TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include <TLegend.h>

#include <TGraphErrors.h>
#include "Fiducial.h"
#include "Run_dependent.h"
#include "constants.h"
#include "global_variables.h"

using namespace std;

int main(int argc, char ** argv){
	if (argc != 4){
		cerr << "Wrong number of arguments. Instead try\n"
			<< "\tsim_pid_proton /path/to/input/file Ebeam_MeV nuclear_target\n\n";
		return -1;	
	}

	// --------------------------------------------------------------------------------------------------
	// Open up the Root file
	TFile * f = new TFile(argv[1]);
	if (f)
		cerr << "Successfully opened file " << argv[1] << "\n";
	else
	{
		cerr << "Could not open file " << argv[1] << "\n\tExiting...\n";
		return -2;
	}

	istringstream iss (argv[2]);
	int in_num_part;
	iss >> in_num_part;

	// --------------------------------------------------------------------------------------------------
	// Open up the tree, and get the important data
	TTree * t = (TTree*)f->Get("data");
	const int nEvents = t->GetEntries();
	const int maxPart = 50;
	int gPart, CCPart, DCPart, ECPart, SCPart;
	int StatDC[maxPart], StatCC[maxPart], StatEC[maxPart], StatSC[maxPart], id_guess[maxPart];
	float STT, W, Yb;
	float Stat[maxPart], EC_in[maxPart], EC_out[maxPart], EC_tot[maxPart], Nphe[maxPart],
	      SC_Time[maxPart], SC_Path[maxPart],
	      Charge[maxPart], Beta[maxPart], mass[maxPart], mom[maxPart], px[maxPart], py[maxPart],
	      pz[maxPart], theta[maxPart], phi[maxPart], targetZ[maxPart], 
	      EC_X[maxPart],EC_Y[maxPart],EC_Z[maxPart];

	t->SetBranchAddress("gPart"    ,&gPart  ); // Number of particles observed (globally) in the event
	t->SetBranchAddress("CCPart"   ,&CCPart ); // Number of particles observed in the Cherenkovs
	t->SetBranchAddress("DCPart"   ,&DCPart ); // Number of particles observed in the Drift Chambers
	t->SetBranchAddress("ECPart"   ,&ECPart ); // Number of particles observed in the ECal
	t->SetBranchAddress("SCPart"   ,&SCPart ); // Number of particles observed in the Scintillation Counters (ToFs)
	t->SetBranchAddress("Stat"     ,Stat    ); // Global status for each particle candidate
	t->SetBranchAddress("StatDC"   ,StatDC  ); // Drift Chamber status for each particle candidate
	t->SetBranchAddress("StatCC"   ,StatCC  ); // Cherenkov status for each particle
	t->SetBranchAddress("StatEC"   ,StatEC  ); // ECal status for each particle
	t->SetBranchAddress("StatSC"   ,StatSC  ); // Scintillation counter status for each particle  
	t->SetBranchAddress("particle" ,id_guess); // Guess at the particle ID made by the recon software (maybe not reliable)
	t->SetBranchAddress("EC_in"    ,EC_in   ); // Inner layer of ECal for each particle
	t->SetBranchAddress("EC_out"   ,EC_out  ); // Outer layer of ECal for each particle
	t->SetBranchAddress("EC_tot"   ,EC_tot  ); // Total energy deposit in the ECal for each particle
	t->SetBranchAddress("Nphe"     ,Nphe    ); // Number of photo-electrons per hit in the Cherenkov detectors
	t->SetBranchAddress("SC_Time"  ,SC_Time ); // Time in the scintillators per particle
	t->SetBranchAddress("SC_Path"  ,SC_Path ); // Path Length per particle
	t->SetBranchAddress("Charge"   ,Charge  ); // Charge per particle
	t->SetBranchAddress("Beta"     ,Beta    ); // Beta per particle
	t->SetBranchAddress("Mass"     ,mass    ); // Mass per particle
	t->SetBranchAddress("Momentum" ,mom     ); // Momentum magnitude per particle
	t->SetBranchAddress("Momentumx",px      ); // Momentum x component per particle
	t->SetBranchAddress("Momentumy",py      ); // Momentum y component per particle
	t->SetBranchAddress("Momentumz",pz      ); // Momentum z component per particle
	t->SetBranchAddress("Theta"    ,theta   ); // Theta per particle
	t->SetBranchAddress("Phi"      ,phi     ); // Phi per particle
	t->SetBranchAddress("TargetZ"  ,targetZ ); // Target Z per particle
	t->SetBranchAddress("STT"      ,&STT    ); // RF-corrected start time.
	t->SetBranchAddress("EC_X"     ,EC_X    ); // x positions of hit in the calorimeter
	t->SetBranchAddress("EC_Y"     ,EC_Y    ); // y positions of hit in the calorimeter
	t->SetBranchAddress("EC_Z"     ,EC_Z    ); // z positions of hit in the calorimeter

	// ---------------------------------------
	// Diagnostic electron histograms

        TH1D * h1_e_Nphe0    = new TH1D("h1_e_Nphe0"    ,"e- before cuts;# photo-electrons in CC;Counts"     ,200,   0.,200.);
        TH1D * h1_e_Nphe1    = new TH1D("h1_e_Nphe1"    ,"e- passing PID;# photo-electrons in CC;Counts"     ,200,   0.,200.);
        TH1D * h1_e_Nphe2    = new TH1D("h1_e_Nphe2"    ,"e- passing PID+fid;# photo-electrons in CC;Counts" ,200,   0.,200.);

        TH1D * h1_e_EC_in0   = new TH1D("h1_e_EC_in0"   ,"e- before cuts;E_in [GeV];Counts"                  ,300,   0.,  1.);
        TH1D * h1_e_EC_in1   = new TH1D("h1_e_EC_in1"   ,"e- passing PID;E_in [GeV];Counts"                  ,300,   0.,  1.);
        TH1D * h1_e_EC_in2   = new TH1D("h1_e_EC_in2"   ,"e- passing PID+fid;E_in [GeV];Counts"              ,300,   0.,  1.);

        TH1D * h1_e_EC_out0  = new TH1D("h1_e_EC_out0"  ,"e- before cuts;E_out [GeV];Counts"                 ,300,   0.,  1.);
        TH1D * h1_e_EC_out1  = new TH1D("h1_e_EC_out1"  ,"e- passing PID;E_out [GeV];Counts"                 ,300,   0.,  1.);
        TH1D * h1_e_EC_out2  = new TH1D("h1_e_EC_out2"  ,"e- passing PID+fid;E_out [GeV];Counts"             ,300,   0.,  1.);

        TH1D * h1_e_EC_tot0   = new TH1D("h1_e_EC_tot0"  ,"e- before cuts;E_tot [GeV];Counts"                ,300,   0.,  1.);
        TH1D * h1_e_EC_tot1   = new TH1D("h1_e_EC_tot1"  ,"e- passing PID;E_tot [GeV];Counts"                ,300,   0.,  1.);
        TH1D * h1_e_EC_tot2   = new TH1D("h1_e_EC_tot2"  ,"e- passing PID+fid;E_tot [GeV];Counts"            ,300,   0.,  1.);

        TH2D * h2_e_phiTheta0= new TH2D("h2_e_phiTheta0","e- before  cuts;#phi [deg];#theta [deg];Counts"      ,300,-100.,380.,300,10.,50.);
        TH2D * h2_e_phiTheta1= new TH2D("h2_e_phiTheta1","e- passing PID;#phi [deg];#theta [deg];Counts"       ,300,-100.,380.,300,10.,50.);
        TH2D * h2_e_phiTheta2= new TH2D("h2_e_phiTheta2","e- passing PID+fid;#phi [deg];#theta [deg];Counts"   ,300,-100.,380.,300,10.,50.);

        TH2D * h2_e_Ein_Eout0= new TH2D("h2_e_Ein_Eout0","e- before cuts;E_in/p;E_out/p;Counts"              ,300,   0., 0.5,300, 0.,0.5);
        TH2D * h2_e_Ein_Eout1= new TH2D("h2_e_Ein_Eout1","e- passing PID;E_in/p;E_out/p;Counts"              ,300,   0., 0.5,300, 0.,0.5);
        TH2D * h2_e_Ein_Eout2= new TH2D("h2_e_Ein_Eout2","e- passing PID+fid;E_in/p;E_out/p;Counts"          ,300,   0., 0.5,300, 0.,0.5);

        TH2D * h2_e_Ein_Eout_0=new TH2D("h2_e_Ein_Eout_0","e- before cuts;E_in;E_out;Counts"                 ,300,   0., 1.5,300, 0.,1.5);
        TH2D * h2_e_Ein_Eout_1=new TH2D("h2_e_Ein_Eout_1","e- passing PID;E_in;E_out;Counts"                 ,300,   0., 1.5,300, 0.,1.5);
        TH2D * h2_e_Ein_Eout_2=new TH2D("h2_e_Ein_Eout_2","e- passing PID+fid;E_in;E_out;Counts"             ,300,   0., 1.5,300, 0.,1.5);

        TH2D * h2_e_xyEC_hit0= new TH2D("h2_e_xyEC_hit0","e- before cuts;ECx [cm];ECy [cm];Counts"           ,300,-400.,400.,300,-400.,400.);
        TH2D * h2_e_xyEC_hit1= new TH2D("h2_e_xyEC_hit1","e- passing PID;ECx [cm];ECy [cm];Counts"           ,300,-400.,400.,300,-400.,400.);
        TH2D * h2_e_xyEC_hit2= new TH2D("h2_e_xyEC_hit2","e- passing PID+fid;ECx [cm];ECy [cm];Counts"       ,300,-400.,400.,300,-400.,400.);

        TH2D * h2_e_p_Etot0  = new TH2D("h2_e_p_Etot0"  ,"e- before cuts;p [GeV];E_tot/p;Counts"             ,300,   0.,  5.,300, 0.,0.7);
        TH2D * h2_e_p_Etot1  = new TH2D("h2_e_p_Etot1"  ,"e- passing PID;p [GeV];E_tot/p;Counts"             ,300,   0.,  5.,300, 0.,0.7);
        TH2D * h2_e_p_Etot2  = new TH2D("h2_e_p_Etot2"  ,"e- passing PID+fid;p [GeV];E_tot/p;Counts"         ,300,   0.,  5.,300, 0.,0.7);

        TH2D * h2_e_p_E0     = new TH2D("h2_e_p_E0"     ,"e- before cuts;p [GeV];E;Counts"                   ,300,   0.,  5.,300, 0.,2.);
        TH2D * h2_e_p_E1     = new TH2D("h2_e_p_E1"     ,"e- passing PID;p [GeV];E;Counts"                   ,300,   0.,  5.,300, 0.,2.);
        TH2D * h2_e_p_E2     = new TH2D("h2_e_p_E2"     ,"e- before PIF+fid;p [GeV];E;Counts"                ,300,   0.,  5.,300, 0.,2.);

        TH2D * h2_e_thetaMom0= new TH2D("e_thetaMom0"     ,"e- before cuts;#theta [deg];Mom [GeV];Counts"       ,300,  10., 50.,300, 0., 6.);
        TH2D * h2_e_thetaMom1= new TH2D("e_thetaMom1"     ,"e- passing PID;#theta [deg];Mom [GeV];Counts"       ,300,  10., 50.,300, 0., 6.);
        TH2D * h2_e_thetaMom2= new TH2D("e_thetaMom2"     ,"e- passing PID+fid;#theta [deg];Mom [GeV];Counts"   ,300,  10., 50.,300, 0., 6.);

        // ---
	TH1D * h1_e_vz_sec10 = new TH1D("h1_e_vz_sec10" ,"e- passing cuts, before vtx corr, sector 1;electron vz [cm]; Counts" ,300, -10., 10.);
        TH1D * h1_e_vz_sec1  = new TH1D("h1_e_vz_sec1"  ,"e- passing cuts,  after vtx corr, sector 1;electron vz [cm]; Counts" ,300, -10., 10.);
        TH1D * h1_e_vz_sec20 = new TH1D("h1_e_vz_sec20" ,"e- passing cuts, before vtx corr, sector 2;electron vz [cm]; Counts" ,300, -10., 10.);
        TH1D * h1_e_vz_sec2  = new TH1D("h1_e_vz_sec2"  ,"e- passing cuts,  after vtx corr, sector 2;electron vz [cm]; Counts" ,300, -10., 10.);
        TH1D * h1_e_vz_sec30 = new TH1D("h1_e_vz_sec30" ,"e- passing cuts, before vtx corr, sector 3;electron vz [cm]; Counts" ,300, -10., 10.);
        TH1D * h1_e_vz_sec3  = new TH1D("h1_e_vz_sec3"  ,"e- passing cuts,  after vtx corr, sector 3;electron vz [cm]; Counts" ,300, -10., 10.);
        TH1D * h1_e_vz_sec40 = new TH1D("h1_e_vz_sec40" ,"e- passing cuts, before vtx corr, sector 4;electron vz [cm]; Counts" ,300, -10., 10.);
        TH1D * h1_e_vz_sec4  = new TH1D("h1_e_vz_sec4"  ,"e- passing cuts,  after vtx corr, sector 4;electron vz [cm]; Counts" ,300, -10., 10.);
        TH1D * h1_e_vz_sec53 = new TH1D("h1_e_vz_sec53" ,"e- passing cuts, before vtx corr, sector 5;electron vz [cm]; Counts" ,300, -10., 10.);
        TH1D * h1_e_vz_sec5  = new TH1D("h1_e_vz_sec5"  ,"e- passing cuts,  after vtx corr, sector 5;electron vz [cm]; Counts" ,300, -10., 10.);
        TH1D * h1_e_vz_sec60 = new TH1D("h1_e_vz_sec60" ,"e- passing cuts, before vtx corr, sector 6;electron vz [cm]; Counts" ,300, -10., 10.);
        TH1D * h1_e_vz_sec6  = new TH1D("h1_e_vz_sec6"  ,"e- passing cuts,  after vtx corr, sector 6;electron vz [cm]; Counts" ,300, -10., 10.);
        TH1D * h1_e_vz0      = new TH1D("h1_e_vz0"      ,"e- passing cuts, before vtx corr;electron vz [cm]; Counts"           ,300, -10., 10.);
        TH1D * h1_e_vz       = new TH1D("h1_e_vz"       ,"e- passing cuts,  after vtx corr;electron vz [cm]; Counts"           ,300, -10., 10.);
        TH2D * h2_e_phiVz0   = new TH2D("h2_e_phiVz0"   ,"e- passing cuts, before vtx corr;#phi [deg];vz [cm];Counts"   ,300,-100.,380.,300,-10,10);
        TH2D * h2_e_phiVz    = new TH2D("h2_e_phiVz"    ,"e- passing cuts,  after vtx corr;#phi [deg];vz [cm];Counts"   ,300,-100.,380.,300,-10,10);
        TH2D * h2_e_thetaVz0 = new TH2D("h2_e_thetaVz0" ,"e- passing cuts, before vtx corr;#theta [deg];vz [cm];Counts" ,300, -10., 60.,300,-11,11);
        TH2D * h2_e_thetaVz  = new TH2D("h2_e_thetaVz"  ,"e- passing cuts,  after vtx corr;#theta [deg];vz [cm];Counts" ,300, -10., 60.,300,-11,11);
        // ---

        TH1D * h1_e_momCor    = new TH1D("h1_e_momCor"    ,"e- passing fid. cuts;p corrected - p[GeV];Counts"          ,300, -.1,.04 );
        TH1D * h1_e_momCor1   = new TH1D("h1_e_momCor1"   ,"e- passing fid. cuts;p corrected/p;Counts"                 ,300,0.97,1.01);
        TH2D * h2_e_momMomCor = new TH2D("h2_e_momMomCor" ,"e- passing fid. cuts;p [GeV];p corrected - p[GeV];Counts"  ,300,  0.,  6.,300,-.1 ,.04 );
        TH2D * h2_e_momMomCor1= new TH2D("h2_e_momMomCor1","e- passing fid. cuts;p [GeV];p corrected/p;Counts"         ,300,  0.,  6.,300,0.97,1.01);

        TH1D * h1_e_momCor_sec1   = new TH1D("h1_e_momCor_sec1"   ,"e- passing PID+fid sec1;p corrected/p;Counts"        ,300,0.97,1.01);
        TH2D * h2_e_momMomCor_sec1= new TH2D("h2_e_momMomCor_sec1","e- passing PID+fid sec1;p [GeV];p corrected/p;Counts",300, 0., 6.,300,0.97,1.01);
        TH1D * h1_e_momCor_sec2   = new TH1D("h1_e_momCor_sec2"   ,"e- passing PID+fid sec2;p corrected/p;Counts"        ,300,0.97,1.01);
        TH2D * h2_e_momMomCor_sec2= new TH2D("h2_e_momMomCor_sec2","e- passing PID+fid sec2;p [GeV];p corrected/p;Counts",300, 0., 6.,300,0.97,1.01);
        TH1D * h1_e_momCor_sec3   = new TH1D("h1_e_momCor_sec3"   ,"e- passing PID+fid sec3;p corrected/p;Counts"        ,300,0.97,1.01);
        TH2D * h2_e_momMomCor_sec3= new TH2D("h2_e_momMomCor_sec3","e- passing PID+fid sec3;p [GeV];p corrected/p;Counts",300, 0., 6.,300,0.97,1.01);
        TH1D * h1_e_momCor_sec4   = new TH1D("h1_e_momCor_sec4"   ,"e- passing PID+fid sec4;p corrected/p;Counts"        ,300,0.97,1.01);
        TH2D * h2_e_momMomCor_sec4= new TH2D("h2_e_momMomCor_sec4","e- passing PID+fid sec4;p [GeV];p corrected/p;Counts",300, 0., 6.,300,0.97,1.01);
        TH1D * h1_e_momCor_sec5   = new TH1D("h1_e_momCor_sec5"   ,"e- passing PID+fid sec5;p corrected/p;Counts"        ,300,0.97,1.01);
        TH2D * h2_e_momMomCor_sec5= new TH2D("h2_e_momMomCor_sec5","e- passing PID+fid sec5;p [GeV];p corrected/p;Counts",300, 0., 6.,300,0.97,1.01);
        TH1D * h1_e_momCor_sec6   = new TH1D("h1_e_momCor_sec6"   ,"e- passing PID+fid sec6;p corrected/p;Counts"        ,300,0.97,1.01);
        TH2D * h2_e_momMomCor_sec6= new TH2D("h2_e_momMomCor_sec6","e- passing PID+fid sec6;p [GeV];p corrected/p;Counts",300, 0., 6.,300,0.97,1.01);

        TH2D * h2_e_vzVzCor   = new TH2D("h2_e_vzVzCor"   ,"e- passing fid. cuts;vz [cm];vz corrected - vz [cm];Counts",300,-20., 20.,300,-1. , 1. );

        // ---------------------------------------
        // Diagnostic positive particle histograms
        TH1D * h1_p_mass     = new TH1D("h1_pos_mass"   ,"+  passing fid. cuts;mass [GeV];Counts"          ,300,  0.,3.5 );
        TH2D * h2_p_pMass    = new TH2D("h2_pos_pMass"  ,"+  passing fid. cuts;p [GeV];mass [GeV];Counts"  ,300,  0., 4.,300, 0.,3.5);
        TH2D * h2_pos_pBeta  = new TH2D("h2_pos_pBeta"  ,"+  passing fid. cuts;p [GeV];#beta;Counts"       ,300,  0., 4.,300, 0.,1.3);

        // ---------------------------------------
        // Diagnostic proton histograms
	TH2D * h2_p_phiTheta0= new TH2D("h2_p_phiTheta0","p before cuts;#phi [deg];#theta [deg];Counts"    ,300,-100.,380.,300,0.,55.);
        TH2D * h2_p_phiTheta1= new TH2D("h2_p_phiTheta1","p passing fid;#phi [deg];#theta [deg];Counts"    ,300,-100.,380.,300,0.,55.);
        TH2D * h2_p_phiTheta2= new TH2D("h2_p_phiTheta2","p passing fid+PID;#phi [deg];#theta [deg];Counts",300,-100.,380.,300,0.,55.);

        TH2D * h2_p_deltaTmom0=new TH2D("h2_p_deltaTmom0","p before cuts;#Delta t [ns];p [GeV];Counts"     ,300,  -7.,  7.,300, 0., 6.);

        TH2D * h2_p_p_momCor0= new TH2D("h2_p_p_momCor0","p passing fid. cuts;p [GeV];p - p_corr [GeV];Counts"     ,100, 0., 2.5,100,-.06,.01);
        TH2D * h2_p_p_momCor1= new TH2D("h2_p_p_momCor1","p passing fid. cuts;p [GeV];p_corr/p;Counts"             ,100, 0., 2.5,100,0.95,1.2);
        TH2D * h2_p_th_pCor0 = new TH2D("h2_p_th_pCor0" ,"p passing fid. cuts;#theta [deg];p - p_corr [GeV];Counts",100, 0., 150,100,-.06,.01);
        TH2D * h2_p_th_pCor1 = new TH2D("h2_p_th_pCor1" ,"p passing fid. cuts;#theta [deg];p_corr/p;Counts"        ,100, 0., 150,100,0.95,1.2);
        TH2D * h2_p_th_p_cor0= new TH2D("h2_p_th_p_cor0","p passing fid. cuts;#theta [deg];p[GeV];p - p_corr [GeV]",100, 0., 150,100,0.  ,2.5);
        TH2D * h2_p_th_p_cor1= new TH2D("h2_p_th_p_cor1","p passing fid. cuts;#theta [deg];p[GeV];p_corr/p"        ,100, 0., 150,100,0.  ,2.5);

        TH2D * h2_p_vzVzCor  = new TH2D("h2_p_vzVzCor"  ,"p passing fid. cuts;vz [cm];vz corrected - vz [cm];Counts",300, -20., 20.,300,-1., 1 );
        TH2D * h2_p_phiVz0   = new TH2D("h2_p_phiVz0"   ,"p passing cuts, before vtx corr;#phi [deg];vz [cm];Counts",300,-100.,380.,300,-10,10 );
        TH2D * h2_p_phiVz    = new TH2D("h2_p_phiVz"    ,"p passing cuts,  after vtx corr;#phi [deg];vz [cm];Counts",300,-100.,380.,300,-10,10 );
        TH2D * h2_p_thetaVz0 = new TH2D("h2_p_thetaVz0" ,"p passing cuts, before vtx corr;#theta [deg];vz [cm];Counts",300, 0.,100.,300,-11,11 );
        TH2D * h2_p_thetaVz  = new TH2D("h2_p_thetaVz"  ,"p passing cuts,  after vtx corr;#theta [deg];vz [cm];Counts",300, 0.,100.,300,-11,11 );
        TH2D * h2_p_pBeta    = new TH2D("h2_p_pBeta"    ,"p passing fid. cuts;p [GeV];#beta;Counts"                 ,300,   0.,  4.,300, 0.,1.3);

	// Histograms to get proton PID parameters
	const float min_x = 0 ;
        const float max_x = 10.0 ;
	const float min   = 0.0;
	const float max   = 6;
	const int nSlices = 15;
	int sdiv = (int) pdeltat_sig_cutrange;

	TH2F * h2_e_p_Etot    = new TH2F("h2_e_p_Etot"    ,";p [GeV];E_tot/p" ,nSlices, min_x, max_x,300, min, max);
	TH1D ** h1_slices     = new TH1D * [nSlices];

	TH2D * h2_p_deltaTmom1= new TH2D("h2_p_deltaTmom1",Form("p+ after %i#sigma pid cut;#Delta t [ns];p [GeV]",sdiv) ,300,  -7.,  7.,300, 0., 6.);
	
	// ---------------------------------------
	// Setting up output tree and branches
	double e_vz_corrected; //, e_vz, e_mom[3], e_phi_mod;
	double p_vz, p_vz_corrected, p_mom_corrected, p_phi_mod;
	double EC_in_cut, el_EC_cut;
	double e_t0,beta_assuming_proton,p_t0,delta_t;
	double corr_px, corr_py, corr_pz;

	TVector3 e_ec_xyz;
	TVector3 T3_e_mom, T3_e_mom_cor, T3_p_mom, u1;

	// --------------------------------------------------------------------------------------------------
	// Obtaining run number and other important parameters
	int NRun;
	int tab_E1, tab_torus, tab_mini;
	string tab_targ;

	tab_E1    = in_num_part;
	tab_torus = 2250   ;
	tab_mini  = 5996   ;
	tab_targ  = argv[3];

	if      (tab_E1 == 4461) NRun = 17908; // Take a random 4.4 GeV run number
	else if (tab_E1 == 2261) NRun = 18201; // Take a random 2.2 GeV run number

	cout << "Ebeam  = " << tab_E1    << endl << "Torus  = " << tab_torus << endl;
	cout << "Mini   = " << tab_mini  << endl << "Target = " << tab_targ  << endl;

	Fiducial fid_params(tab_E1,tab_torus,tab_mini,tab_targ, false);  // Create an instance of the Fiducial Class
	Run_dependent run_dependent_corrections(NRun);       // Create an instance of the Run_dependent Class

	// --------------------------------------------------------------------------------------------------
        // Open up the output file
        TFile * outfile = new TFile(Form("el_Epratio_mom_%i_sim.root",tab_E1),"RECREATE");

	// Values for some cuts
	if      (tab_E1 == 4461){
		EC_in_cut = 0.055; //GeV (Values for Energy deposited in EC inner layer cut)
		el_EC_cut = 0.330; //GeV (Values for Enough total energy in the EC cut)
	}
	else if (tab_E1 == 2261){
		EC_in_cut = 0.060; //GeV (Values for Energy deposited in EC inner layer cut)
		el_EC_cut = -999.; // No cut in this case (Values for Enough total energy in the EC cut)
	}
	else {
		cout << "Error: Check skim_tree and add parameters for Ebeam = " << tab_E1 << endl;
		exit(-2);
	}

	// --------------------------------------------------------------------------------------------------
	// Loop over events
	for (int event=0; event < nEvents ; event++)
	{
		if (event % 100000 == 0){cerr << "Working on event " << event << " out of " << nEvents << "\n";}

		t->GetEvent(event);	

		if (gPart <= 0) continue; // Ignore events that have no particle candidates	

		// --------------------------------------------------------------------------------------------------
		// Sector index for electrons
		int e_sect = (int)(phi[0]+30)/60;
		if (e_sect>5) e_sect = 5;
		if (e_sect<0) e_sect = 0;

		// --------------------------------------------------------------------------------------------------
		double el_cand_EC = TMath::Max(EC_in[0] + EC_out[0], EC_tot[0]); // Define the electron candidate energy in the EC
		T3_e_mom.SetXYZ(px[0],py[0],pz[0]); // Electron momentum expressed in a TVector3
		e_vz_corrected = targetZ[0]+fid_params.vz_corr(T3_e_mom);
		e_ec_xyz.SetXYZ(EC_X[0],EC_Y[0],EC_Z[0]);

		// ---------------------------------------------------------------------------------------
		// Electron general cuts
		if (!(			(StatEC[0] > 0) && // EC status is good for the electron candidate
					(StatDC[0] > 0) && // DC status is good for the electron candidate
					(StatCC[0] > 0) && // CC status is good for the electron candidate
					(StatSC[0] > 0) && // SC status is good for the electron candidate
					(Charge[0] < 0)    // Electron candidate curvature direction is negative
		     ))
		{continue;}
		// ---------------------------------------------------------------------------------------

		h1_e_Nphe0     -> Fill( Nphe        [0]  );
		h1_e_EC_in0    -> Fill( EC_in       [0]  );
		h1_e_EC_out0   -> Fill( EC_out      [0]  );
		h1_e_EC_tot0   -> Fill( EC_tot      [0]  );
		h2_e_phiTheta0 -> Fill( phi         [0], theta        [0]);
		h2_e_Ein_Eout0 -> Fill( EC_in[0]/mom[0], EC_out[0]/mom[0]);
		h2_e_Ein_Eout_0-> Fill( EC_in[0]       , EC_out[0]       );
		h2_e_p_Etot0   -> Fill( mom         [0], EC_tot[0]/mom[0]);
		h2_e_xyEC_hit0 -> Fill( EC_X        [0], EC_Y         [0]);
		h2_e_p_E0      -> Fill( mom         [0], el_cand_EC      );
		h2_e_thetaMom0 -> Fill( theta       [0], mom          [0]);

		// ---------------------------------------------------------------------------------------
		//Electron particle Identification
		if (!(                  (EC_in [0] > EC_in_cut) &&      // Electron candidate has enough energy deposit in inner layer of EC    
					(el_cand_EC > el_EC_cut) &&     // Enough total energy in the EC
					(fid_params.in_e_EoverP(el_cand_EC/mom[0],mom[0],epratio_sig_cutrange)) // Electron PID (E/p)
		     ))
		{continue;}

		// ---------------------------------------
		// Additional cut for 2GeV data:
		/*
		double el_sccc_dt = SC_Time[0] - CC_Time[0] - (SC_Path[0] - CC_Path[0])/(c_m_s*ns_to_s*100.);

		if(			(tab_E1==2261)&&(	
					CC_Chi2[0]>=0.1 ||
					el_sccc_dt < sc_cc_dt_cut_sect[e_sect] ||
					sqrt(mom[0]*mom[0]+me*me)>tab_E1/1000.
					))
		{continue;}	
		*/
		// ---------------------------------------------------------------------------------------
		// If the event made it here, the electron candidate passed all PID cuts
		h1_e_Nphe1     -> Fill( Nphe        [0]  );
		h1_e_EC_in1    -> Fill( EC_in       [0]  );
		h1_e_EC_out1   -> Fill( EC_out      [0]  );
		h1_e_EC_tot1   -> Fill( EC_tot      [0]  );
		h2_e_phiTheta1 -> Fill( phi         [0], theta        [0]);
		h2_e_Ein_Eout1 -> Fill( EC_in[0]/mom[0], EC_out[0]/mom[0]);
		h2_e_Ein_Eout_1-> Fill( EC_in[0]       , EC_out[0]       );
		h2_e_p_Etot1   -> Fill( mom         [0], EC_tot[0]/mom[0]);
		h2_e_xyEC_hit1 -> Fill( EC_X        [0], EC_Y         [0]);
		h2_e_p_E1      -> Fill( mom         [0], el_cand_EC      );
		h2_e_thetaMom1 -> Fill( theta       [0], mom          [0]);

		// ---------------------------------------------------------------------------------------
		// Electron Fiducial cuts
		if (!fid_params.e_inFidRegion(T3_e_mom)) continue; // Electron theta-phi cut
		if (!fid_params.CutUVW_e(e_ec_xyz)       ) continue; // Cuts on edges of calorimeter (u>60, v<360, w<400);

		// ---------------------------------------------------------------------------------------
		// If electron passes all cuts, then momentum-correct it (only works for theta > 16 deg):
		if (            ((tab_E1==4461)&&(180./M_PI*T3_e_mom.Theta()>16.))||
				( tab_E1==2261)
		   )
			T3_e_mom_cor = fid_params.eMomentumCorrection(T3_e_mom);
		else 	T3_e_mom_cor = T3_e_mom;

		// Fill some diagnostic histograms
		h1_e_Nphe2     -> Fill(Nphe  [0]);
		h1_e_EC_in2    -> Fill(EC_in [0]);
		h1_e_EC_out2   -> Fill(EC_out[0]);
		h1_e_EC_tot2   -> Fill(EC_tot[0]);
		h2_e_xyEC_hit2 -> Fill(EC_X[0],EC_Y[0]);
		h2_e_thetaMom2 -> Fill(theta[0],mom[0]);
		h1_e_momCor    -> Fill(T3_e_mom_cor.Mag()-T3_e_mom.Mag());
		h1_e_momCor1   -> Fill(T3_e_mom_cor.Mag()/T3_e_mom.Mag());
		h2_e_momMomCor -> Fill(T3_e_mom.Mag(),T3_e_mom_cor.Mag()-T3_e_mom.Mag());
		h2_e_momMomCor1-> Fill(T3_e_mom.Mag(),T3_e_mom_cor.Mag()/T3_e_mom.Mag());
		h2_e_vzVzCor   -> Fill(targetZ[0],e_vz_corrected-targetZ[0]);
		h2_e_Ein_Eout2 -> Fill(EC_in[0]/mom[0],EC_out[0]/mom[0]);
		h2_e_Ein_Eout_2-> Fill(EC_in  [0]     , EC_out[0]      );
		h2_e_p_Etot2   -> Fill(mom    [0],EC_tot[0]/mom[0]);
		h2_e_p_E2      -> Fill(mom    [0],el_cand_EC      );
		h2_e_phiTheta2 -> Fill(phi    [0],theta  [0]      );
		h2_e_phiVz0    -> Fill(phi    [0],targetZ[0]      );
		h2_e_phiVz     -> Fill(phi    [0],e_vz_corrected  );
		h2_e_thetaVz0  -> Fill(theta  [0],targetZ[0]      );
		h2_e_thetaVz   -> Fill(theta  [0],e_vz_corrected  );
		h1_e_vz0       -> Fill(targetZ[0]                 );
		h1_e_vz        -> Fill(e_vz_corrected             );
	
		// --------------------------------------------------------------------------------------------------
		// Loop over events looking for other particles
		for (int i=1 ; i<gPart ; i++)
		{	
			T3_p_mom.SetXYZ(px[i],py[i],pz[i]);
			u1 = T3_p_mom.Unit();

			e_t0 = SC_Time[0] - SC_Path[0]/c_cm_ns;

			beta_assuming_proton = mom[i]/sqrt(mom[i]*mom[i] + mP*mP);
			p_t0 = SC_Time[i] - SC_Path[i]/(beta_assuming_proton * c_cm_ns);
			delta_t = p_t0 - e_t0;

			// ------------------------------------------------------------------------------------------
			// Test if positive particle
			if(             (StatSC[i] > 0) && 		// SC status is good for the positive candidate
					(StatDC[i] > 0) &&              // DC status is good for the positive candidate
					(Stat  [i] > 0) &&		// Global status is good for the positive candidate
					(Charge[i] > 0) 		// Charge is positive
			  )
			{

				h2_p_phiTheta0    -> Fill(phi[i]     ,theta[i]);
				h2_p_deltaTmom0   -> Fill(delta_t    ,mom  [i]); 	

				// Passing positive hadron fiducial cuts
				if(fid_params.pFiducialCut(T3_p_mom)){

					// Positive particle vertex (_z) correction
					p_vz_corrected = targetZ[i]+fid_params.vz_corr(T3_p_mom);

					h1_p_mass         -> Fill(mass[i]);
					h2_p_pMass        -> Fill(mom [i]    ,mass [i]);
					h2_pos_pBeta      -> Fill(mom [i]    ,Beta [i]);
					h2_p_phiTheta1    -> Fill(phi [i]    ,theta[i]);
					h2_p_deltaTmom1   -> Fill(delta_t    ,mom  [i]);
					
					// --------------------------------------------------------------------
					// Look specifically for protons
					if((id_guess[i] == 2212 ) &&       // Guess at the particle ID is good for the proton candidate
							(fid_params.in_p_deltaT(delta_t, mom[i], pdeltat_sig_cutrange)) // Proton PID (delta T vs p)
					  ){
						if(run_dependent_corrections.ProtonMomCorrection_He3_4Cell(T3_p_mom,p_vz_corrected) != -1)
							p_mom_corrected=run_dependent_corrections.ProtonMomCorrection_He3_4Cell(T3_p_mom,p_vz_corrected);
						else	p_mom_corrected=mom[i];	

						corr_px = p_mom_corrected*u1.X();
						corr_py = p_mom_corrected*u1.Y();
						corr_pz = p_mom_corrected*u1.Z();

						T3_p_mom.SetXYZ(corr_px,corr_py,corr_pz);
						
						/*
						h2_p_deltaTmom2-> Fill(delta_t   ,mom                   [i]);
						h2_p_phiTheta2 -> Fill(phi    [i],theta                 [i]);
						h2_p_vzVzCor   -> Fill(targetZ[i],p_vz_corrected-targetZ[i]);
						h2_p_p_momCor0 -> Fill(mom    [i],mom [i]-p_mom_corrected  );
						h2_p_p_momCor1 -> Fill(mom    [i],p_mom_corrected/mom   [i]);
						h2_p_th_pCor0  -> Fill(theta  [i],mom [i]-p_mom_corrected  );
						h2_p_th_pCor1  -> Fill(theta  [i],p_mom_corrected/mom   [i]);
						h2_p_th_p_cor0 -> Fill(theta  [i],mom[i],mom[i]-p_mom_corrected);
						h2_p_th_p_cor1 -> Fill(theta  [i],mom[i],p_mom_corrected/mom[i]);
						h2_p_phiVz0    -> Fill(phi    [i],targetZ               [i]);
						h2_p_phiVz     -> Fill(phi    [i],p_vz_corrected           );
						h2_p_thetaVz0  -> Fill(theta  [i],targetZ               [i]);
						h2_p_thetaVz   -> Fill(theta  [i],p_vz_corrected	   );
						h2_p_pBeta     -> Fill(mom    [i],Beta                  [i]);
						*/
					}
				}
			}
		}
	}
	cerr << "Finished with the event loop...\n";

	gStyle->SetOptStat(0);

	// -------------------------------------------------------------------------------------------------- 
	// Doing fits
	TF1 * f_gaus  = new TF1("f_gaus" ,"gaus",min,max);
	TF1 * f_mean  = new TF1("f_mean" ,"pol5",min,max);
	TF1 * f_sig   = new TF1("f_sig"  ,"pol5",min,max);

	double bin_ctr [nSlices];
	double fit_mean[nSlices];
	double fit_sd  [nSlices];

	for(int i = 0 ; i < nSlices ; i++){
		bin_ctr  [i] = h2_e_p_Etot -> GetXaxis() -> GetBinCenter(i+1);
                h1_slices[i] = h2_e_p_Etot -> ProjectionY(Form("%i",i),i+1,i+2);
		h1_slices[i] -> SetTitle(Form("Bin center: p = %.2f GeV",bin_ctr[i]));

                h1_slices[i] -> Fit(f_gaus,"","",min,max);

                fit_mean[i] = f_gaus -> GetParameter(1);
                fit_sd  [i] = f_gaus -> GetParameter(2);
        }
	
	TGraphErrors * g_fit_mean = new TGraphErrors(nSlices,bin_ctr,fit_mean,0,0);
	g_fit_mean -> SetTitle("mean");
	g_fit_mean -> GetXaxis() -> SetTitle("p [GeV]");
	g_fit_mean -> GetYaxis() -> SetTitle("#mu");

	TGraphErrors * g_fit_sd   = new TGraphErrors(nSlices,bin_ctr,fit_sd  ,0,0);
	g_fit_sd -> SetTitle("standard deviation");
        g_fit_sd -> GetXaxis() -> SetTitle("p [GeV]");
        g_fit_sd -> GetYaxis() -> SetTitle("#sigma");

	// Fitting mean and standard deviation with polynomials
	g_fit_mean -> Fit(f_mean ,"","",min_x,max_x);
	g_fit_sd   -> Fit(f_sig  ,"","",min_x,max_x);
	
	// --------------------------------------------------------------------------------------------------
        // Loop over events
	double delta_t_lo_limit, delta_t_up_limit, eff_mom;

	cout << "Looping over events a second time to check the cuts" << endl;
        for (int event=0; event < nEvents ; event++)
        {
                if (event % 100000 == 0){cerr << "Working on event " << event << " out of " << nEvents << "\n";}
                t->GetEvent(event);
                if (gPart <= 0) continue; // Ignore events that have no particle candidates
                double el_cand_EC = TMath::Max(EC_in[0] + EC_out[0], EC_tot[0]); // Define the electron candidate energy in the EC      

                // ---------------------------------------------------------------------------------------
                // Electron general cuts
                if (!(                  (StatEC[0] > 0) && // EC status is good for the electron candidate
                                        (StatDC[0] > 0) && // DC status is good for the electron candidate
                                        (StatCC[0] > 0) && // CC status is good for the electron candidate
                                        (StatSC[0] > 0) && // SC status is good for the electron candidate
                                        (Charge[0] < 0)    // Electron candidate curvature direction is negative
                     ))
                {continue;}
                // ---------------------------------------------------------------------------------------
		eff_mom = mom[0];
		if(eff_mom > max_x) eff_mom = max_x;

		delta_t_lo_limit = f_mean->Eval(eff_mom) - epratio_sig_cutrange * f_sig->Eval(eff_mom);
                delta_t_up_limit = f_mean->Eval(eff_mom) + epratio_sig_cutrange * f_sig->Eval(eff_mom);

		if((EC_tot[0]/mom[0]>delta_t_lo_limit)&&(EC_tot[0]/mom[0]<delta_t_up_limit)){
                	h2_e_p_Etot1   -> Fill( mom         [0], EC_tot[0]/mom[0]);
		}
        }
        cerr << "Finished with the event loop...\n";

	// -------------------------------------------------------------------------------------------------- 
	TCanvas *c1  = new TCanvas("c1" );	h1_e_Nphe0      -> Draw();
	TCanvas *c2  = new TCanvas("c2" );	h1_e_EC_in0     -> Draw();
	TCanvas *c3  = new TCanvas("c3" );	h1_e_EC_out0    -> Draw();
	TCanvas *c4  = new TCanvas("c4" );	h1_e_EC_tot0    -> Draw();
	TCanvas *c5  = new TCanvas("c5" );	h2_e_thetaMom0  -> Draw("COLZ"); 
	TCanvas *c6  = new TCanvas("c6" );	h2_e_Ein_Eout0  -> Draw("COLZ");
	TCanvas *c7  = new TCanvas("c7" );	h2_e_Ein_Eout_0 -> Draw("COLZ");
	TCanvas *c8  = new TCanvas("c8" );	h2_e_xyEC_hit0  -> Draw("COLZ");
	TCanvas *c9  = new TCanvas("c9" );	h2_e_p_Etot0    -> Draw("COLZ");
	TCanvas *c10 = new TCanvas("c10");	h2_e_p_E0       -> Draw("COLZ");
	TCanvas *c11 = new TCanvas("c11");	h2_e_phiTheta0  -> Draw("COLZ");
	TCanvas *c12 = new TCanvas("c12");	h2_e_p_Etot     -> Draw("COLZ");

	TCanvas *c13 = new TCanvas("c13");
	c13 -> Divide(nSlices/3,3);
	for(int i = 0 ; i < nSlices ; i++){	c13 -> cd(i+1);	h1_slices[i] -> Draw();}

	TCanvas *c14 = new TCanvas("c14");
        c14 -> Divide(2,1);
	c14 -> cd(1);	g_fit_mean -> Draw("AL");
	c14 -> cd(2);	g_fit_sd   -> Draw("AL");

	TCanvas *c15 = new TCanvas("c15");
        c15 -> Divide(2,1);
        c15 -> cd(1);	h2_e_p_Etot0 -> Draw("COLZ");
	c15 -> cd(2);	h2_e_p_Etot1 -> Draw("COLZ");

	// --------------------------------------------------------------------------------------------------
	// Print histograms on a pdf file
	c1  -> Print(Form("./e2a_maps_%d.pdf(",tab_E1) ,"pdf");
	c2  -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c3  -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c4  -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c5  -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");   
	c6  -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c7  -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c8  -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c9  -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c10 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c11 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c12 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c13 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c14 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c15 -> Print(Form("./e2a_maps_%d.pdf)",tab_E1) ,"pdf");
	
	// --------------------------------------------------------------------------------------------------

	// --------------------------------------------------------------------------------------------------
	// Write the output file
	outfile->cd();	
	f_mean ->Write();
	f_sig  ->Write();

	// Clean up
	f->Close();
	outfile->Close();

	return 0;
}
