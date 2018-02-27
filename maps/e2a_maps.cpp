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

#include "Fiducial.h"
#include "Run_dependent.h"
#include "constants.h"
#include "global_variables.h"

using namespace std;

// Difference between positive corrected z vertex and electron corrected z vertex
const double pos_z_cut_min = -2;
const double pos_z_cut_max =  2;

// Parameters for 2GeV data cuts
const double sc_cc_dt_cut_sect[6]={-2,-5,-8,-8,-2,2};

int main(int argc, char ** argv){
	if (argc != 5){
		cerr << "Wrong number of arguments. Instead try\n"
			<< "\te2a_maps /path/to/input/file /path/to/output/file Ebeam_MeV nuclear_target\n\n";
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

	istringstream iss (argv[3]);
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
	      SC_Time[maxPart], SC_Path[maxPart], CC_Time[maxPart], CC_Path[maxPart],
	      Charge[maxPart], Beta[maxPart], mass[maxPart], mom[maxPart], px[maxPart], py[maxPart],
	      pz[maxPart], theta[maxPart], phi[maxPart], targetZ[maxPart], theta_pq[maxPart],
	      EC_X[maxPart],EC_Y[maxPart],EC_Z[maxPart], EC_U[maxPart],EC_V[maxPart],EC_W[maxPart], 
	      CC_Chi2[maxPart];

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
	t->SetBranchAddress("CC_Time"  ,CC_Time ); // Time in the cherenkov per particle
	t->SetBranchAddress("CC_Path"  ,CC_Path ); // Path Length per particle
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
	t->SetBranchAddress("EC_U"     ,EC_U    ); // u positions of hit in the calorimeter
	t->SetBranchAddress("EC_V"     ,EC_V    ); // v positions of hit in the calorimeter
	t->SetBranchAddress("EC_W"     ,EC_W    ); // w positions of hit in the calorimeter
	t->SetBranchAddress("CC_Chi2"  ,CC_Chi2 ); // angle between CC hit and nearest SC hit (in rad)

	// Generated variables
	int num_g, id_g[maxPart];
	float z_g[maxPart], theta_g[maxPart], phi_g[maxPart], p_g[maxPart];

	t->SetBranchAddress("Number_g"  ,&num_g );
	t->SetBranchAddress("particle_g",id_g   );
	t->SetBranchAddress("TargetZ_g" ,z_g    );
	t->SetBranchAddress("Theta_g"   ,theta_g);
	t->SetBranchAddress("Phi_g"     ,phi_g  );
	t->SetBranchAddress("Momentum_g",p_g    );

	//t->SetBranchAddress("W"        ,&W      ); // Hadronic mass
	//t->SetBranchAddress("Yb"       ,&Yb     ); // Y-scaling variable

	/* 
	   the kinematic variables: Q2, W, Xb, Nu, Yb are calculated incorrectly in the particle_data.root file
	   (since they assume the eg2c beam energy as hard-coded in the Tidentificator library).
	   So, need to calculate them by hand.
	 */

	// --------------------------------------------------------------------------------------------------
	// Open up the output file
	TFile * outfile = new TFile(argv[2],"RECREATE");

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
	TH1D * h1_e_vz0      = new TH1D("h1_e_vz0"      ,"e- passing cuts, before vtx corr;electron vz [cm]; Counts"           ,300, -10., 10.);
	TH1D * h1_e_vz       = new TH1D("h1_e_vz"       ,"e- passing cuts,  after vtx corr;electron vz [cm]; Counts"           ,300, -10., 10.);
	TH2D * h2_e_phiVz0   = new TH2D("h2_e_phiVz0"   ,"e- passing cuts, before vtx corr;#phi [deg];vz [cm];Counts"   ,300,-100.,380.,300,-10,10);
	TH2D * h2_e_thetaVz0 = new TH2D("h2_e_thetaVz0" ,"e- passing cuts, before vtx corr;#theta [deg];vz [cm];Counts" ,300, -10., 60.,300,-11,11);

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

	TH2D * h2_p_deltaTmom0=new TH2D("h2_p_deltaTmom0","p before cuts;#Delta t [ns];p [GeV];Counts"     ,300,  -7.,  7.,300, 0., 3.);
	TH2D * h2_p_deltaTmom1=new TH2D("h2_p_deltaTmom1","p passing fid;#Delta t [ns];p [GeV];Counts"     ,300,  -7.,  7.,300, 0., 3.);
	TH2D * h2_p_deltaTmom2=new TH2D("h2_p_deltaTmom2","p passing fid+PID;#Delta t [ns];p [GeV];Counts" ,300,  -7.,  7.,300, 0., 3.);

	TH2D * h2_p_phiVz0   = new TH2D("h2_p_phiVz0"   ,"p passing cuts, before vtx corr;#phi [deg];vz [cm];Counts",300,-100.,380.,300,-10,10 );
	TH2D * h2_p_thetaVz0 = new TH2D("h2_p_thetaVz0" ,"p passing cuts, before vtx corr;#theta [deg];vz [cm];Counts",300, 0.,100.,300,-11,11 );

	TH2D * h2_p_pBeta    = new TH2D("h2_p_pBeta"    ,"p passing fid. cuts;p [GeV];#beta;Counts"                 ,300,   0.,  4.,300, 0.,1.3);

	// ---------------------------------------
	// Diagnostic pi+ histograms
	TH2D * h2_pip_pBeta     = new TH2D("h2_pip_pBeta"     ,"#pi+ passing fid. cuts;p [GeV];#beta;Counts"      ,300,   0.,  4.,300, 0.,1.3);
	TH2D * h2_pip_deltaTmom0= new TH2D("h2_pip_deltaTmom0","#pi+ before cuts;#Delta t [ns];p [GeV];Counts"    ,300,  -7.,  7.,300, 0.,5.0);
	TH2D * h2_pip_deltaTmom1= new TH2D("h2_pip_deltaTmom1","#pi+ passing fid;#Delta t [ns];p [GeV];Counts"    ,300,  -7.,  7.,300, 0.,5.0);
	TH2D * h2_pip_deltaTmom2= new TH2D("h2_pip_deltaTmom2","#pi+ passing fid+PID;#Delta t [ns];p [GeV];Counts",300,  -7.,  7.,300, 0.,5.0);

	// ---------------------------------------
	// Diagnostic pi- histograms
	TH2D * h2_pim_pBeta     = new TH2D("h2_pim_pBeta"     ,"#pi- passing fid. cuts;p [GeV];#beta;Counts"      ,300,   0.,  4.,300, 0.,1.3);
	TH2D * h2_pim_deltaTmom0= new TH2D("h2_pim_deltaTmom0","#pi- before cuts;#Delta t [ns];p [GeV];Counts"    ,300,  -7.,  7.,300, 0.,5.0);
	TH2D * h2_pim_deltaTmom1= new TH2D("h2_pim_deltaTmom1","#pi- passing fid;#Delta t [ns];p [GeV];Counts"    ,300,  -7.,  7.,300, 0.,5.0);
	TH2D * h2_pim_deltaTmom2= new TH2D("h2_pim_deltaTmom2","#pi- passing fid+PID;#Delta t [ns];p [GeV];Counts",300,  -7.,  7.,300, 0.,5.0);

	// ---------------------------------------
	// Histograms for maps
	Float_t phi_low[] = {-30, 30,  90, 150, 210, 270};
	Float_t phi_hi [] = { 30, 90, 150, 210, 270, 330};

	TH1F * h1a = new TH1F("h1a", "All Generated Electrons: Z Vertex;Z_{e} [cm]"         ,100,-10,10);
	TH1F * h1b = new TH1F("h1b", "Accepted Electrons: Generated Z Vertex;Z_{e} [cm]"    ,100,-10,10);
	TH1F * h1c = new TH1F("h1c", "Accepted Electrons: Reconstructed Z Vertex;Z_{e} [cm]",100,-10,10);

	TH2F * h2a = new TH2F("h2a","All Generated Electrons: #theta vs. #phi;#phi [deg];#theta [deg]"         ,370,-35,335,180,0,90);
	TH2F * h2b = new TH2F("h2b","Accepted Electrons: Generated #theta vs. #phi;#phi [deg];#theta [deg]"    ,370,-35,335,180,0,90);
	TH2F * h2c = new TH2F("h2c","Accepted Electrons: Reconstructed #theta vs. #phi;#phi [deg];#theta [deg]",370,-35,335,180,0,90);

	TH2F * h3a = new TH2F("h3a","All Generated Electrons: Momentum vs. #theta;#theta [deg];Momentum [GeV]"         ,180,0,90,200,-0.05,5.05);
	TH2F * h3b = new TH2F("h3b","Accepted Electrons: Generated Momentum vs. #theta;#theta [deg];Momentum [GeV]"    ,180,0,90,200,-0.05,5.05);
	TH2F * h3c = new TH2F("h3c","Accepted Electrons: Reconstructed Momentum vs. #theta;#theta [deg];Momentum [GeV]",180,0,90,200,-0.05,5.05);

	TH1F * hp1a = new TH1F("hp1a", "Generated Protons: Z Vertex;Z_{p} [cm]"             ,100,-10,10);
	TH1F * hp1b = new TH1F("hp1b", "Accepted Protons: Generated Z Vertex;Z_{p} [cm]"    ,100,-10,10);
	TH1F * hp1c = new TH1F("hp1c", "Accepted Protons: Reconstructed Z Vertex;Z_{p} [cm]",100,-10,10);

	TH2F * hp2a = new TH2F("hp2a","Generated Protons: #theta vs. #phi;#phi [deg];#theta [deg]"             ,370,-35,335,180,0,180);
	TH2F * hp2b = new TH2F("hp2b","Accepted Protons: Generated #theta vs. #phi;#phi [deg];#theta [deg]"    ,370,-35,335,180,0,180);
	TH2F * hp2c = new TH2F("hp2c","Accepted Protons: Reconstructed #theta vs. #phi;#phi [deg];#theta [deg]",370,-35,335,180,0,180);

	TH2F * hp3a = new TH2F("hp3a","Generated Protons: Momentum vs. #theta;#theta [deg];Momentum [GeV/c]"             ,180,0,180,200,-0.05,3.75);
	TH2F * hp3b = new TH2F("hp3b","Accepted Protons: Generated Momentum vs. #theta;#theta [deg];Momentum [GeV/c]"    ,180,0,180,200,-0.05,3.75);
	TH2F * hp3c = new TH2F("hp3c","Accepted Protons: Reconstructed Momentum vs. #theta;#theta [deg];Momentum [GeV/c]",180,0,180,200,-0.05,3.75);

	TH2F * h4a [6], * h4b [6];
	TH2F * hp4a[6], * hp4b[6];

	for(int i=0; i<6; i++) {
		h4a [i]= new TH2F(Form("h4a[%d]" ,i),Form("e #theta vs #phi: 2.0 GeV < P < 3.5 GeV, Sect %d;#phi [deg];#theta [deg]",i+1),30,phi_low[i],phi_hi[i],15,10,55);
		h4b [i]= new TH2F(Form("h4b[%d]" ,i),Form("e #theta vs #phi: 2.0 GeV < P < 3.5 GeV, Sect %d;#phi [deg];#theta [deg]",i+1),30,phi_low[i],phi_hi[i],15,10,55);
		hp4a[i]= new TH2F(Form("hp4a[%d]",i),Form("p #theta vs #phi: 1.0 GeV < P < 2.5 GeV, Sect %d;#phi [deg];#theta [deg]",i+1),30,phi_low[i],phi_hi[i],35,10,150);
		hp4b[i]= new TH2F(Form("hp4b[%d]",i),Form("p #theta vs #phi: 1.0 GeV < P < 2.5 GeV, Sect %d;#phi [deg];#theta [deg]",i+1),30,phi_low[i],phi_hi[i],35,10,150);
	}

	// ---------------------------------------
	// Setting up output tree and branches
	double e_vz, e_vz_corrected, e_mom[3], e_phi_mod;
	double p_vz, p_vz_corrected, p_mom_corrected, p_phi_mod;
	double EC_in_cut, el_EC_cut;
	double e_t0,beta_assuming_proton,p_t0,delta_t,beta_assuming_pip,pip_t0,pip_delta_t;
	double corr_px, corr_py, corr_pz, n_px, n_py, n_pz, n_p, EC_Path_corr, Beta_corr;

	TVector3 e_ec_xyz, n_ec_xyz;
	TVector3 T3_e_mom, T3_e_mom_cor, T3_p_mom, u1;

	int nParticles;
	int nProtons, nNeutrons, nPiplus, nPiminus, nPi0; 

	// --------------------------------------------------------------------------------------------------
	// Obtaining run number and other important parameters
	//t->GetEvent(0);
	int NRun, tab_E1, tab_torus, tab_mini;
	string tab_targ;

	tab_E1    = in_num_part;
	tab_torus = 2250   ;
	tab_mini  = 5996   ;
	tab_targ  = argv[4];

	if      (tab_E1 == 4461) NRun = 17908; // Take a random 4.4 GeV run number
	else if (tab_E1 == 2261) NRun = 18201; // Take a random 2.2 GeV run number

	cout << "Run    = " << NRun      << endl;
	cout << "Ebeam  = " << tab_E1    << endl;
	cout << "Torus  = " << tab_torus << endl;
	cout << "Mini   = " << tab_mini  << endl;
	cout << "Target = " << tab_targ  << endl;

	Fiducial     fid_params    (tab_E1,tab_torus,tab_mini,tab_targ, false);	// Create an instance of the Fiducial     Class	
	Run_dependent run_dependent_corrections(NRun);				// Create an instance of the Run_dependent Class

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

		if(num_g>0 && id_g[0]==11){ //This should be always true 

			//Generated Kinematic Variables
			h1a->Fill(z_g    [0]           );
			h2a->Fill(phi_g  [0],theta_g[0]);
			h3a->Fill(theta_g[0],p_g[0]    );

			for(int j = 0 ; j < 6 ; j++){
				if(p_g[0]>2.0 && p_g[0]<3.5 && phi_g[0]>phi_low[j] && phi_g[0]<phi_hi[j])
					h4a[j]->Fill(phi_g[0],theta_g[0]);
			}
		}

		if (gPart <= 0) continue; // Ignore events that have no particle candidates

		nParticles = 0;
		nProtons   = 0;     
		nNeutrons  = 0;
		nPiplus    = 0;
		nPiminus   = 0;
		nPi0       = 0;

		// --------------------------------------------------------------------------------------------------
		// Sector index for electrons
		int e_sect = (int)(phi[0]+30)/60;
		if (e_sect>5) e_sect = 5;
		if (e_sect<0) e_sect = 0;

		// --------------------------------------------------------------------------------------------------
		double el_cand_EC = TMath::Max(EC_in[0] + EC_out[0], EC_tot[0]); // Define the electron candidate energy in the EC
		T3_e_mom.SetXYZ(px[0],py[0],pz[0]); // Electron momentum expressed in a TVector3
		e_vz_corrected = targetZ[0];
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
		double el_sccc_dt = SC_Time[0] - CC_Time[0] - (SC_Path[0] - CC_Path[0])/(c_m_s*ns_to_s*100.);

		if(			(tab_E1==2261)&&(	
					CC_Chi2[0]>=0.1 ||
					el_sccc_dt < sc_cc_dt_cut_sect[e_sect] ||
					sqrt(mom[0]*mom[0]+me*me)>tab_E1/1000.
					))
		{continue;}	

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
		if (!fid_params.CutUVW_e(e_ec_xyz)     ) continue; // Cuts on edges of calorimeter (u>60, v<360, w<400);

		// ---------------------------------------------------------------------------------------
		// If electron passes all cuts, then momentum-correct it (only works for theta > 16 deg):
		T3_e_mom_cor = T3_e_mom;	

		nParticles++;	

		// Fill some diagnostic histograms	
		h1_e_Nphe2     -> Fill(Nphe  [0]);
		h1_e_EC_in2    -> Fill(EC_in [0]);
		h1_e_EC_out2   -> Fill(EC_out[0]);
		h1_e_EC_tot2   -> Fill(EC_tot[0]);
		h2_e_xyEC_hit2 -> Fill(EC_X[0],EC_Y[0]);
		h2_e_thetaMom2 -> Fill(theta[0],mom[0]);	
		h2_e_Ein_Eout2 -> Fill(EC_in[0]/mom[0],EC_out[0]/mom[0]);
		h2_e_Ein_Eout_2-> Fill(EC_in  [0]     , EC_out[0]      );
		h2_e_p_Etot2   -> Fill(mom    [0],EC_tot[0]/mom[0]);
		h2_e_p_E2      -> Fill(mom    [0],el_cand_EC      );
		h2_e_phiTheta2 -> Fill(phi    [0],theta  [0]      );
		h2_e_phiVz0    -> Fill(phi    [0],targetZ[0]      );	
		h2_e_thetaVz0  -> Fill(theta  [0],targetZ[0]      );	
		h1_e_vz0       -> Fill(targetZ[0]                 );
		h1_e_vz        -> Fill(e_vz_corrected             );

		h1b->Fill(z_g    [0]);
		h1c->Fill(targetZ[0]);
		h2b->Fill(phi_g  [0],theta_g[0]);
		h2c->Fill(phi    [0],theta  [0]);
		h3b->Fill(theta_g[0],p_g    [0]);
		h3c->Fill(theta  [0],mom    [0]);

		for(int j = 0 ; j < 6 ; j++){
			if(p_g[0]>2.0 && p_g[0]<3.5 && phi_g[0]>phi_low[j] && phi_g[0]<phi_hi[j])
				h4b[j]->Fill(phi_g[0],theta_g[0]);
		}

		// --------------------------------------------------------------------------------------------------
		//Fill Generated Proton Histograms only if an electron is reconstructed, 
		//since we need a reconstructed electron for the proton timing,
		//and want to calculate the proton acceptance independently of the electron.
		if(num_g>1 && id_g[1]==2212){ //This will be always true when only an electon and proton are generated for every event
			hp1a->Fill(z_g[1]);
			hp2a->Fill(phi_g[1],theta_g[1]);
			hp3a->Fill(theta_g[1],p_g[1]);

			for(Int_t j=0;j<6;j++){
				if(p_g[1]>1.0 && p_g[1]<2.5 && phi_g[1]>phi_low[j] && phi_g[1]<phi_hi[j])
					hp4a[j]->Fill(phi_g[1],theta_g[1]);
			}
		}
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

			beta_assuming_pip = mom[i]/sqrt(mom[i]*mom[i] + mpc*mpc);
			pip_t0 = SC_Time[i] - SC_Path[i]/(beta_assuming_pip * c_cm_ns);
			pip_delta_t = pip_t0 - e_t0;
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
				h2_pip_deltaTmom0 -> Fill(pip_delta_t,mom  [i]);

				// Passing positive hadron fiducial cuts
				if(fid_params.pFiducialCut(T3_p_mom)){

					// Positive particle vertex (_z) correction
					p_vz_corrected = targetZ[i];

					h1_p_mass         -> Fill(mass[i]);
					h2_p_pMass        -> Fill(mom [i]    ,mass [i]);
					h2_pos_pBeta      -> Fill(mom [i]    ,Beta [i]);
					h2_p_phiTheta1    -> Fill(phi [i]    ,theta[i]);
					h2_p_deltaTmom1   -> Fill(delta_t    ,mom  [i]);
					h2_pip_deltaTmom1 -> Fill(pip_delta_t,mom  [i]);
					// --------------------------------------------------------------------
					// Look specifically for protons
					if((id_guess[i] == 2212 ) &&       // Guess at the particle ID is good for the proton candidate
							(fid_params.in_p_deltaT(delta_t, mom[i], pdeltat_sig_cutrange)) // Proton PID (delta T vs p)
					  ){
						p_mom_corrected=mom[i];	

						corr_px = p_mom_corrected*u1.X();
						corr_py = p_mom_corrected*u1.Y();
						corr_pz = p_mom_corrected*u1.Z();

						T3_p_mom.SetXYZ(corr_px,corr_py,corr_pz);	

						h2_p_deltaTmom2-> Fill(delta_t   ,mom                   [i]);
						h2_p_phiTheta2 -> Fill(phi    [i],theta                 [i]);	
						h2_p_phiVz0    -> Fill(phi    [i],targetZ               [i]);	
						h2_p_thetaVz0  -> Fill(theta  [i],targetZ               [i]);	
						h2_p_pBeta     -> Fill(mom    [i],Beta                  [i]);

						hp1b->Fill(z_g    [1]); //Careful w/ array elements here
						hp1c->Fill(targetZ[i]);
						hp2b->Fill(phi_g  [1],theta_g[1]);
						hp2c->Fill(phi    [i],theta  [i]);
						hp3b->Fill(theta_g[1],p_g    [1]);
						hp3c->Fill(theta_g[i],p_g    [i]);

						for(int k = 0 ; k < 6 ; k++){
							if(p_g[1]>1.0 && p_g[1]<2.5 && phi_g[1]>phi_low[k] && phi_g[1]<phi_hi[k])
								hp4b[k]->Fill(phi_g[i],theta_g[i]);
						}

						nProtons++;
						nParticles++;
					}
					// --------------------------------------------------------------------
					// Look specifically for pions +
					else if((id_guess[i] == 211)&&       // Guess at the particle ID is good for the pion+ candidate
							(fid_params.in_pip_deltaT(pip_delta_t, mom[i], pipdeltat_sig_cutrange)) // Pi+ PID
					       ){	

						h2_pip_pBeta      -> Fill(mom[i]          ,Beta [i]);	
						h2_pip_deltaTmom2 -> Fill(pip_delta_t     ,mom  [i]);

						nPiplus++;
						nParticles++;
					}
					// --------------------------------------------------------------------
				}
			}
			// ------------------------------------------------------------------------------------------
			// Test if negative particle
			else if(        (StatSC[i] > 0) &&              // SC status is good for the positive candidate
					(StatDC[i] > 0) &&              // DC status is good for the positive candidate
					(Stat  [i] > 0) &&              // Global status is good for the positive candidate
					(Charge[i] < 0)                 // Charge is negative
			       )
			{
				h2_pim_deltaTmom0 -> Fill(pip_delta_t,mom  [i]);

				// No fiducial cuts in place just yet
				h2_pim_deltaTmom1 -> Fill(pip_delta_t,mom  [i]);

				if((id_guess[i] == -211)&& // Guess at the particle ID is good for the pion- candidate
						(fid_params.in_pim_deltaT(pip_delta_t, mom[i], pipdeltat_sig_cutrange)) // Pi- PID
				  ){
					h2_pim_pBeta      -> Fill(mom[i]          ,Beta [i]);
					h2_pim_deltaTmom2 -> Fill(pip_delta_t     ,mom  [i]);

					nPiminus++;
					nParticles++;
				}
			}
		}

	}
	cerr << "Finished with the event loop...\n";

	gStyle->SetOptStat(0);
	gStyle->SetPaintTextFormat("3.2g");

	// --------------------------------------------------------------------------------------------------
	// Calculating the maps
	TH2F * h4r [6]; //ratio of accepted to generated
	TH2F * hp4r[6]; //ratio of accepted to generated

	for(int i = 0 ; i < 6 ; i++){
		h4r[i]  = (TH2F*)h4b[i]->Clone(Form("h4r[%d]", i));
		h4r[i] ->Divide(h4a[i]);
		hp4r[i] = (TH2F*)hp4b[i]->Clone(Form("hp4r[%d]", i));
		hp4r[i]->Divide(hp4a[i]);
	}

	TCanvas *c4[6];
	for(int i = 0 ; i < 6 ; i++){
		c4[i] = new TCanvas(Form("c4[%d]",i));
		h4r[i]->Draw("col text");
		c4[i]->Modified();c4[i]->Update();
	}

	TCanvas *c8[6];
        for(int i = 0 ; i < 6 ; i++){
                c8[i] = new TCanvas(Form("c8[%d]",i));
                hp4r[i]->Draw("col text");
                c8[i]->Modified();c8[i]->Update();
        }

	// --------------------------------------------------------------------------------------------------
	// Editing histograms
	h1_e_Nphe1  -> SetLineColor(2);
	h1_e_Nphe2  -> SetLineColor(8);
	h1_e_EC_in1 -> SetLineColor(2);
	h1_e_EC_in2 -> SetLineColor(8);
	h1_e_EC_out1-> SetLineColor(2);
	h1_e_EC_out2-> SetLineColor(8);
	h1_e_EC_tot1-> SetLineColor(2);
	h1_e_EC_tot2-> SetLineColor(8);
	h1_e_vz     -> SetLineColor(2);
	// ---
	TLegend * leg = new TLegend(0.5,0.5,0.8,0.8);
	leg -> AddEntry(h1_e_Nphe0,"before cuts"       );
	leg -> AddEntry(h1_e_Nphe1,"after PID cuts"    );
	leg -> AddEntry(h1_e_Nphe2,"after PID+Fiducial");
	// -------------------------------------------------------------------------------------------------- 
	TCanvas *c1 = new TCanvas("c1");
	c1->Divide(1,3);
	c1->cd(1);h1a->Draw();
	c1->cd(2);h1b->Draw();
	c1->cd(3);h1c->Draw();

	TCanvas *c2 = new TCanvas("c2");
	c2->Divide(1,3);
	c2->cd(1);h2a->Draw("colz");
	c2->cd(2);h2b->Draw("colz");
	c2->cd(3);h2c->Draw("colz");

	TCanvas *c3 = new TCanvas("c3");
	c3->Divide(1,3);
	c3->cd(1);h3a->Draw("colz");
	c3->cd(2);h3b->Draw("colz");
	c3->cd(3);h3c->Draw("colz");

	TCanvas *c5 = new TCanvas("c5");
	c5->Divide(1,3);
	c5->cd(1);hp1a->Draw();
	c5->cd(2);hp1b->Draw();
	c5->cd(3);hp1c->Draw();

	TCanvas *c6 = new TCanvas("c6");
	c6->Divide(1,3);
	c6->cd(1);hp2a->Draw("colz");
	c6->cd(2);hp2b->Draw("colz");
	c6->cd(3);hp2c->Draw("colz");

	TCanvas *c7 = new TCanvas("c7");
	c7->Divide(1,3);
	c7->cd(1);hp3a->Draw("colz");
	c7->cd(2);hp3b->Draw("colz");
	c7->cd(3);hp3c->Draw("colz");

	TCanvas *c9 = new TCanvas("c9");
	h1_e_Nphe0 -> Draw();
	h1_e_Nphe1 -> Draw("same");
	h1_e_Nphe2 -> Draw("same");
	leg        -> Draw("same");

	TCanvas *c10 = new TCanvas("c10");
	h1_e_EC_in0 -> Draw();
	h1_e_EC_in1 -> Draw("same");
	h1_e_EC_in2 -> Draw("same");
	leg         -> Draw("same");

	TCanvas *c11 = new TCanvas("c11");
	h1_e_EC_out0 -> Draw();
	h1_e_EC_out1 -> Draw("same");
	h1_e_EC_out2 -> Draw("same");
	leg          -> Draw("same");

	TCanvas *c12 = new TCanvas("c12");
	h1_e_EC_tot0 -> Draw();
	h1_e_EC_tot1 -> Draw("same");
	h1_e_EC_tot2 -> Draw("same");
	leg          -> Draw("same");

	TCanvas *c13 = new TCanvas("c13");
	c13 -> Divide(3,1);
	c13 -> cd(1);	h2_e_thetaMom0 -> Draw("COLZ"); 
	c13 -> cd(2);	h2_e_thetaMom1 -> Draw("COLZ");
	c13 -> cd(3);	h2_e_thetaMom2 -> Draw("COLZ");

	TCanvas *c14 = new TCanvas("c14");
	c14 -> Divide(3,1);
	c14 -> cd(1);	h2_e_Ein_Eout0 -> Draw("COLZ");
	c14 -> cd(2);	h2_e_Ein_Eout1 -> Draw("COLZ");
	c14 -> cd(3);	h2_e_Ein_Eout2 -> Draw("COLZ");

	TCanvas *c15 = new TCanvas("c15");
	c15 -> Divide(3,1);
	c15 -> cd(1);    h2_e_Ein_Eout_0 -> Draw("COLZ");
	c15 -> cd(2);    h2_e_Ein_Eout_1 -> Draw("COLZ");
	c15 -> cd(3);    h2_e_Ein_Eout_2 -> Draw("COLZ");

	TCanvas *c16 = new TCanvas("c16");
	c16 -> Divide(3,1);
	c16 -> cd(1);	h2_e_xyEC_hit0 -> Draw("COLZ");
	c16 -> cd(2);	h2_e_xyEC_hit1 -> Draw("COLZ");
	c16 -> cd(3);	h2_e_xyEC_hit2 -> Draw("COLZ");

	TCanvas *c17 = new TCanvas("c17");
	c17 -> Divide(3,1);
	c17 -> cd(1);	h2_e_p_Etot0 -> Draw("COLZ");
	c17 -> cd(2);	h2_e_p_Etot1 -> Draw("COLZ");
	c17 -> cd(3);	h2_e_p_Etot2 -> Draw("COLZ");

	TCanvas *c18 = new TCanvas("c18");
	c18 -> Divide(3,1);
	c18 -> cd(1);	h2_e_p_E0 -> Draw("COLZ");
	c18 -> cd(2);	h2_e_p_E1 -> Draw("COLZ");
	c18 -> cd(3);	h2_e_p_E2 -> Draw("COLZ");

	TCanvas *c19 = new TCanvas("c19");	
	h2_e_phiTheta0 -> Draw("COLZ");

	TCanvas *c20 = new TCanvas("c20");
	h2_e_phiTheta1 -> Draw("COLZ");
	
	TCanvas *c21 = new TCanvas("c21");
	h2_e_phiTheta2 -> Draw("COLZ");

	TCanvas *c22 = new TCanvas("c22");
	h1_e_vz0 -> Draw();
	h1_e_vz  -> Draw("same");

	TCanvas *c23 = new TCanvas("c23");
	h2_e_phiVz0 -> Draw("COLZ");

	TCanvas *c24 = new TCanvas("c24");
	h2_e_thetaVz0 -> Draw("COLZ");

	TCanvas *c25 = new TCanvas("c25");
	h2_pos_pBeta -> Draw("COLZ");

	TCanvas *c26 = new TCanvas("c26");
	h2_p_pBeta -> Draw("COLZ");

	TCanvas *c27 = new TCanvas("c27");
	c27 -> Divide(2,1);
	c27 -> cd(1);	h1_p_mass ->Draw();
	c27 -> cd(2);	h2_p_pMass->Draw("COLZ");

	TCanvas *c28 = new TCanvas("c28");
	h2_p_phiVz0 -> Draw("COLZ");

	TCanvas *c29 = new TCanvas("c29");
	h2_p_thetaVz0 -> Draw("COLZ");

	TCanvas *c30 = new TCanvas("c30");
	h2_p_deltaTmom0 -> Draw("COLZ");

	TCanvas *c31 = new TCanvas("c31");
	h2_p_deltaTmom1 -> Draw("COLZ");

	TCanvas *c32 = new TCanvas("c32");
	h2_p_deltaTmom2 -> Draw("COLZ");

	TCanvas *c33 = new TCanvas("c33");
	h2_p_phiTheta0 -> Draw("COLZ");

	TCanvas *c34 = new TCanvas("c34");
	h2_p_phiTheta1 -> Draw("COLZ");

	TCanvas *c35 = new TCanvas("c35");
	h2_p_phiTheta2 -> Draw("COLZ");

	TCanvas *c36 = new TCanvas("c36");
	h2_pip_deltaTmom0 -> Draw("COLZ");

	TCanvas *c37 = new TCanvas("c37");
	h2_pip_deltaTmom1 -> Draw("COLZ");

	TCanvas *c38 = new TCanvas("c38");
	h2_pip_deltaTmom2 -> Draw("COLZ");

	TCanvas *c39 = new TCanvas("c39");
	h2_pim_deltaTmom0 -> Draw("COLZ");

	TCanvas *c40 = new TCanvas("c40");
	h2_pim_deltaTmom1 -> Draw("COLZ");

	TCanvas *c41 = new TCanvas("c41");
	h2_pim_deltaTmom2 -> Draw("COLZ");

	// --------------------------------------------------------------------------------------------------
	// Print histograms on a pdf file

	c1  -> Print(Form("./e2a_maps_%d.pdf(",tab_E1) ,"pdf");
	c2  -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c3  -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");	
	for(int i = 0 ; i < 6 ; i++) c4[i]->Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c5  -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        c6  -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        c7  -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        for(int i = 0 ; i < 6 ; i++) c8[i]->Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c9  -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c10 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c11 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c12 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c13 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");   
	c14 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c15 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c16 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c17 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c18 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c19 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c20 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c21 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c22 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c23 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c24 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c25 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c26 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c27 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c28 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c29 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c30 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c31 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c32 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c33 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c34 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c35 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c36 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c37 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c38 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf"); 
	c39 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c40 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
	c41 -> Print(Form("./e2a_maps_%d.pdf)",tab_E1) ,"pdf");

	// --------------------------------------------------------------------------------------------------

	// --------------------------------------------------------------------------------------------------
	// Write the output file
	outfile->cd();	
	h1_e_EC_in0       ->Write();
	h1_e_EC_in1       ->Write();
	h1_e_EC_in2       ->Write();
	h1_e_EC_out0      ->Write();
	h1_e_EC_out1      ->Write();
	h1_e_EC_out2      ->Write();
	h1_e_EC_tot0      ->Write();
	h1_e_EC_tot1      ->Write();
	h1_e_EC_tot2      ->Write();
	h2_e_thetaMom0    ->Write();
	h2_e_thetaMom1    ->Write();
	h2_e_thetaMom2    ->Write();	
	h2_e_Ein_Eout0    ->Write();
	h2_e_Ein_Eout1    ->Write();
	h2_e_Ein_Eout2    ->Write();
	h2_e_xyEC_hit0    ->Write();
	h2_e_xyEC_hit1    ->Write();
	h2_e_xyEC_hit2    ->Write();
	h2_e_p_Etot0      ->Write();
	h2_e_p_Etot1      ->Write();
	h2_e_p_Etot2      ->Write();
	h2_e_p_E0         ->Write();
	h2_e_p_E1         ->Write();
	h2_e_p_E2         ->Write();
	h1_e_Nphe0        ->Write();
	h1_e_Nphe1        ->Write();
	h1_e_Nphe2        ->Write();
	h2_e_phiTheta0    ->Write();
	h2_e_phiTheta1    ->Write();
	h2_e_phiTheta2    ->Write();
	h1_e_vz0          ->Write();
	h1_e_vz           ->Write();
	h2_e_phiVz0       ->Write();	
	h2_e_thetaVz0     ->Write();	
	// ---
	h2_pos_pBeta      ->Write();
	h2_p_pBeta        ->Write();
	h2_pip_pBeta      ->Write();
	// ---
	h1_p_mass         ->Write();
	h2_p_pMass        ->Write();	
	h2_p_phiVz0       ->Write(); 
	h2_p_thetaVz0     ->Write();	
	h2_p_deltaTmom0   ->Write();
	h2_p_deltaTmom1   ->Write();
	h2_p_deltaTmom2   ->Write();
	h2_p_phiTheta0    ->Write();
	h2_p_phiTheta1    ->Write();
	h2_p_phiTheta2    ->Write();
	// ---
	h2_pip_deltaTmom0 ->Write();
	h2_pip_deltaTmom1 ->Write();
	h2_pip_deltaTmom2 ->Write();

	// Clean up
	f->Close();
	outfile->Close();

	return 0;
}
