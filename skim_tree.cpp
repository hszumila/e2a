#include <fstream>
#include <iostream>
#include <cmath>

#include"TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include <TLegend.h>

#include "Fiducial.h"
#include "Run_dependent.h"
#include "constants.h"

using namespace std;

// Range of momentum for which we have good fits for PID
const double epratio_sig_cutrange=3.; // +/- sigma for PID cut
const double pdeltat_sig_cutrange=3.;
const double pipdeltat_sig_cutrange=1.5;

// Difference between positive corrected z vertex and electron corrected z vertex
const double pos_z_cut_min = -2;
const double pos_z_cut_max =  2;

// Parameters for 2GeV data cuts
const double sc_cc_dt_cut_sect[6]={-2,-5,-8,-8,-2,2};

int main(int argc, char ** argv)
{
	if (argc != 3)
	{
		cerr << "Wrong number of arguments. Instead try\n"
			<< "\tskim_tree /path/to/input/file /path/to/output/file\n";
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
	// --------------------------------------------------------------------------------------------------
	// Open up the tree, and get the important data
	TTree * t = (TTree*)f->Get("data");
	const int nEvents = t->GetEntries();
	const int maxPart = 50;
	int gPart, CCPart, DCPart, ECPart, SCPart, NRun;
	int StatDC[maxPart], StatCC[maxPart], StatEC[maxPart], StatSC[maxPart], id_guess[maxPart];
	float Xb, STT, Q2, W, Nu, Yb;
	float Stat[maxPart], EC_in[maxPart], EC_out[maxPart], EC_tot[maxPart], Nphe[maxPart],
	      SC_Time[maxPart], SC_Path[maxPart], CC_Time[maxPart], CC_Path[maxPart],
	      charge[maxPart], beta[maxPart], mass[maxPart], mom[maxPart], px[maxPart], py[maxPart],
	      pz[maxPart], theta[maxPart], phi[maxPart], targetZ[maxPart], theta_pq[maxPart],
	      EC_X[maxPart],EC_Y[maxPart],EC_Z[maxPart], CC_Chi2[maxPart];
	t->SetBranchAddress("NRun"     ,&NRun   ); // Run number
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
	t->SetBranchAddress("Charge"   ,charge  ); // Charge per particle
	t->SetBranchAddress("Beta"     ,beta    ); // Beta per particle
	t->SetBranchAddress("Mass"     ,mass    ); // Mass per particle
	t->SetBranchAddress("Momentum" ,mom     ); // Momentum magnitude per particle
	t->SetBranchAddress("Momentumx",px      ); // Momentum x component per particle
	t->SetBranchAddress("Momentumy",py      ); // Momentum y component per particle
	t->SetBranchAddress("Momentumz",pz      ); // Momentum z component per particle
	t->SetBranchAddress("Theta"    ,theta   ); // Theta per particle
	t->SetBranchAddress("Phi"      ,phi     ); // Phi per particle
	t->SetBranchAddress("TargetZ"  ,targetZ ); // Target Z per particle
	t->SetBranchAddress("Thetapq"  ,theta_pq); // Angle wrt to q vector per particle
	t->SetBranchAddress("STT"      ,&STT    ); // RF-corrected start time.
	t->SetBranchAddress("EC_X"     ,EC_X    ); // x positions of hit in the calorimeter
	t->SetBranchAddress("EC_Y"     ,EC_Y    ); // y positions of hit in the calorimeter
	t->SetBranchAddress("EC_Z"     ,EC_Z    ); // z positions of hit in the calorimeter
	t->SetBranchAddress("CC_Chi2"  ,CC_Chi2 ); // angle between CC hit and nearest SC hit (in rad)
	//t->SetBranchAddress("BjorkenX" ,&Xb     ); // Bjorken X
	//t->SetBranchAddress("Q2"       ,&Q2     ); // Momentum transfer
	//t->SetBranchAddress("W"        ,&W      ); // Hadronic mass
	//t->SetBranchAddress("Nu"       ,&Nu     ); // Energy transfer
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
	TH1D * hist_e_Nphe0    = new TH1D("hist_e_Nphe0"    ,"e- before cuts;# photo-electrons in CC;Counts"     ,200,   0.,200.);
	TH1D * hist_e_Nphe1    = new TH1D("hist_e_Nphe1"    ,"e- passing PID;# photo-electrons in CC;Counts"     ,200,   0.,200.);
	TH1D * hist_e_Nphe2    = new TH1D("hist_e_Nphe2"    ,"e- passing PID+fid;# photo-electrons in CC;Counts" ,200,   0.,200.);

	TH1D * hist_e_EC_in0   = new TH1D("hist_e_EC_in0"   ,"e- before cuts;E_in [GeV];Counts"                  ,300,   0.,  1.);
	TH1D * hist_e_EC_in1   = new TH1D("hist_e_EC_in1"   ,"e- passing PID;E_in [GeV];Counts"                  ,300,   0.,  1.);
	TH1D * hist_e_EC_in2   = new TH1D("hist_e_EC_in2"   ,"e- passing PID+fid;E_in [GeV];Counts"              ,300,   0.,  1.);

	TH1D * hist_e_EC_out0  = new TH1D("hist_e_EC_out0"  ,"e- before cuts;E_out [GeV];Counts"                 ,300,   0.,  1.);
	TH1D * hist_e_EC_out1  = new TH1D("hist_e_EC_out1"  ,"e- passing PID;E_out [GeV];Counts"                 ,300,   0.,  1.);
	TH1D * hist_e_EC_out2  = new TH1D("hist_e_EC_out2"  ,"e- passing PID+fid;E_out [GeV];Counts"             ,300,   0.,  1.);

	TH1D * hist_e_EC_tot0   = new TH1D("hist_e_EC_tot0"  ,"e- before cuts;E_tot [GeV];Counts"                ,300,   0.,  1.);
	TH1D * hist_e_EC_tot1   = new TH1D("hist_e_EC_tot1"  ,"e- passing PID;E_tot [GeV];Counts"                ,300,   0.,  1.);
	TH1D * hist_e_EC_tot2   = new TH1D("hist_e_EC_tot2"  ,"e- passing PID+fid;E_tot [GeV];Counts"            ,300,   0.,  1.);

	TH2D * hist_e_phiTheta0= new TH2D("hist_e_phiTheta0","e- before  cuts;#phi [deg];#theta [deg];Counts"      ,300,-100.,380.,300,10.,50.);
	TH2D * hist_e_phiTheta1= new TH2D("hist_e_phiTheta1","e- passing PID;#phi [deg];#theta [deg];Counts"       ,300,-100.,380.,300,10.,50.);
	TH2D * hist_e_phiTheta2= new TH2D("hist_e_phiTheta2","e- passing PID+fid;#phi [deg];#theta [deg];Counts"   ,300,-100.,380.,300,10.,50.);

	TH2D * hist_e_Ein_Eout0= new TH2D("hist_e_Ein_Eout0","e- before cuts;E_in/p;E_out/p;Counts"              ,300,   0., 0.5,300, 0.,0.5);
	TH2D * hist_e_Ein_Eout1= new TH2D("hist_e_Ein_Eout1","e- passing PID;E_in/p;E_out/p;Counts"              ,300,   0., 0.5,300, 0.,0.5);
	TH2D * hist_e_Ein_Eout2= new TH2D("hist_e_Ein_Eout2","e- passing PID+fid;E_in/p;E_out/p;Counts"          ,300,   0., 0.5,300, 0.,0.5);

	TH2D * hist_e_xyEC_hit0= new TH2D("hist_e_xyEC_hit0","e- before cuts;ECx [cm];ECy [cm];Counts"           ,300,-400.,400.,300,-400.,400.);
	TH2D * hist_e_xyEC_hit1= new TH2D("hist_e_xyEC_hit1","e- passing PID;ECx [cm];ECy [cm];Counts"           ,300,-400.,400.,300,-400.,400.);
	TH2D * hist_e_xyEC_hit2= new TH2D("hist_e_xyEC_hit2","e- passing PID+fid;ECx [cm];ECy [cm];Counts"       ,300,-400.,400.,300,-400.,400.);

	TH2D * hist_e_p_Etot0  = new TH2D("hist_e_p_Etot0"  ,"e- before cuts;p [GeV];E_tot/p;Counts"             ,300,   0.,  5.,300, 0.,0.7);
	TH2D * hist_e_p_Etot1  = new TH2D("hist_e_p_Etot1"  ,"e- passing PID;p [GeV];E_tot/p;Counts"             ,300,   0.,  5.,300, 0.,0.7);
	TH2D * hist_e_p_Etot2  = new TH2D("hist_e_p_Etot2"  ,"e- passing PID+fid;p [GeV];E_tot/p;Counts"         ,300,   0.,  5.,300, 0.,0.7);

	TH2D * hist_e_p_E0     = new TH2D("hist_e_p_E0"     ,"e- before cuts;p [GeV];E;Counts"                   ,300,   0.,  5.,300, 0.,2.);
	TH2D * hist_e_p_E1     = new TH2D("hist_e_p_E1"     ,"e- passing PID;p [GeV];E;Counts"                   ,300,   0.,  5.,300, 0.,2.);
	TH2D * hist_e_p_E2     = new TH2D("hist_e_p_E2"     ,"e- before PIF+fid;p [GeV];E;Counts"                ,300,   0.,  5.,300, 0.,2.);

	TH2D * hist_e_thetaMom0= new TH2D("e_thetaMom0"     ,"e- before cuts;#theta [deg];Mom [GeV];Counts"       ,300,  10., 50.,300, 0., 6.);
	TH2D * hist_e_thetaMom1= new TH2D("e_thetaMom1"     ,"e- passing PID;#theta [deg];Mom [GeV];Counts"       ,300,  10., 50.,300, 0., 6.);
	TH2D * hist_e_thetaMom2= new TH2D("e_thetaMom2"     ,"e- passing PID+fid;#theta [deg];Mom [GeV];Counts"   ,300,  10., 50.,300, 0., 6.);

	// ---
	TH1D * hist_e_vz_sec10 = new TH1D("hist_e_vz_sec10" ,"e- passing cuts, before vtx corr, sector 1;electron vz [cm]; Counts" ,300, -10., 10.);
	TH1D * hist_e_vz_sec1  = new TH1D("hist_e_vz_sec1"  ,"e- passing cuts,  after vtx corr, sector 1;electron vz [cm]; Counts" ,300, -10., 10.);
	TH1D * hist_e_vz_sec20 = new TH1D("hist_e_vz_sec20" ,"e- passing cuts, before vtx corr, sector 2;electron vz [cm]; Counts" ,300, -10., 10.);
	TH1D * hist_e_vz_sec2  = new TH1D("hist_e_vz_sec2"  ,"e- passing cuts,  after vtx corr, sector 2;electron vz [cm]; Counts" ,300, -10., 10.);
	TH1D * hist_e_vz_sec30 = new TH1D("hist_e_vz_sec30" ,"e- passing cuts, before vtx corr, sector 3;electron vz [cm]; Counts" ,300, -10., 10.);
	TH1D * hist_e_vz_sec3  = new TH1D("hist_e_vz_sec3"  ,"e- passing cuts,  after vtx corr, sector 3;electron vz [cm]; Counts" ,300, -10., 10.);
	TH1D * hist_e_vz_sec40 = new TH1D("hist_e_vz_sec40" ,"e- passing cuts, before vtx corr, sector 4;electron vz [cm]; Counts" ,300, -10., 10.);
	TH1D * hist_e_vz_sec4  = new TH1D("hist_e_vz_sec4"  ,"e- passing cuts,  after vtx corr, sector 4;electron vz [cm]; Counts" ,300, -10., 10.);
	TH1D * hist_e_vz_sec50 = new TH1D("hist_e_vz_sec50" ,"e- passing cuts, before vtx corr, sector 5;electron vz [cm]; Counts" ,300, -10., 10.);
	TH1D * hist_e_vz_sec5  = new TH1D("hist_e_vz_sec5"  ,"e- passing cuts,  after vtx corr, sector 5;electron vz [cm]; Counts" ,300, -10., 10.);
	TH1D * hist_e_vz_sec60 = new TH1D("hist_e_vz_sec60" ,"e- passing cuts, before vtx corr, sector 6;electron vz [cm]; Counts" ,300, -10., 10.);
	TH1D * hist_e_vz_sec6  = new TH1D("hist_e_vz_sec6"  ,"e- passing cuts,  after vtx corr, sector 6;electron vz [cm]; Counts" ,300, -10., 10.);
	TH1D * hist_e_vz0      = new TH1D("hist_e_vz0"      ,"e- passing cuts, before vtx corr;electron vz [cm]; Counts"           ,300, -10., 10.);
	TH1D * hist_e_vz       = new TH1D("hist_e_vz"       ,"e- passing cuts,  after vtx corr;electron vz [cm]; Counts"           ,300, -10., 10.);
	TH2D * hist_e_phiVz0   = new TH2D("hist_e_phiVz0"   ,"e- passing cuts, before vtx corr;#phi [deg];vz [cm];Counts"   ,300,-100.,380.,300,-10,10);
	TH2D * hist_e_phiVz    = new TH2D("hist_e_phiVz"    ,"e- passing cuts,  after vtx corr;#phi [deg];vz [cm];Counts"   ,300,-100.,380.,300,-10,10);
	TH2D * hist_e_thetaVz0 = new TH2D("hist_e_thetaVz0" ,"e- passing cuts, before vtx corr;#theta [deg];vz [cm];Counts" ,300, -10., 60.,300,-11,11);
	TH2D * hist_e_thetaVz  = new TH2D("hist_e_thetaVz"  ,"e- passing cuts,  after vtx corr;#theta [deg];vz [cm];Counts" ,300, -10., 60.,300,-11,11);
	// ---

	TH1D * hist_e_momCor   = new TH1D("hist_e_momCor"   ,"e- passing fid. cuts;p corrected - p[GeV];Counts"          ,300, -.1,.04);
	TH2D * hist_e_momMomCor= new TH2D("hist_e_momMomCor","e- passing fid. cuts;p [GeV];p corrected - p[GeV];Counts"  ,300,  0.,  6.,300,-.1,.04);
	TH2D * hist_e_vzVzCor  = new TH2D("hist_e_vzVzCor"  ,"e- passing fid. cuts;vz [cm];vz corrected - vz [cm];Counts",300,-20., 20.,300,-1., 1.);

	// ---------------------------------------
	// Diagnostic positive particle histograms
	TH1D * hist_p_mass     = new TH1D("hist_pos_mass"   ,"+  passing fid. cuts;mass [GeV];Counts"          ,300,  0.,3.5 );
	TH2D * hist_p_pMass    = new TH2D("hist_pos_pMass"  ,"+  passing fid. cuts;p [GeV];mass [GeV];Counts"  ,300,  0., 4.,300, 0.,3.5);
	TH2D * hist_pos_pBeta  = new TH2D("hist_pos_pBeta"  ,"+  passing fid. cuts;p [GeV];#beta;Counts"       ,300,  0., 4.,300, 0.,1.3);

	// ---------------------------------------
	// Diagnostic proton histograms
	TH2D * hist_p_phiTheta0= new TH2D("hist_p_phiTheta0","p before cuts;#phi [deg];#theta [deg];Counts"    ,300,-100.,380.,300,0.,55.);
	TH2D * hist_p_phiTheta1= new TH2D("hist_p_phiTheta1","p passing fid;#phi [deg];#theta [deg];Counts"    ,300,-100.,380.,300,0.,55.);
	TH2D * hist_p_phiTheta2= new TH2D("hist_p_phiTheta2","p passing fid+PID;#phi [deg];#theta [deg];Counts",300,-100.,380.,300,0.,55.);

	TH2D * hist_p_deltaTmom0=new TH2D("hist_p_deltaTmom0","p before cuts;#Delta t [ns];p [GeV];Counts"     ,300,  -7.,  7.,300, 0., 3.);
	TH2D * hist_p_deltaTmom1=new TH2D("hist_p_deltaTmom1","p passing fid;#Delta t [ns];p [GeV];Counts"     ,300,  -7.,  7.,300, 0., 3.);
	TH2D * hist_p_deltaTmom2=new TH2D("hist_p_deltaTmom2","p passing fid+PID;#Delta t [ns];p [GeV];Counts" ,300,  -7.,  7.,300, 0., 3.);

	TH2D * hist_p_p_momCor = new TH2D("hist_p_p_momCor" ,"p passing fid. cuts;p [GeV];p - p_corr [GeV];Counts" ,100,   0.,  3.,100,-.1,.01);
	TH2D * hist_p_vzVzCor  = new TH2D("hist_p_vzVzCor"  ,"p passing fid. cuts;vz [cm];vz corrected - vz [cm];Counts",300, -20., 20.,300,-1., 1 );
	TH2D * hist_p_phiVz0   = new TH2D("hist_p_phiVz0"   ,"p passing cuts, before vtx corr;#phi [deg];vz [cm];Counts",300,-100.,380.,300,-10,10 );
	TH2D * hist_p_phiVz    = new TH2D("hist_p_phiVz"    ,"p passing cuts,  after vtx corr;#phi [deg];vz [cm];Counts",300,-100.,380.,300,-10,10 );
	TH2D * hist_p_thetaVz0 = new TH2D("hist_p_thetaVz0" ,"p passing cuts, before vtx corr;#theta [deg];vz [cm];Counts",300, 0.,100.,300,-11,11 );
	TH2D * hist_p_thetaVz  = new TH2D("hist_p_thetaVz"  ,"p passing cuts,  after vtx corr;#theta [deg];vz [cm];Counts",300, 0.,100.,300,-11,11 );
	TH2D * hist_p_pBeta    = new TH2D("hist_p_pBeta"    ,"p passing fid. cuts;p [GeV];#beta;Counts"                 ,300,   0.,  4.,300, 0.,1.3);

	// ---------------------------------------
	// Diagnostic pi+ histograms
	TH2D * hist_pip_pBeta     = new TH2D("hist_pip_pBeta"     ,"pi+ passing fid. cuts;p [GeV];#beta;Counts"       ,300,   0.,  4.,300, 0.,1.3);
	TH2D * hist_pip_deltaTmom0= new TH2D("hist_pip_deltaTmom0","#pi+ before cuts;#Delta t [ns];p [GeV];Counts"    ,300,  -7.,  7.,300, 0.,5.0);
	TH2D * hist_pip_deltaTmom1= new TH2D("hist_pip_deltaTmom1","#pi+ passing fid;#Delta t [ns];p [GeV];Counts"    ,300,  -7.,  7.,300, 0.,5.0);
	TH2D * hist_pip_deltaTmom2= new TH2D("hist_pip_deltaTmom2","#pi+ passing fid+PID;#Delta t [ns];p [GeV];Counts",300,  -7.,  7.,300, 0.,5.0);

	// ---------------------------------------
	// Diagnostic neutral histograms
	TH2D * hist_neu_pBeta     = new TH2D("hist_neu_pBeta"     ,"neutral passing fid. cuts;p [GeV];#beta;Counts"   ,300,   0.,  4.,300, 0.,1.3);

	// ---------------------------------------
	// Diagnostic neutron histograms
	TH2D * hist_n_phiTheta0= new TH2D("hist_n_phiTheta0","n before cuts;#phi [deg];#theta [deg];Counts"     ,300,-100.,380.,300,0.,55.);
	TH2D * hist_n_phiTheta1= new TH2D("hist_n_phiTheta1","n passing fid;#phi [deg];#theta [deg];Counts"     ,300,-100.,380.,300,0.,55.);
	TH2D * hist_n_phiTheta2= new TH2D("hist_n_phiTheta2","n passing fid+PID;#phi [deg];#theta [deg];Counts" ,300,-100.,380.,300,0.,55.);

	TH2D * hist_n_pBeta    = new TH2D("hist_n_pBeta"    ,"n passing fid. cuts;p [GeV];#beta;Counts"         ,300,   0.,  4.,300, 0.,1.3);

	// ---------------------------------------
	// Setting up output tree and branches
	TTree * outtree = new TTree("T","Skimmed tree");
	double e_vz, e_vz_corrected, e_mom[3], e_phi_mod;
	double p_vz, p_vz_corrected, p_mom_corrected, p_phi_mod;
	double EC_in_cut, el_EC_cut;
	double e_t0,beta_assuming_proton,p_t0,delta_t,beta_assuming_pip,pip_t0,pip_delta_t;

	TVector3 e_ec_xyz;
	TVector3 T3_e_mom, T3_e_mom_cor, T3_p_mom;

	int nParticles;
	int nProtons, nNeutrons, nPiplus, nPiminus, nPi0;
	int Part_type    [maxPart];
	double vtx_z_unc [maxPart], vtx_z_cor[maxPart];
	double mom_x     [maxPart], mom_y    [maxPart], mom_z    [maxPart]; 
	double e_deltat  [maxPart];

	/* 
	   =========================
	   Part_type values	(http://www.star.bnl.gov/public/comp/simu/newsite/gstar/kumacs/NewParticle.html)
	   -11  = electron
	   2212 = proton
	   2112 = neutron
	   211  = pi+
	   -211 = pi-
	   111  = pi0
	   =========================
	 */ 

	outtree->Branch("nParticles", &nParticles, "nParticles/I"            );
	outtree->Branch("nProtons"  , &nProtons  , "nProtons/I"              );
	outtree->Branch("nNeutrons" , &nNeutrons , "nNeutrons/I"             );
	outtree->Branch("nPiplus"   , &nPiplus   , "nPiplus/I"               );
	outtree->Branch("Part_type" ,  Part_type , "Part_type[nParticles]/I" );
	//outtree->Branch("e_vz"      , &e_vz      , "e_vz/D"                  );
	outtree->Branch("vtx_z_unc" ,  vtx_z_unc , "vtx_z_unc[nParticles]/D" );
	outtree->Branch("vtx_z_cor" ,  vtx_z_cor , "vtx_z_cor[nParticles]/D" );
	//outtree->Branch("e_mom"     , e_mom      , "e_mom[3]/D"             );
	outtree->Branch("mom_x"     ,  mom_x     , "mom_x[nParticles]/D"     );
	outtree->Branch("mom_y"     ,  mom_y     , "mom_y[nParticles]/D"     );
	outtree->Branch("mom_z"     ,  mom_z     , "mom_z[nParticles]/D"     );
	outtree->Branch("e_deltat"  ,  e_deltat  , "e_deltat[nParticles]/D"  );

	// --------------------------------------------------------------------------------------------------
	// Obtaining run number and other important parameters
	t->GetEvent(0);
	int tab_run, tab_E1, tab_torus, tab_mini;
	string tab_targ;

	char param_file_name[256];
	string homedir = string(getenv("HOME"));
	sprintf(param_file_name,"%s/.e2a/run_table.dat",homedir.c_str());
	ifstream run_table;
	run_table.open(param_file_name);
	do{
		run_table >> tab_run  ;
		run_table >> tab_E1   ;
		run_table >> tab_torus;
		run_table >> tab_mini ;
		run_table >> tab_targ ;
	}
	while(tab_run != NRun);
	cout << "Run    = " << tab_run   << endl;
	cout << "Ebeam  = " << tab_E1    << endl;
	cout << "Torus  = " << tab_torus << endl;
	cout << "Mini   = " << tab_mini  << endl;
	cout << "Target = " << tab_targ  << endl;

	Fiducial fid_params(tab_E1,tab_torus,tab_mini,tab_targ);  // Create an instance of the Fiducial Class
	Run_dependent run_dependent_corrections(NRun);       // Create an instance of the Run_dependent Class

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
		e_vz_corrected = targetZ[0]+fid_params.vz_corr(T3_e_mom);
		e_ec_xyz.SetXYZ(EC_X[0],EC_Y[0],EC_Z[0]);

		// ---------------------------------------------------------------------------------------
		// Electron general cuts
		if (!(			(StatEC[0] > 0) && // EC status is good for the electron candidate
					(StatDC[0] > 0) && // DC status is good for the electron candidate
					(StatCC[0] > 0) && // CC status is good for the electron candidate
					(StatSC[0] > 0) && // SC status is good for the electron candidate
					(charge[0] < 0)    // Electron candidate curvature direction is negative
		     ))
		{continue;}
		// ---------------------------------------------------------------------------------------

		hist_e_Nphe0     -> Fill( Nphe        [0]  );
		hist_e_EC_in0    -> Fill( EC_in       [0]  );
		hist_e_EC_out0   -> Fill( EC_out      [0]  );
		hist_e_EC_tot0   -> Fill( EC_tot      [0]  );
		hist_e_phiTheta0 -> Fill( phi         [0], theta        [0]);
		hist_e_Ein_Eout0 -> Fill( EC_in[0]/mom[0], EC_out[0]/mom[0]);
		hist_e_p_Etot0   -> Fill( mom         [0], EC_tot[0]/mom[0]);
		hist_e_xyEC_hit0 -> Fill( EC_X        [0], EC_Y         [0]);
		hist_e_p_E0      -> Fill( mom         [0], el_cand_EC      );
		hist_e_thetaMom0 -> Fill( theta       [0], mom          [0]);

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
		hist_e_Nphe1     -> Fill( Nphe        [0]  );
		hist_e_EC_in1    -> Fill( EC_in       [0]  );
		hist_e_EC_out1   -> Fill( EC_out      [0]  );
		hist_e_EC_tot1   -> Fill( EC_tot      [0]  );
		hist_e_phiTheta1 -> Fill( phi         [0], theta        [0]);
		hist_e_Ein_Eout1 -> Fill( EC_in[0]/mom[0], EC_out[0]/mom[0]);
		hist_e_p_Etot1   -> Fill( mom         [0], EC_tot[0]/mom[0]);
		hist_e_xyEC_hit1 -> Fill( EC_X        [0], EC_Y         [0]);
		hist_e_p_E1      -> Fill( mom         [0], el_cand_EC      );
		hist_e_thetaMom1 -> Fill( theta       [0], mom          [0]);

		// ---------------------------------------------------------------------------------------
		// Electron Fiducial cuts
		if (!fid_params.e_inFidRegion(T3_e_mom)) continue; // Electron theta-phi cut
		if (!fid_params.CutUVW(e_ec_xyz)       ) continue; // Cuts on edges of calorimeter (u>60, v<360, w<400);

		// ---------------------------------------------------------------------------------------
		// If electron passes all cuts, then momentum-correct it (only works for theta > 16 deg):
		if (            ((tab_E1==4461)&&(180./M_PI*T3_e_mom.Theta()>16.))||
				( tab_E1==2261)
		   )
			T3_e_mom_cor = fid_params.eMomentumCorrection(T3_e_mom);
		else 	T3_e_mom_cor = T3_e_mom;

		nParticles++;
		Part_type[0] = -11;
		e_deltat [0] =   0;
		mom_x    [0] = T3_e_mom_cor.X();
		mom_y    [0] = T3_e_mom_cor.Y();
		mom_z    [0] = T3_e_mom_cor.Z();
		vtx_z_unc[0] = targetZ[0];
		vtx_z_cor[0] = e_vz_corrected;

		// If we get to here, then the electron passed fiducial cuts
		// Fill some diagnostic histograms
		hist_e_Nphe2    -> Fill(Nphe  [0]);
		hist_e_EC_in2   -> Fill(EC_in [0]);
		hist_e_EC_out2  -> Fill(EC_out[0]);
		hist_e_EC_tot2  -> Fill(EC_tot[0]);
		hist_e_xyEC_hit2-> Fill(EC_X[0],EC_Y[0]);
		hist_e_thetaMom2-> Fill(theta[0],mom[0]);
		hist_e_momCor   -> Fill(T3_e_mom_cor.Mag()-T3_e_mom.Mag());
		hist_e_momMomCor-> Fill(T3_e_mom.Mag(),T3_e_mom_cor.Mag()-T3_e_mom.Mag());
		hist_e_vzVzCor  -> Fill(targetZ[0],e_vz_corrected-targetZ[0]);
		hist_e_Ein_Eout2-> Fill(EC_in[0]/mom[0],EC_out[0]/mom[0]);
		hist_e_p_Etot2  -> Fill(mom    [0],EC_tot[0]/mom[0]);
		hist_e_p_E2     -> Fill(mom    [0],el_cand_EC      );
		hist_e_phiTheta2-> Fill(phi    [0],theta  [0]      );
		hist_e_phiVz0   -> Fill(phi    [0],targetZ[0]      );
		hist_e_phiVz    -> Fill(phi    [0],e_vz_corrected  );
		hist_e_thetaVz0 -> Fill(theta  [0],targetZ[0]      );
		hist_e_thetaVz  -> Fill(theta  [0],e_vz_corrected  );
		hist_e_vz0      -> Fill(targetZ[0]                 );
		hist_e_vz       -> Fill(e_vz_corrected             );
		if      (e_sect==0) {hist_e_vz_sec10 -> Fill(targetZ[0]);	hist_e_vz_sec1 -> Fill(e_vz_corrected);}
		else if (e_sect==1) {hist_e_vz_sec20 -> Fill(targetZ[0]);	hist_e_vz_sec2 -> Fill(e_vz_corrected);}
		else if (e_sect==2) {hist_e_vz_sec30 -> Fill(targetZ[0]);	hist_e_vz_sec3 -> Fill(e_vz_corrected);}
		else if (e_sect==3) {hist_e_vz_sec40 -> Fill(targetZ[0]);	hist_e_vz_sec4 -> Fill(e_vz_corrected);}
		else if (e_sect==4) {hist_e_vz_sec50 -> Fill(targetZ[0]);	hist_e_vz_sec5 -> Fill(e_vz_corrected);}
		else if (e_sect==5) {hist_e_vz_sec60 -> Fill(targetZ[0]);	hist_e_vz_sec6 -> Fill(e_vz_corrected);}
		else {cout << "Something is wrong with the definition of sectors" << endl;}

		// --------------------------------------------------------------------------------------------------
		// Loop over events looking for other particles
		for (int i=1 ; i<gPart ; i++)
		{
			Part_type[i] = 0;

			T3_p_mom.SetXYZ(px[i],py[i],pz[i]);
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
					(charge[i] > 0) 		// Charge is positive
			  )
			{
				hist_p_phiTheta0    -> Fill(phi[i]     ,theta[i]);
				hist_p_deltaTmom0   -> Fill(delta_t    ,mom  [i]); 
				hist_pip_deltaTmom0 -> Fill(pip_delta_t,mom  [i]);

				// Passing positive hadron fiducial cuts
				if(fid_params.pFiducialCut(T3_p_mom)){
					// Positive particle vertex (_z) correction
					p_vz_corrected = targetZ[i]+fid_params.vz_corr(T3_p_mom);

					hist_p_mass         -> Fill(mass[i]);
					hist_p_pMass        -> Fill(mom [i]    ,mass [i]);
					hist_pos_pBeta      -> Fill(mom [i]    ,beta [i]);
					hist_p_phiTheta1    -> Fill(phi [i]    ,theta[i]);
					hist_p_deltaTmom1   -> Fill(delta_t    ,mom  [i]);
					hist_pip_deltaTmom1 -> Fill(pip_delta_t,mom  [i]);
					// --------------------------------------------------------------------
					// Look specifically for protons
					if((id_guess[i] == 2212 ) &&       // Guess at the particle ID is good for the proton candidate
							(fid_params.in_p_deltaT(delta_t, mom[i], pdeltat_sig_cutrange)) // Proton PID (delta T vs p)
					  ){
						if(run_dependent_corrections.ProtonMomCorrection_He3_4Cell(T3_p_mom,p_vz_corrected) != -1)
							p_mom_corrected=run_dependent_corrections.ProtonMomCorrection_He3_4Cell(T3_p_mom,p_vz_corrected);
						else	p_mom_corrected=mom[i];	
		
						Part_type[nParticles] = 2212;
						e_deltat [nParticles] = delta_t;
						mom_x    [nParticles] = T3_p_mom.X();
						mom_y    [nParticles] = T3_p_mom.Y();
						mom_z    [nParticles] = T3_p_mom.Z();
						vtx_z_unc[nParticles] = targetZ  [i];	
						vtx_z_cor[nParticles] = p_vz_corrected;

						hist_p_deltaTmom2-> Fill(delta_t            ,mom [i]);
						hist_p_phiTheta2 -> Fill(phi    [i],theta[i]);
						hist_p_vzVzCor   -> Fill(targetZ[i],p_vz_corrected-targetZ[i]);
						hist_p_p_momCor  -> Fill(mom    [i],mom [i]-p_mom_corrected  );
						hist_p_phiVz0    -> Fill(phi    [i],targetZ[i]               );
						hist_p_phiVz     -> Fill(phi    [i],p_vz_corrected           );
						hist_p_thetaVz0  -> Fill(theta  [i],targetZ[i]               );
						hist_p_thetaVz   -> Fill(theta  [i],p_vz_corrected        	   );
						hist_p_pBeta     -> Fill(mom    [i],beta [i]                 );
						
						nProtons++;
						nParticles++;
					}
					// --------------------------------------------------------------------
					// Look specifically for pions +
					else if(id_guess[i] == 211       // Guess at the particle ID is good for the pion+ candidate
							// PID cuts for pions here
					       ){
						Part_type [nParticles] = 211;
						e_deltat  [nParticles] = pip_delta_t;
						mom_x     [nParticles] = T3_p_mom.X();
						mom_y     [nParticles] = T3_p_mom.Y();
						mom_z     [nParticles] = T3_p_mom.Z();
						vtx_z_unc [nParticles] = targetZ  [i];
						vtx_z_cor [nParticles] = p_vz_corrected;

						hist_pip_pBeta      -> Fill(mom[i] ,beta [i]);	
						hist_pip_deltaTmom2 -> Fill(pip_delta_t     ,mom  [i]);
					
						nPiplus++;
						nParticles++;
					}

					// --------------------------------------------------------------------
				}
			}
			// ------------------------------------------------------------------------------------------
			// Test if neutral particle
			else if(        StatEC[i] > 0 && // EC status is good for the positive candidate
					StatDC[i] <=0 && // DC status is good for the positive candidate
					StatSC[i] <=0 && // SC status is good for the positive candidate
					Stat  [i] > 0 && // Global status is good for the positive candidate
					charge[i] ==0    // Charge is neutral
			       )
			{
				hist_n_phiTheta0 -> Fill(phi[i],theta[i]);
				if(fid_params.pFiducialCut(T3_p_mom)){

					hist_neu_pBeta   -> Fill(mom [i],beta [i]);
					hist_n_phiTheta1 -> Fill(phi [i],theta[i]);

					// --------------------------------------------------------------------
					// Look specifically for neutrons 
					if(             beta[i] < 0.95 &&
							id_guess[i] == 2112 // Guess at the particle ID is good for the neutron candidate
					  )
					{

						Part_type[nParticles] = 2112;

						hist_n_phiTheta2 -> Fill(phi[i],theta[i]);
						hist_n_pBeta     -> Fill(mom[i],beta [i]);
						
						nNeutrons++;
						nParticles++;
					
					}

				}

			}
			// ------------------------------------------------------------------------------------------
			// Test if negative particle
			else if(        (StatSC[i] > 0) &&              // SC status is good for the positive candidate
					(StatDC[i] > 0) &&              // DC status is good for the positive candidate
					(Stat  [i] > 0) &&              // Global status is good for the positive candidate
					(charge[i] < 0)                 // Charge is negative
			       )
			{

			}

		}

		// --------------------------------------------------------------------------------------------------
		// Prep the output tree
		//e_vz     = e_vz_corrected;
		//e_mom[0] = px[0];
		//e_mom[1] = py[0];
		//e_mom[2] = pz[0];


		// Fill the output tree
		outtree->Fill();
	}
	cerr << "Finished with the event loop...\n";

	// --------------------------------------------------------------------------------------------------
	// Editing histograms
	hist_e_Nphe1  -> SetLineColor(2);
	hist_e_Nphe2  -> SetLineColor(8);
	hist_e_EC_in1 -> SetLineColor(2);
	hist_e_EC_in2 -> SetLineColor(8);
	hist_e_EC_out1-> SetLineColor(2);
	hist_e_EC_out2-> SetLineColor(8);
	hist_e_EC_tot1-> SetLineColor(2);
	hist_e_EC_tot2-> SetLineColor(8);
	hist_e_vz     -> SetLineColor(2);
	// ---
	TLegend * leg = new TLegend(0.5,0.5,0.8,0.8);
	leg -> AddEntry(hist_e_Nphe0,"before cuts"       );
	leg -> AddEntry(hist_e_Nphe1,"after PID cuts"    );
	leg -> AddEntry(hist_e_Nphe2,"after PID+Fiducial");
	// --------------------------------------------------------------------------------------------------
	TCanvas *c1 = new TCanvas("c1");
	hist_e_Nphe0 -> Draw();
	hist_e_Nphe1 -> Draw("same");
	hist_e_Nphe2 -> Draw("same");
	leg          -> Draw("same");

	TCanvas *c2 = new TCanvas("c2");
	hist_e_EC_in0 -> Draw();
	hist_e_EC_in1 -> Draw("same");
	hist_e_EC_in2 -> Draw("same");
	leg           -> Draw("same");

	TCanvas *c3 = new TCanvas("c3");
	hist_e_EC_out0 -> Draw();
	hist_e_EC_out1 -> Draw("same");
	hist_e_EC_out2 -> Draw("same");
	leg            -> Draw("same");

	TCanvas *c4 = new TCanvas("c4");
	hist_e_EC_tot0 -> Draw();
	hist_e_EC_tot1 -> Draw("same");
	hist_e_EC_tot2 -> Draw("same");
	leg            -> Draw("same");

	TCanvas *c5 = new TCanvas("c5");
	c5 -> Divide(3,1);
	c5 -> cd(1);	hist_e_thetaMom0 -> Draw("COLZ"); 
	c5 -> cd(2);	hist_e_thetaMom1 -> Draw("COLZ");
	c5 -> cd(3);	hist_e_thetaMom2 -> Draw("COLZ");

	TCanvas *c6 = new TCanvas("c6");
	c6 -> Divide(2,1);
	c6 -> cd(1);	hist_e_momCor -> Draw();
	c6 -> cd(2);	hist_e_momMomCor -> Draw("COLZ");

	TCanvas *c7 = new TCanvas("c7");
	hist_e_vzVzCor -> Draw("COLZ");

	TCanvas *c8 = new TCanvas("c8");
	c8 -> Divide(3,1);
	c8 -> cd(1);	hist_e_Ein_Eout0 -> Draw("COLZ");
	c8 -> cd(2);	hist_e_Ein_Eout1 -> Draw("COLZ");
	c8 -> cd(3);	hist_e_Ein_Eout2 -> Draw("COLZ");

	TCanvas *c9 = new TCanvas("c9");
	c9 -> Divide(3,1);
	c9 -> cd(1);	hist_e_xyEC_hit0 -> Draw("COLZ");
	c9 -> cd(2);	hist_e_xyEC_hit1 -> Draw("COLZ");
	c9 -> cd(3);	hist_e_xyEC_hit2 -> Draw("COLZ");

	TCanvas *c10 = new TCanvas("c10");
	c10 -> Divide(3,1);
	c10 -> cd(1);	hist_e_p_Etot0 -> Draw("COLZ");
	c10 -> cd(2);	hist_e_p_Etot1 -> Draw("COLZ");
	c10 -> cd(3);	hist_e_p_Etot2 -> Draw("COLZ");

	TCanvas *c11 = new TCanvas("c11");
	c11 -> Divide(3,1);
	c11 -> cd(1);	hist_e_p_E0 -> Draw("COLZ");
	c11 -> cd(2);	hist_e_p_E1 -> Draw("COLZ");
	c11 -> cd(3);	hist_e_p_E2 -> Draw("COLZ");

	TCanvas *c12 = new TCanvas("c12");	
	hist_e_phiTheta0 -> Draw("COLZ");
	TCanvas *c13 = new TCanvas("c13");
	hist_e_phiTheta1 -> Draw("COLZ");
	TCanvas *c14 = new TCanvas("c14");
	hist_e_phiTheta2 -> Draw("COLZ");

	TCanvas *c15 = new TCanvas("c15");
	hist_e_vz0 -> Draw();
	hist_e_vz  -> Draw("same");

	TCanvas *c16 = new TCanvas("c16");
	c16 -> Divide(3,2);
	c16 -> cd(1);	hist_e_vz_sec10 -> Draw();
	c16 -> cd(2);	hist_e_vz_sec20 -> Draw();
	c16 -> cd(3);	hist_e_vz_sec30 -> Draw();
	c16 -> cd(4);	hist_e_vz_sec40 -> Draw();
	c16 -> cd(5);	hist_e_vz_sec50 -> Draw();
	c16 -> cd(6);	hist_e_vz_sec60 -> Draw();

	TCanvas *c17 = new TCanvas("c17");
	c17 -> Divide(3,2);
	c17 -> cd(1);	hist_e_vz_sec1 -> Draw();
	c17 -> cd(2);	hist_e_vz_sec2 -> Draw();
	c17 -> cd(3);	hist_e_vz_sec3 -> Draw();
	c17 -> cd(4);	hist_e_vz_sec4 -> Draw();
	c17 -> cd(5);	hist_e_vz_sec5 -> Draw();
	c17 -> cd(6);	hist_e_vz_sec6 -> Draw();

	TCanvas *c18 = new TCanvas("c18");
	c18 -> Divide(2,1);
	c18 -> cd(1);	hist_e_phiVz0 -> Draw("COLZ");
	c18 -> cd(2);	hist_e_phiVz  -> Draw("COLZ");

	TCanvas *c19 = new TCanvas("c19");
	c19 -> Divide(2,1);
	c19 -> cd(1);	hist_e_thetaVz0 -> Draw("COLZ");
	c19 -> cd(2);	hist_e_thetaVz  -> Draw("COLZ");

	TCanvas *c20 = new TCanvas("c20");
	hist_pos_pBeta -> Draw("COLZ");

	TCanvas *c21 = new TCanvas("c21");
	hist_p_pBeta -> Draw("COLZ");

	TCanvas *c22 = new TCanvas("c22");
	c22 -> Divide(2,1);
	c22 -> cd(1);	hist_p_mass ->Draw();
	c22 -> cd(2);	hist_p_pMass->Draw("COLZ");

	TCanvas *c23 = new TCanvas("c23");
	hist_p_p_momCor -> Draw("COLZ");

	TCanvas *c24 = new TCanvas("c24");
	hist_p_vzVzCor -> Draw("COLZ");

	TCanvas *c25 = new TCanvas("c25");
	c25 -> Divide(2,1);
	c25 -> cd(1);	hist_p_phiVz0 -> Draw("COLZ");
	c25 -> cd(2);	hist_p_phiVz  -> Draw("COLZ");

	TCanvas *c26 = new TCanvas("c26");
	c26 -> Divide(2,1);
	c26 -> cd(1);	hist_p_thetaVz0 -> Draw("COLZ");
	c26 -> cd(2);	hist_p_thetaVz  -> Draw("COLZ");

	TCanvas *c27 = new TCanvas("c27");
	hist_p_deltaTmom0 -> Draw("COLZ");

	TCanvas *c28 = new TCanvas("c28");
	hist_p_deltaTmom1 -> Draw("COLZ");

	TCanvas *c29 = new TCanvas("c29");
	hist_p_deltaTmom2 -> Draw("COLZ");

	TCanvas *c30 = new TCanvas("c30");
	hist_p_phiTheta0 -> Draw("COLZ");

	TCanvas *c31 = new TCanvas("c31");
	hist_p_phiTheta1 -> Draw("COLZ");

	TCanvas *c32 = new TCanvas("c32");
	hist_p_phiTheta2 -> Draw("COLZ");

	TCanvas *c33 = new TCanvas("c33");
	hist_pip_deltaTmom0 -> Draw("COLZ");

	TCanvas *c34 = new TCanvas("c34");
	hist_pip_deltaTmom1 -> Draw("COLZ");

	TCanvas *c35 = new TCanvas("c35");
	hist_pip_deltaTmom2 -> Draw("COLZ");

	TCanvas *c36 = new TCanvas("c36");
	hist_neu_pBeta -> Draw("COLZ");

	TCanvas *c37 = new TCanvas("c37");
	hist_n_pBeta -> Draw("COLZ");

	TCanvas *c38 = new TCanvas("c38");
	hist_n_phiTheta0 -> Draw("COLZ");

	TCanvas *c39 = new TCanvas("c39");
	hist_n_phiTheta1 -> Draw("COLZ");

	TCanvas *c40 = new TCanvas("c40");
	hist_n_phiTheta2 -> Draw("COLZ");
	// --------------------------------------------------------------------------------------------------
	// Print histograms on a pdf file
	c1  -> Print("plots.pdf(","pdf");
	c2  -> Print("plots.pdf" ,"pdf");
	c3  -> Print("plots.pdf" ,"pdf");
	c4  -> Print("plots.pdf" ,"pdf");
	c5  -> Print("plots.pdf" ,"pdf");
	c6  -> Print("plots.pdf" ,"pdf");
	c7  -> Print("plots.pdf" ,"pdf");
	c8  -> Print("plots.pdf" ,"pdf");
	c9  -> Print("plots.pdf" ,"pdf");
	c10 -> Print("plots.pdf" ,"pdf");
	c11 -> Print("plots.pdf" ,"pdf");
	c12 -> Print("plots.pdf" ,"pdf");
	c13 -> Print("plots.pdf" ,"pdf");
	c14 -> Print("plots.pdf" ,"pdf");
	c15 -> Print("plots.pdf" ,"pdf");
	c16 -> Print("plots.pdf" ,"pdf");
	c17 -> Print("plots.pdf" ,"pdf");
	c18 -> Print("plots.pdf" ,"pdf");
	c19 -> Print("plots.pdf" ,"pdf");
	c20 -> Print("plots.pdf" ,"pdf");
	c21 -> Print("plots.pdf" ,"pdf");
	c22 -> Print("plots.pdf" ,"pdf");
	c23 -> Print("plots.pdf" ,"pdf");
	c24 -> Print("plots.pdf" ,"pdf");
	c25 -> Print("plots.pdf" ,"pdf");
	c26 -> Print("plots.pdf" ,"pdf");
	c27 -> Print("plots.pdf" ,"pdf");
	c28 -> Print("plots.pdf" ,"pdf");
	c29 -> Print("plots.pdf" ,"pdf");
	c30 -> Print("plots.pdf" ,"pdf");
	c31 -> Print("plots.pdf" ,"pdf");
	c32 -> Print("plots.pdf" ,"pdf");
	c33 -> Print("plots.pdf" ,"pdf");
	c34 -> Print("plots.pdf" ,"pdf");
	c35 -> Print("plots.pdf" ,"pdf");
	c36 -> Print("plots.pdf" ,"pdf");
	c37 -> Print("plots.pdf" ,"pdf");
	c38 -> Print("plots.pdf" ,"pdf");
	c39 -> Print("plots.pdf" ,"pdf");
	c40 -> Print("plots.pdf)","pdf");

	// --------------------------------------------------------------------------------------------------
	// Write the output file
	outfile->cd();
	outtree->Write();
	hist_e_EC_in0       ->Write();
	hist_e_EC_in1       ->Write();
	hist_e_EC_in2       ->Write();
	hist_e_EC_out0      ->Write();
	hist_e_EC_out1      ->Write();
	hist_e_EC_out2      ->Write();
	hist_e_EC_tot0      ->Write();
	hist_e_EC_tot1      ->Write();
	hist_e_EC_tot2      ->Write();
	hist_e_thetaMom0    ->Write();
	hist_e_thetaMom1    ->Write();
	hist_e_thetaMom2    ->Write();
	hist_e_momCor       ->Write();
	hist_e_momMomCor    ->Write();
	hist_e_vzVzCor      ->Write();
	hist_e_Ein_Eout0    ->Write();
	hist_e_Ein_Eout1    ->Write();
	hist_e_Ein_Eout2    ->Write();
	hist_e_xyEC_hit0    ->Write();
	hist_e_xyEC_hit1    ->Write();
	hist_e_xyEC_hit2    ->Write();
	hist_e_p_Etot0      ->Write();
	hist_e_p_Etot1      ->Write();
	hist_e_p_Etot2      ->Write();
	hist_e_p_E0         ->Write();
	hist_e_p_E1         ->Write();
	hist_e_p_E2         ->Write();
	hist_e_Nphe0        ->Write();
	hist_e_Nphe1        ->Write();
	hist_e_Nphe2        ->Write();
	hist_e_phiTheta0    ->Write();
	hist_e_phiTheta1    ->Write();
	hist_e_phiTheta2    ->Write();
	hist_e_vz0          ->Write();
	hist_e_vz           ->Write();
	hist_e_vz_sec10     ->Write();
	hist_e_vz_sec1      ->Write();
	hist_e_vz_sec20     ->Write();
	hist_e_vz_sec2      ->Write();
	hist_e_vz_sec30     ->Write();
	hist_e_vz_sec3      ->Write();
	hist_e_vz_sec40     ->Write();
	hist_e_vz_sec4      ->Write();
	hist_e_vz_sec50     ->Write();
	hist_e_vz_sec5      ->Write();
	hist_e_vz_sec60     ->Write();
	hist_e_vz_sec6      ->Write();
	hist_e_phiVz0       ->Write();
	hist_e_phiVz        ->Write();
	hist_e_thetaVz0     ->Write();
	hist_e_thetaVz      ->Write();
	// ---
	hist_pos_pBeta      ->Write();
	hist_p_pBeta        ->Write();
	hist_pip_pBeta      ->Write();
	// ---
	hist_p_mass         ->Write();
	hist_p_pMass        ->Write();	
	hist_p_p_momCor     ->Write();
	hist_p_vzVzCor      ->Write();
	hist_p_phiVz0       ->Write(); 
	hist_p_phiVz        ->Write();
	hist_p_thetaVz0     ->Write();
	hist_p_thetaVz      ->Write();
	hist_p_deltaTmom0   ->Write();
	hist_p_deltaTmom1   ->Write();
	hist_p_deltaTmom2   ->Write();
	hist_p_phiTheta0    ->Write();
	hist_p_phiTheta1    ->Write();
	hist_p_phiTheta2    ->Write();
	// ---
	hist_pip_deltaTmom0 ->Write();
	hist_pip_deltaTmom1 ->Write();
	hist_pip_deltaTmom2 ->Write();

	// Clean up
	f->Close();
	outfile->Close();

	return 0;
}
