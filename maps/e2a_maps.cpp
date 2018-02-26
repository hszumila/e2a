#include <fstream>
#include <iostream>
#include <cmath>
#include <sstream>

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
	// Setting up output tree and branches
	TTree * outtree = new TTree("T","Skimmed tree");
	double e_vz, e_vz_corrected, e_mom[3], e_phi_mod;
	double p_vz, p_vz_corrected, p_mom_corrected, p_phi_mod;
	double EC_in_cut, el_EC_cut;
	double e_t0,beta_assuming_proton,p_t0,delta_t,beta_assuming_pip,pip_t0,pip_delta_t;
	double corr_px, corr_py, corr_pz, n_px, n_py, n_pz, n_p, EC_Path_corr, Beta_corr;

	TVector3 e_ec_xyz, n_ec_xyz;
	TVector3 T3_e_mom, T3_e_mom_cor, T3_p_mom, u1;

	int nParticles;
	int nProtons, nNeutrons, nPiplus, nPiminus, nPi0;
	int Part_type    [maxPart];
	double Nu, Q2, Xb, Nu_unc, Q2_unc, Xb_unc, t0;
	double vtx_z_unc [maxPart], vtx_z_cor[maxPart];
	double mom_x     [maxPart], mom_y    [maxPart], mom_z  [maxPart]; 
	double e_deltat  [maxPart];
	int    stat_sc   [maxPart], stat_ec  [maxPart], stat_dc[maxPart];
	double sc_time   [maxPart], sc_path  [maxPart];
	double charge    [maxPart], beta     [maxPart];
	/* 
	   =========================
	   Part_type values	(http://www.star.bnl.gov/public/comp/simu/newsite/gstar/kumacs/NewParticle.html)
	   -11  = electron
	   2212 = proton
	   2112 = neutron
	   +211 = pi+
	   -211 = pi-
	   +111 = pi0
	   =========================
	 */ 

	outtree->Branch("nParticles", &nParticles, "nParticles/I"            );
	outtree->Branch("nProtons"  , &nProtons  , "nProtons/I"              );
	outtree->Branch("nNeutrons" , &nNeutrons , "nNeutrons/I"             );
	outtree->Branch("nPiplus"   , &nPiplus   , "nPiplus/I"               );
	outtree->Branch("nPiminus"  , &nPiminus  , "nPiminus/I"              );
	outtree->Branch("t0"        , &t0        , "t0/D"                    );
	outtree->Branch("Nu"        , &Nu        , "Nu/D"                    );
	outtree->Branch("Q2"        , &Q2        , "Q2/D"                    );
	outtree->Branch("Xb"        , &Xb        , "Xb/D"                    );
	outtree->Branch("charge"    ,  charge    , "charge[nParticles]/D"    );
	outtree->Branch("beta"      ,  beta      , "beta[nParticles]/D"      );
	outtree->Branch("Part_type" ,  Part_type , "Part_type[nParticles]/I" );
	outtree->Branch("vtx_z_unc" ,  vtx_z_unc , "vtx_z_unc[nParticles]/D" );
	outtree->Branch("vtx_z_cor" ,  vtx_z_cor , "vtx_z_cor[nParticles]/D" );
	outtree->Branch("mom_x"     ,  mom_x     , "mom_x[nParticles]/D"     );
	outtree->Branch("mom_y"     ,  mom_y     , "mom_y[nParticles]/D"     );
	outtree->Branch("mom_z"     ,  mom_z     , "mom_z[nParticles]/D"     );
	outtree->Branch("e_deltat"  ,  e_deltat  , "e_deltat[nParticles]/D"  );
	outtree->Branch("stat_sc"   ,  stat_sc   , "stat_sc[nParticles]/I"   );
	outtree->Branch("stat_dc"   ,  stat_dc   , "stat_dc[nParticles]/I"   );
	outtree->Branch("stat_ec"   ,  stat_ec   , "stat_ec[nParticles]/I"   );
	outtree->Branch("sc_time"   ,  sc_time   , "sc_time[nParticles]/D"   );
	outtree->Branch("sc_path"   ,  sc_path   , "sc_path[nParticles]/D"   );

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

		// If we get to here, then the electron passed fiducial cuts
		t0 = STT;

		// Without electron momentum correction
		Nu_unc = tab_E1/1000. - T3_e_mom.Mag(); //Energy transfer
		Q2_unc = 4.*tab_E1/1000.*T3_e_mom.Mag()*sin(T3_e_mom.Theta()/2.)*sin(T3_e_mom.Theta()/2.);      //4-momentum transfer^2
		Xb_unc = Q2_unc / (2*mP*Nu_unc);    //Bjorken scaling variable

		// With electron momentum correction (these are saved in the tree)
		Nu = tab_E1/1000. - T3_e_mom_cor.Mag();	//Energy transfer
		Q2 = 4.*tab_E1/1000.*T3_e_mom_cor.Mag()*sin(T3_e_mom_cor.Theta()/2.)*sin(T3_e_mom_cor.Theta()/2.);	//4-momentum transfer^2
		Xb = Q2 / (2*mP*Nu);	//Bjorken scaling variable

		nParticles++;
		Part_type[0] = -11;
		e_deltat [0] =   0;
		mom_x    [0] = T3_e_mom_cor.X();
		mom_y    [0] = T3_e_mom_cor.Y();
		mom_z    [0] = T3_e_mom_cor.Z();
		vtx_z_unc[0] = targetZ[0];
		vtx_z_cor[0] = e_vz_corrected;
		stat_sc  [0] = StatSC [0];
		stat_dc  [0] = StatDC [0];
		stat_ec  [0] = StatEC [0];
		sc_time  [0] = SC_Time[0];
		sc_path  [0] = SC_Path[0];	
		charge   [0] = Charge [0];
		beta     [0] = Beta   [0];

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
	
		// --------------------------------------------------------------------------------------------------
		// Loop over events looking for other particles

		for (int i=1 ; i<gPart ; i++)
		{
			Part_type[i] = 0;

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

						Part_type[nParticles] = 2212;
						e_deltat [nParticles] = delta_t;
						mom_x    [nParticles] = T3_p_mom.X();
						mom_y    [nParticles] = T3_p_mom.Y();
						mom_z    [nParticles] = T3_p_mom.Z();
						vtx_z_unc[nParticles] = targetZ  [i];	
						vtx_z_cor[nParticles] = p_vz_corrected;
						stat_sc  [nParticles] = StatSC [i];
						stat_dc  [nParticles] = StatDC [i];
						stat_ec  [nParticles] = StatEC [i];
						sc_time  [nParticles] = SC_Time[i];
						sc_path  [nParticles] = SC_Path[i];	
						charge   [nParticles] = Charge [i];
						beta     [nParticles] = Beta   [i];

						h2_p_deltaTmom2-> Fill(delta_t   ,mom                   [i]);
						h2_p_phiTheta2 -> Fill(phi    [i],theta                 [i]);	
						h2_p_phiVz0    -> Fill(phi    [i],targetZ               [i]);	
						h2_p_thetaVz0  -> Fill(theta  [i],targetZ               [i]);	
						h2_p_pBeta     -> Fill(mom    [i],Beta                  [i]);

						nProtons++;
						nParticles++;
					}
					// --------------------------------------------------------------------
					// Look specifically for pions +
					else if((id_guess[i] == 211)&&       // Guess at the particle ID is good for the pion+ candidate
							(fid_params.in_pip_deltaT(pip_delta_t, mom[i], pipdeltat_sig_cutrange)) // Pi+ PID
					       ){
						Part_type[nParticles] = 211;
						e_deltat [nParticles] = pip_delta_t;
						mom_x    [nParticles] = T3_p_mom.X();
						mom_y    [nParticles] = T3_p_mom.Y();
						mom_z    [nParticles] = T3_p_mom.Z();
						vtx_z_unc[nParticles] = targetZ  [i];
						vtx_z_cor[nParticles] = p_vz_corrected;
						stat_sc  [nParticles] = StatSC [i];
						stat_dc  [nParticles] = StatDC [i];
						stat_ec  [nParticles] = StatEC [i];
						sc_time  [nParticles] = SC_Time[i];
						sc_path  [nParticles] = SC_Path[i];		
						charge   [nParticles] = Charge [i];
						beta     [nParticles] = Beta   [i];

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

					Part_type[nParticles] = -211;
					e_deltat [nParticles] = pip_delta_t;
					mom_x    [nParticles] = T3_p_mom.X();
					mom_y    [nParticles] = T3_p_mom.Y();
					mom_z    [nParticles] = T3_p_mom.Z();
					vtx_z_unc[nParticles] = targetZ  [i];
					vtx_z_cor[nParticles] = p_vz_corrected;
					stat_sc  [nParticles] = StatSC [i];
					stat_dc  [nParticles] = StatDC [i];
					stat_ec  [nParticles] = StatEC [i];
					sc_time  [nParticles] = SC_Time[i];
					sc_path  [nParticles] = SC_Path[i];	
					charge   [nParticles] = Charge [i];
					beta     [nParticles] = Beta   [i];

					h2_pim_pBeta      -> Fill(mom[i]          ,Beta [i]);
					h2_pim_deltaTmom2 -> Fill(pip_delta_t     ,mom  [i]);

					nPiminus++;
					nParticles++;
				}
			}

		}

		// --------------------------------------------------------------------------------------------------
		// Fill the output tree
		outtree->Fill();

	}
	cerr << "Finished with the event loop...\n";

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
	h1_e_Nphe0 -> Draw();
	h1_e_Nphe1 -> Draw("same");
	h1_e_Nphe2 -> Draw("same");
	leg          -> Draw("same");

	TCanvas *c2 = new TCanvas("c2");
	h1_e_EC_in0 -> Draw();
	h1_e_EC_in1 -> Draw("same");
	h1_e_EC_in2 -> Draw("same");
	leg           -> Draw("same");

	TCanvas *c3 = new TCanvas("c3");
	h1_e_EC_out0 -> Draw();
	h1_e_EC_out1 -> Draw("same");
	h1_e_EC_out2 -> Draw("same");
	leg            -> Draw("same");

	TCanvas *c4 = new TCanvas("c4");
	h1_e_EC_tot0 -> Draw();
	h1_e_EC_tot1 -> Draw("same");
	h1_e_EC_tot2 -> Draw("same");
	leg            -> Draw("same");

	TCanvas *c5 = new TCanvas("c5");
	c5 -> Divide(3,1);
	c5 -> cd(1);	h2_e_thetaMom0 -> Draw("COLZ"); 
	c5 -> cd(2);	h2_e_thetaMom1 -> Draw("COLZ");
	c5 -> cd(3);	h2_e_thetaMom2 -> Draw("COLZ");

	TCanvas *c13 = new TCanvas("c13");
	c13 -> Divide(3,1);
	c13 -> cd(1);	h2_e_Ein_Eout0 -> Draw("COLZ");
	c13 -> cd(2);	h2_e_Ein_Eout1 -> Draw("COLZ");
	c13 -> cd(3);	h2_e_Ein_Eout2 -> Draw("COLZ");

	TCanvas *c14 = new TCanvas("c14");
	c14 -> Divide(3,1);
	c14 -> cd(1);    h2_e_Ein_Eout_0 -> Draw("COLZ");
	c14 -> cd(2);    h2_e_Ein_Eout_1 -> Draw("COLZ");
	c14 -> cd(3);    h2_e_Ein_Eout_2 -> Draw("COLZ");

	TCanvas *c15 = new TCanvas("c15");
	c15 -> Divide(3,1);
	c15 -> cd(1);	h2_e_xyEC_hit0 -> Draw("COLZ");
	c15 -> cd(2);	h2_e_xyEC_hit1 -> Draw("COLZ");
	c15 -> cd(3);	h2_e_xyEC_hit2 -> Draw("COLZ");

	TCanvas *c16 = new TCanvas("c16");
	c16 -> Divide(3,1);
	c16 -> cd(1);	h2_e_p_Etot0 -> Draw("COLZ");
	c16 -> cd(2);	h2_e_p_Etot1 -> Draw("COLZ");
	c16 -> cd(3);	h2_e_p_Etot2 -> Draw("COLZ");

	TCanvas *c17 = new TCanvas("c17");
	c17 -> Divide(3,1);
	c17 -> cd(1);	h2_e_p_E0 -> Draw("COLZ");
	c17 -> cd(2);	h2_e_p_E1 -> Draw("COLZ");
	c17 -> cd(3);	h2_e_p_E2 -> Draw("COLZ");

	TCanvas *c18 = new TCanvas("c18");	
	h2_e_phiTheta0 -> Draw("COLZ");
	TCanvas *c19 = new TCanvas("c19");
	h2_e_phiTheta1 -> Draw("COLZ");
	TCanvas *c20 = new TCanvas("c20");
	h2_e_phiTheta2 -> Draw("COLZ");

	TCanvas *c21 = new TCanvas("c21");
	h1_e_vz0 -> Draw();
	h1_e_vz  -> Draw("same");

	TCanvas *c24 = new TCanvas("c24");
	h2_e_phiVz0 -> Draw("COLZ");
	
	TCanvas *c25 = new TCanvas("c25");
	h2_e_thetaVz0 -> Draw("COLZ");

	TCanvas *c26 = new TCanvas("c26");
	h2_pos_pBeta -> Draw("COLZ");

	TCanvas *c27 = new TCanvas("c27");
	h2_p_pBeta -> Draw("COLZ");

	TCanvas *c28 = new TCanvas("c28");
	c28 -> Divide(2,1);
	c28 -> cd(1);	h1_p_mass ->Draw();
	c28 -> cd(2);	h2_p_pMass->Draw("COLZ");

	TCanvas *c33 = new TCanvas("c33");
	h2_p_phiVz0 -> Draw("COLZ");

	TCanvas *c34 = new TCanvas("c34");
	h2_p_thetaVz0 -> Draw("COLZ");

	TCanvas *c35 = new TCanvas("c35");
	h2_p_deltaTmom0 -> Draw("COLZ");

	TCanvas *c36 = new TCanvas("c36");
	h2_p_deltaTmom1 -> Draw("COLZ");

	TCanvas *c37 = new TCanvas("c37");
	h2_p_deltaTmom2 -> Draw("COLZ");

	TCanvas *c38 = new TCanvas("c38");
	h2_p_phiTheta0 -> Draw("COLZ");

	TCanvas *c39 = new TCanvas("c39");
	h2_p_phiTheta1 -> Draw("COLZ");

	TCanvas *c40 = new TCanvas("c40");
	h2_p_phiTheta2 -> Draw("COLZ");

	TCanvas *c41 = new TCanvas("c41");
	h2_pip_deltaTmom0 -> Draw("COLZ");

	TCanvas *c42 = new TCanvas("c42");
	h2_pip_deltaTmom1 -> Draw("COLZ");

	TCanvas *c43 = new TCanvas("c43");
	h2_pip_deltaTmom2 -> Draw("COLZ");

	TCanvas *c53 = new TCanvas("c53");
	h2_pim_deltaTmom0 -> Draw("COLZ");

	TCanvas *c54 = new TCanvas("c54");
	h2_pim_deltaTmom1 -> Draw("COLZ");

	TCanvas *c55 = new TCanvas("c55");
	h2_pim_deltaTmom2 -> Draw("COLZ");

	// --------------------------------------------------------------------------------------------------
        // Print histograms on a pdf file

        c1  -> Print(Form("./e2a_maps_%d.pdf(",tab_E1) ,"pdf"); 
        c2  -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        c3  -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        c4  -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        c5  -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");   
        c13 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        c14 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        c15 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        c16 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        c17 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        c18 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        c19 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        c20 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        c21 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        c24 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        c25 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        c26 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        c27 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        c28 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        c33 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        c34 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        c35 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        c36 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        c37 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        c38 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        c39 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        c40 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        c41 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        c42 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        c43 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf"); 
        c53 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        c54 -> Print(Form("./e2a_maps_%d.pdf" ,tab_E1) ,"pdf");
        c55 -> Print(Form("./e2a_maps_%d.pdf)",tab_E1) ,"pdf");

        // --------------------------------------------------------------------------------------------------

	// --------------------------------------------------------------------------------------------------
	// Write the output file
	outfile->cd();
	outtree->Write();
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
