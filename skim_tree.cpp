#include <iostream>
#include <cmath>

#include"TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TF1.h"

#include "Fiducial.h"
#include "Run_dependent.h"
#include "constants.h"

using namespace std;

// Range of momentum for which we have good fits for PID
const double epratio_sig_cutrange=3.; // +/- sigma for PID cut
const double pdeltat_sig_cutrange=3.;

// Target position definition
const double min_Z = -15.; // Extremely wide settings for now. 
const double max_Z =  15.; 

// Difference between positive corrected z vertex and electron corrected z vertex
const double pos_z_cut_min = -2;
const double pos_z_cut_max =  2;

int main(int argc, char ** argv)
{
	if (argc != 3)
	{
		cerr << "Wrong number of arguments. Instead try\n"
			<< "\tskim_tree /path/to/input/file /path/to/output/file\n";
		return -1;	
	}

	// Create an instance of the Fiducial Class to store important calibration params
	Fiducial fid_params(4461, 2250, 6000); // Hard-code the beam energy and toroid currents for now.
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
	// Getting run number from filename
	char run_numc[10];
	strncpy(run_numc, argv[1]+22, 5);
	run_numc[5] = '\0';
	int run_numi = atoi(run_numc);
	cout << "Run number: " << run_numi  << endl;	

	// Create an instance of the Run_dependent Class to store important calibration params
	Run_dependent run_dependent_corrections(run_numi); // Hard-code the run-number from now
	// --------------------------------------------------------------------------------------------------
	// Open up the tree, and get the important data
	TTree * t = (TTree*)f->Get("data");
	const int nEvents = t->GetEntries();
	const int maxPart = 50;
	int gPart, CCPart, DCPart, ECPart, SCPart;
	int StatDC[maxPart], StatCC[maxPart], StatEC[maxPart], StatSC[maxPart], id_guess[maxPart];
	float Xb, STT, Q2, W, Nu, Yb;
	float Stat[maxPart], EC_in[maxPart], EC_out[maxPart], EC_tot[maxPart], Nphe[maxPart], SC_Time[maxPart],
	      SC_Path[maxPart], charge[maxPart], beta[maxPart], mass[maxPart], mom[maxPart], px[maxPart], py[maxPart],
	      pz[maxPart], theta[maxPart], phi[maxPart], targetZ[maxPart], theta_pq[maxPart];
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
	// Diagnostic histograms
	TH1D * hist_e_Nphe0    = new TH1D("hist_e_Nphe0"    ,"e- before  cuts;# photo-electrons in CC;Counts"       ,100,   0.,200.);
	TH2D * hist_e_Ein_Eout0= new TH2D("hist_e_Ein_Eout0","e- before  cuts;E_in/p;E_out/p;Counts"                ,100,   0., 0.5,100, 0.,0.5);
	TH2D * hist_e_p_Etot0  = new TH2D("hist_e_p_Etot0"  ,"e- before  cuts;p [GeV];E_tot/p;Counts"               ,100,   0.,  5.,100, 0.,0.7);
	TH2D * hist_e_phiTheta0= new TH2D("hist_e_phiTheta0","e- before  cuts;Phi [deg];Theta [deg];Counts"         ,100,-100.,380.,100,10.,50.);
	// --
	TH1D * hist_e_Nphe     = new TH1D("hist_e_Nphe"     ,"e- passing cuts;# photo-electrons in CC;Counts"       ,100,   0.,200.);
	TH2D * hist_e_Ein_Eout = new TH2D("hist_e_Ein_Eout" ,"e- passing cuts;E_in/p;E_out/p;Counts"                ,100,   0., 0.5,100, 0.,0.5);
	TH2D * hist_e_p_Etot   = new TH2D("hist_e_p_Etot"   ,"e- passing cuts;p [GeV];E_tot/p;Counts"               ,100,   0.,  5.,100, 0.,0.7);
	TH2D * hist_e_phiTheta = new TH2D("hist_e_phiTheta" ,"e- passing cuts;Phi [deg];Theta [deg];Counts"         ,100,-100.,380.,100,10.,50.);
	// ---
	TH1D * hist_e_vz_sec10 = new TH1D("hist_e_vz_sec10" ,"e- passing cuts, before vtx corr, sector 1;electron vz [cm]; Counts"    ,100, -10., 10.);
	TH1D * hist_e_vz_sec1  = new TH1D("hist_e_vz_sec1"  ,"e- passing cuts,  after vtx corr, sector 1;electron vz [cm]; Counts"    ,100, -10., 10.);
	TH1D * hist_e_vz_sec20 = new TH1D("hist_e_vz_sec20" ,"e- passing cuts, before vtx corr, sector 2;electron vz [cm]; Counts"    ,100, -10., 10.);
	TH1D * hist_e_vz_sec2  = new TH1D("hist_e_vz_sec2"  ,"e- passing cuts,  after vtx corr, sector 2;electron vz [cm]; Counts"    ,100, -10., 10.);
	TH1D * hist_e_vz_sec30 = new TH1D("hist_e_vz_sec30" ,"e- passing cuts, before vtx corr, sector 3;electron vz [cm]; Counts"    ,100, -10., 10.);
	TH1D * hist_e_vz_sec3  = new TH1D("hist_e_vz_sec3"  ,"e- passing cuts,  after vtx corr, sector 3;electron vz [cm]; Counts"    ,100, -10., 10.);
	TH1D * hist_e_vz_sec40 = new TH1D("hist_e_vz_sec40" ,"e- passing cuts, before vtx corr, sector 4;electron vz [cm]; Counts"    ,100, -10., 10.);
	TH1D * hist_e_vz_sec4  = new TH1D("hist_e_vz_sec4"  ,"e- passing cuts,  after vtx corr, sector 4;electron vz [cm]; Counts"    ,100, -10., 10.);
	TH1D * hist_e_vz_sec50 = new TH1D("hist_e_vz_sec50" ,"e- passing cuts, before vtx corr, sector 5;electron vz [cm]; Counts"    ,100, -10., 10.);
	TH1D * hist_e_vz_sec5  = new TH1D("hist_e_vz_sec5"  ,"e- passing cuts,  after vtx corr, sector 5;electron vz [cm]; Counts"    ,100, -10., 10.);
	TH1D * hist_e_vz_sec60 = new TH1D("hist_e_vz_sec60" ,"e- passing cuts, before vtx corr, sector 6;electron vz [cm]; Counts"    ,100, -10., 10.);
	TH1D * hist_e_vz_sec6  = new TH1D("hist_e_vz_sec6"  ,"e- passing cuts,  after vtx corr, sector 6;electron vz [cm]; Counts"    ,100, -10., 10.);
	TH2D * hist_e_phiVz0   = new TH2D("hist_e_phiVz0"   ,"e- passing cuts, before vtx corr; phi [deg];vz [cm];Counts"   ,100,-100.,380.,100,-10,10.);
	TH2D * hist_e_phiVz    = new TH2D("hist_e_phiVz"    ,"e- passing cuts,  after vtx corr; phi [deg];vz [cm];Counts"   ,100,-100.,380.,100,-10,10.);
	// ---
	TH2D * hist_e_thetaVz0 = new TH2D("hist_e_thetaVz0" ,"e- passing cuts, before vtx corr; theta [deg];vz [cm];Counts" ,100, -10., 60.,100,-11,11);
	TH2D * hist_e_thetaVz  = new TH2D("hist_e_thetaVz"  ,"e- passing cuts,  after vtx corr; theta [deg];vz [cm];Counts" ,100, -10., 60.,100,-11,11);
	// ---
	TH2D * hist_e_thetaMom = new TH2D("e_thetaMom"      ,"e- passing fid. cuts;Theta [deg];Mom [GeV];Counts"         , 40,  10., 50., 60, 0., 6.);
	TH2D * hist_e_xQ2      = new TH2D("e_xQ2"           ,"e- passing fid. cuts;x;Q2 [GeV^2];Counts"                  , 40,   0.,  2., 40, 0.,10.);
	TH2D * hist_e_momMomCor= new TH2D("hist_e_momMomCor","e- passing fid. cuts;p [GeV];p corrected - p[GeV];Counts"  , 60,   0.,  6., 60,-.1, .1);
	TH2D * hist_e_vzVzCor  = new TH2D("hist_e_vzVzCor"  ,"e- passing fid. cuts;vz [cm];vz corrected - vz [cm];Counts",100, -20., 20.,100,-1., 1.);
	// ---
	TH1D * hist_p_mass     = new TH1D("hist_p_mass"     ,"+  passing fid. cuts;mass [GeV];Counts"                    ,100,   0.,3.5);
	TH2D * hist_p_pMass    = new TH2D("hist_p_pMass"    ,"+  passing fid. cuts;p [GeV];mass [GeV];Counts"            ,100,   0., 4.,100,  0.,3.5);
	TH2D * hist_p_phiTheta = new TH2D("hist_p_phiTheta" ,"p  passing fid. cuts;Phi [deg];Theta [deg];Counts"         ,100,-100.,380.,100,10.,50.);
	TH2D * hist_p_deltaTmom= new TH2D("hist_p_deltaTmom","p  passing fid. cuts;deltaT;p [GeV];Counts"                , 40,   0.,  7., 40, 0., 5.);
	TH2D * hist_p_pBeta    = new TH2D("hist_p_pBeta"    ,"+  passing fid. cuts;p [GeV];#beta;Counts"                 ,100,   0.,  4.,100, 0.,1.3);
	// ---------------------------------------

	TTree * outtree = new TTree("T","Skimmed tree");
	double e_vz, e_vz_corrected, e_mom[3];
	TVector3 T3_e_mom, T3_e_mom_cor, T3_p_mom;
	int nProtons;
	outtree->Branch("e_vz",&e_vz,"e_vz/D");
	outtree->Branch("e_mom",e_mom,"e_mom[3]/D");
	outtree->Branch("nProtons",&nProtons,"nProtons/I");
	// --------------------------------------------------------------------------------------------------
	// Loop over events
	for (int event=0; event < nEvents ; event++)
	{
		if (event % 100000 == 0)
		{
			cerr << "Working on event " << event << " out of " << nEvents << "\n";
		}

		t->GetEvent(event);

		if (gPart <= 0) continue; // Ignore events that have no particle candidates

		// --------------------------------------------------------------------------------------------------
		// Sector index for electrons
		int e_sect = (int)(phi[0]+30)/60;
		if (e_sect>5) e_sect = 5;
		if (e_sect<0) e_sect = 0;
		// --------------------------------------------------------------------------------------------------
		// Here's where we do electron fiducial cuts

		// Define the electron candidate energy in the EC
		double el_cand_EC = TMath::Max(EC_in[0] + EC_out[0], EC_tot[0]);

		// Electron momentum expressed in a TVector3
		T3_e_mom.SetXYZ(px[0],py[0],pz[0]);	

		// Electron vertex (_z) correction
		e_vz_corrected = targetZ[0]+fid_params.vz_corr(180./M_PI*T3_e_mom.Phi(),180./M_PI*T3_e_mom.Theta());

		// ---
		hist_e_Ein_Eout0 -> Fill(EC_in[0]/mom[0],EC_out[0]/mom[0]);
		hist_e_p_Etot0   -> Fill(mom[0],EC_tot[0]/mom[0]);
		hist_e_Nphe0     -> Fill(Nphe[0]);
		hist_e_phiTheta0 -> Fill(phi[0],theta[0]);
		// ---

		// Decide if the electron passes fiducial cuts
		if (!( (StatEC[0] > 0) && 				// EC status is good for the electron candidate
					(StatDC[0] > 0) &&              // DC status is good for the electron candidate
					(StatCC[0] > 0) && 		// CC status is good for the electron candidate
					(StatSC[0] > 0) && 		// SC status is good for the electron candidate
					(charge[0] < 0) && 		// Electron candidate curvature direction is negative
					(EC_in[0] > 0.055) && 		// Electron candidate has enough energy deposit in inner layer of EC
					(el_cand_EC > 0.33) && 		// Enough total energy in the EC
					(fid_params.in_e_EoverP(el_cand_EC/mom[0],mom[0],epratio_sig_cutrange)) &&	// Electron PID (E/p)
					(180./M_PI*T3_e_mom.Theta()>15.) &&						// Theta > 15 deg
					(fid_params.inFidRegion(T3_e_mom,charge[0])) &&					// Electron theta-phi cut
					(e_vz_corrected > min_Z) && 	// Vertex is within the target region
					(e_vz_corrected < max_Z)
		     ))
		{
			continue;
		}

		// If electron passes all cuts, then momentum-correct it (only works for theta > 16 deg):
		if (180./M_PI*T3_e_mom.Theta()>16.) T3_e_mom_cor = fid_params.eMomentumCorrection(T3_e_mom);
		else 				    T3_e_mom_cor = T3_e_mom;

		// If we get to here, then the electron passed fiducial cuts
		// Fill some diagnostic histograms
		hist_e_thetaMom -> Fill(theta[0],mom[0]);
		hist_e_xQ2      -> Fill(Xb,Q2);
		hist_e_momMomCor-> Fill(T3_e_mom.Mag(),T3_e_mom_cor.Mag()-T3_e_mom.Mag());
		hist_e_vzVzCor  -> Fill(targetZ[0],e_vz_corrected-targetZ[0]);
		hist_e_Ein_Eout -> Fill(EC_in[0]/mom[0],EC_out[0]/mom[0]);
		hist_e_p_Etot   -> Fill(mom[0],EC_tot[0]/mom[0]);
		hist_e_Nphe     -> Fill(Nphe[0]);
		hist_e_phiTheta -> Fill(phi[0],theta[0]);
		hist_e_phiVz0   -> Fill(phi[0],targetZ[0]);
		hist_e_phiVz    -> Fill(phi[0],e_vz_corrected);
		hist_e_thetaVz0 -> Fill(theta[0],targetZ[0]);
		hist_e_thetaVz  -> Fill(theta[0],e_vz_corrected);
		if      (e_sect==0) {hist_e_vz_sec10 -> Fill(targetZ[0]);	hist_e_vz_sec1 -> Fill(e_vz_corrected);}
		else if (e_sect==1) {hist_e_vz_sec20 -> Fill(targetZ[0]);	hist_e_vz_sec2 -> Fill(e_vz_corrected);}
		else if (e_sect==2) {hist_e_vz_sec30 -> Fill(targetZ[0]);	hist_e_vz_sec3 -> Fill(e_vz_corrected);}
		else if (e_sect==3) {hist_e_vz_sec40 -> Fill(targetZ[0]);	hist_e_vz_sec4 -> Fill(e_vz_corrected);}
		else if (e_sect==4) {hist_e_vz_sec50 -> Fill(targetZ[0]);	hist_e_vz_sec5 -> Fill(e_vz_corrected);}
		else if (e_sect==5) {hist_e_vz_sec60 -> Fill(targetZ[0]);	hist_e_vz_sec6 -> Fill(e_vz_corrected);}
		else {cout << "Something is wrong with the definition of sectors" << endl;}
		// --------------------------------------------------------------------------------------------------
		// Loop over events looking for protons
		nProtons=0;      
		for (int i=1 ; i<gPart ; i++)
		{
			T3_p_mom.SetXYZ(px[i],py[i],pz[i]);
			double beta_assuming_proton = mom[i]/sqrt(mom[i]*mom[i] + mP*mP);
			double p_t0 = SC_Time[i] - SC_Path[i]/(beta_assuming_proton * c_cm_ns);
			double e_t0 = SC_Time[0] - SC_Path[0]/c_cm_ns;
			double delta_t = p_t0 - e_t0;

			// Need to do proton vertex correction here
			// And then need to cut based on the vtx difference with respect to electron vtx

			// Test if positive particle
			if( (StatSC[i] > 0) && 				// SC status is good for the proton candidate
					(StatDC[i] > 0) &&              // DC status is good for the proton candidate
					(Stat[i] > 0 )  &&		// Global status is good for the proton candidate
					(fid_params.pFiducialCut(T3_p_mom)) // proton theta-phi cut
			  )
			{
				hist_p_mass      -> Fill(mass[i]);
				hist_p_pMass     -> Fill(mom [i],mass[i]);
				hist_p_pBeta     -> Fill(mom[i],beta [i]);

				// Look specifically for protons (not any positive particle)
				if((id_guess[i] == 2212 ) &&       // Guess at the particle ID is good for the proton candidate
						(fid_params.in_p_deltaT(delta_t, mom[i], pdeltat_sig_cutrange)) // Proton PID (delta T vs p)
				  ){
					nProtons++;	
					hist_p_deltaTmom -> Fill(delta_t,mom [i]);
					hist_p_phiTheta  -> Fill(phi[i],theta[i]);

					/*
					   if(ProtonMomCorrection_He3_4Cell(V4_prot_uncorr,prot_vert_corr) != -1){
					   prot_mom_corr=ProtonMomCorrection_He3_4Cell(V4_prot_uncorr,prot_vert_corr);}
					   else
					   {prot_mom_corr=p[ind_p];}
					 */

				}
			}
		}

		// --------------------------------------------------------------------------------------------------
		// Prep the output tree
		e_vz     = e_vz_corrected;
		e_mom[0] = px[0];
		e_mom[1] = py[0];
		e_mom[2] = pz[0];

		// Fill the output tree
		outtree->Fill();
	}
	cerr << "Finished with the event loop...\n";

	// --------------------------------------------------------------------------------------------------
	// Write the output file
	outfile->cd();
	outtree->Write();
	hist_e_thetaMom     ->Write();
	hist_e_xQ2          ->Write();
	hist_e_momMomCor    ->Write();
	hist_e_vzVzCor      ->Write();
	hist_p_deltaTmom    ->Write();
	hist_p_phiTheta     ->Write();
	hist_e_Ein_Eout0    ->Write();
	hist_e_Ein_Eout     ->Write();
	hist_e_p_Etot0      ->Write();
	hist_e_p_Etot       ->Write();
	hist_e_Nphe0        ->Write();
	hist_e_Nphe         ->Write();
	hist_e_phiTheta0    ->Write();
	hist_e_phiTheta     ->Write();
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
	hist_p_pBeta        ->Write();
	hist_p_mass         ->Write();
	hist_p_pMass        ->Write();

	// Clean up
	f->Close();
	outfile->Close();

	return 0;
}
