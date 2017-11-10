#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TF1.h"

#include "Fiducial.h"
#include "constants.h"

using namespace std;

// Range of momentum for which we have good fits for PID
const double epratio_sig_cutrange=3.; // +/- sigma for PID cut
const double pdeltat_sig_cutrange=3.;

// Target position definition
const double min_Z = -15.; // Extremely wide settings for now. 
const double max_Z = 15.; 

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

  // Open up the Root file
  TFile * f = new TFile(argv[1]);
  if (f)
    cerr << "Successfully opened file " << argv[1] << "\n";
  else
    {
      cerr << "Could not open file " << argv[1] << "\n\tExiting...\n";
      return -2;
    }

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
  t->SetBranchAddress("gPart",&gPart); // Number of particles observed (globally) in the event
  t->SetBranchAddress("CCPart",&CCPart); // Number of particles observed in the Cherenkovs
  t->SetBranchAddress("DCPart",&DCPart); // Number of particles observed in the Drift Chambers
  t->SetBranchAddress("ECPart",&ECPart); // Number of particles observed in the ECal
  t->SetBranchAddress("SCPart",&SCPart); // Number of particles observed in the Scintillation Counters (ToFs)
  t->SetBranchAddress("Stat",Stat); // Global status for each particle candidate
  t->SetBranchAddress("StatDC",StatDC); // Drift Chamber status for each particle candidate
  t->SetBranchAddress("StatCC",StatCC); // Cherenkov status for each particle
  t->SetBranchAddress("StatEC",StatEC); // ECal status for each particle
  t->SetBranchAddress("StatSC",StatSC); // Scintillation counter status for each particle  
  t->SetBranchAddress("particle",id_guess); // Guess at the particle ID made by the recon software (maybe not reliable)
  t->SetBranchAddress("EC_in",EC_in); // Inner layer of ECal for each particle
  t->SetBranchAddress("EC_out",EC_out); // Outer layer of ECal for each particle
  t->SetBranchAddress("EC_tot",EC_tot); // Total energy deposit in the ECal for each particle
  t->SetBranchAddress("Nphe",Nphe); // Number of photo-electrons per hit in the Cherenkov detectors
  t->SetBranchAddress("SC_Time",SC_Time); // Time in the scintillators per particle
  t->SetBranchAddress("SC_Path",SC_Path); // Path Length per particle
  t->SetBranchAddress("Charge",charge); // Charge per particle
  t->SetBranchAddress("Beta",beta); // Beta per particle
  t->SetBranchAddress("Mass",mass); // Mass per particle
  t->SetBranchAddress("Momentum",mom); // Momentum magnitude per particle
  t->SetBranchAddress("Momentumx",px); // Momentum x component per particle
  t->SetBranchAddress("Momentumx",py); // Momentum y component per particle
  t->SetBranchAddress("Momentumx",pz); // Momentum z component per particle
  t->SetBranchAddress("Theta",theta); // Theta per particle
  t->SetBranchAddress("Phi",phi); // Phi per particle
  t->SetBranchAddress("TargetZ",targetZ); // Target Z per particle
  t->SetBranchAddress("Thetapq",theta_pq); // Angle wrt to q vector per particle
  t->SetBranchAddress("BjorkenX",&Xb); // Bjorken X
  t->SetBranchAddress("STT",&STT); // RF-corrected start time.
  t->SetBranchAddress("Q2",&Q2); // Momentum transfer
  t->SetBranchAddress("W",&W); // Hadronic mass
  t->SetBranchAddress("Nu",&Nu); // Energy transfer
  t->SetBranchAddress("Yb",&Yb); // Y-scaling variable

  // Open up the output file
  TFile * outfile = new TFile(argv[2],"RECREATE");
  TH2D * hist_e_thetaMom = new TH2D("e_thetaMom","e- passing fid. cuts;Theta [deg];Mom [GeV];Counts",40,10.,50.,60,0.,6.);
  TH2D * hist_e_xQ2 = new TH2D("e_xQ2","e- passing fid. cuts;x;Q2 [GeV^2];Counts",40,0.,2.,40,0.,10.);
  TTree * outtree = new TTree("T","Skimmed tree");
  double e_vz, e_mom[3];
  int nProtons;
  outtree->Branch("e_vz",&e_vz,"e_vz/D");
  outtree->Branch("e_mom",e_mom,"e_mom[3]/D");
  outtree->Branch("nProtons",&nProtons,"nProtons/I");

  // Loop over events
  for (int event=0; event < nEvents ; event++)
    {
      if (event % 100000 == 0)
	{
	  cerr << "Working on event " << event << " out of " << nEvents << "\n";
	}

      t->GetEvent(event);

      if (gPart <= 0) continue; // Ignore events that have no particle candidates

      // ************************************************************************
      // Here's where we do electron fiducial cuts

      // Define the electron candidate energy in the EC
      double el_cand_EC = TMath::Max(EC_in[0] + EC_out[0], EC_tot[0]);

      // Decide if the electron passes fiducial cuts
      if (!( (StatEC[0] > 0) && // EC status is good for the electron candidate
	     (StatCC[0] > 0) && // CC status is good for the electron candidate
	     (StatSC[0] > 0) && // SC status is good for the electron candidate
	     (charge[0] < 0) && // Electron candidate curvature direction is negative
	     (EC_in[0] > 0.055) && // Electron candidate has enough energy deposit in inner layer of EC
	     (el_cand_EC > 0.33) && // Enough total energy in the EC
	     (fid_params.in_e_EoverP(el_cand_EC/mom[0],mom[0],epratio_sig_cutrange)) &&	     
	     (targetZ[0] > min_Z) && // Vertex is within the target region
	     (targetZ[0] < max_Z)
	     ))
	{
	  continue;
	}

      // If we get to here, then the electron passed fiducial cuts
      // Fill some diagnostic histograms
      hist_e_thetaMom->Fill(theta[0],mom[0]);
      hist_e_xQ2->Fill(Xb,Q2);

      // Loop over events looking for protons
      nProtons=0;      
      for (int i=1 ; i<gPart ; i++)
	{
	  double beta_assuming_proton = mom[i]/sqrt(mom[i]*mom[i] + mP*mP);
	  double p_t0 = SC_Time[i] - SC_Path[i]/(beta_assuming_proton * c_cm_ns);
	  double e_t0 = SC_Time[0] - SC_Path[0]/c_cm_ns;
	  double delta_t = p_t0 - e_t0;

	  // Test if proton
	  if( (StatSC[i] > 0) && 
	      (Stat[i] > 0 )  &&
	      (id_guess[i] == 2212 ) &&
	      (fid_params.in_p_deltaT(delta_t, mom[i], pdeltat_sig_cutrange))
	      )
	    {
	      // Then we have a proton
	      nProtons++;
	    }
	}

      // Prep the output tree
      e_vz=targetZ[0];
      e_mom[0] = px[0];
      e_mom[1] = py[0];
      e_mom[2] = pz[0];

      // Fill the output tree
      outtree->Fill();
    }
  cerr << "Finished with the event loop...\n";

  // Write the output file
  outfile->cd();
  outtree->Write();
  hist_e_thetaMom->Write();
  hist_e_xQ2->Write();

  // Clean up
  f->Close();
  outfile->Close();

  return 0;
}
