#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TTree.h"

using namespace std;

int main(int argc, char ** argv)
{
  if (argc < 2)
    {
      cerr << "Not enough arguments. Give me a root file I can open.\n";
      return -1;	
    }

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
  const int maxPart = 30;
  int gPart, CCPart, DCPart, ECPart, SCPart;
  int StatCC[maxPart], StatDC[maxPart], StatEC[maxPart], StatSC[maxPart], id_guess[maxPart];
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
  t->SetBranchAddress("Nphe",Nphe); // Number of photons per particle???
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
  t->SetBranchAddress("STT",&STT); // STT
  t->SetBranchAddress("Q2",&Q2); // Momentum transfer
  t->SetBranchAddress("W",&W); // Hadronic mass
  t->SetBranchAddress("Nu",&Nu); // Energy transfer
  t->SetBranchAddress("Yb",&Yb); // Y-scaling variable ???
  
  // Loop over events
  for (int event=0; event < nEvents ; event++)
    {
      t->GetEvent(event);
      
      cout << "Reading event " << event << " which has " << gPart << " particles\n";
      cout << "\tElectron candidate: " << mom[0] << " GeV   " << theta[0] << " deg.   Q^2 = " << Q2 << " GeV^2 \n";      
      // Loop over the particles in the event
      for (int particle=0 ; particle < gPart ; particle++)
	{
	  cout << "\t\tParticle " << particle << " had type " << id_guess[particle] 
	       << "   " << mom[particle] << " GeV    " << theta[particle] << " deg\n";
	}
    }

  f->Close();
  return 0;
}
