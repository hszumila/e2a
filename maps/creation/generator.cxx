#include "Riostream.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TStyle.h"
#include "TRint.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TSystem.h"

#include <cstdio>
#include <cstdlib>

using namespace std;

int main(int argc, char **argv){

#ifdef WITHRINT
  TRint *myapp = new TRint("RootSession",&argc,argv,NULL,0);
#else
  TApplication *myapp = new TApplication("myapp",0,0);
#endif

  ofstream outfile;
  outfile.open ("./mctk_uniform.txt");

  //Number of Events to Generate
  Int_t num_p_bins = 100;
  Int_t num_cost_bins = 100;
  Int_t num_generated = 1;
  Int_t nEntries = num_p_bins*num_cost_bins*num_generated;
  cout << "Total Number of Entries = " << nEntries << endl;
  
  //Set number of particles to generate (order:e,p,pi+,pi-)
  Int_t top_num = 2;

  //Make sure seeds are different when submitting multiple jobs
  if(argc!=2){
    cout<<"Will Not Sleep Before Running!"<<endl;
    cout<<"Please Input a Single Number to Engage Sleep Function!"<<endl;
  }
  else{
    Int_t time = atoi(argv[1]);
    if(time>100)time-= 100; //Since Runs start at 101
    time*=3;
    cout<<"Sleeping for "<<time<<" seconds"<<endl;
    gSystem->Exec(Form("sleep %d",time));
  }
  Int_t run_number = atoi(argv[1]);
  //Get Random Number Generator
  TRandom3 *gRandom = new TRandom3(0); //use time as a random seed

  cout << "The Random Number seed is "<<gRandom->GetSeed()<<endl;

  //Set Variables
  Float_t cx[50],cy[50],cz[50],mom_tot[50];
  Int_t pid[50];Float_t mass[50];Int_t charge[50];

  //These Don't Change
  Int_t pid_e = 11; Float_t mass_e = 0.0005; Int_t charge_e = -1; //electron
  Int_t pid_p = 2212; Float_t mass_p = 0.9383; Int_t charge_p = 1; //proton
  Int_t pid_pp = 211; Float_t mass_pp = 0.1396; Int_t charge_pp = 1; //piplus
  Int_t pid_pm = -211; Float_t mass_pm = 0.1396; Int_t charge_pm = -1; //piminums
  Float_t x = 0.000, y = 0.000, z = 0.000 ;
  Float_t t_off = 0.000;
  Int_t flag = 0;

  //Set temp variables
  Double_t cost , phi; //i.e. Cos(theta),Phi
  Float_t px, py, pz;

  //Loop over Data
  for(Int_t p_bin=0;p_bin<num_p_bins;p_bin++){
    for(Int_t cost_bin=0;cost_bin<num_cost_bins;cost_bin++){
      for(Int_t gen=0;gen<num_generated;gen++){
	if (cost_bin >= 67)
	  {
	    //Get Info for the electron
	    mom_tot[0] = gRandom->Uniform(.05*p_bin, .05*p_bin + .05);
	    //To generate uniformly, we should do the following: "cos(theta) = 1 - 2*Uniform[0,1]"
	    //This is the same as "cos(theta) = Uniform[-1,1]
	    //For electrons restrict to forward hemisphere (0-71 degrees), since no electron detectors in back
	    cost = gRandom->Uniform(-1+.02*cost_bin, -.98+.02*cost_bin);
	    phi =  gRandom->Uniform(2*TMath::Pi()*(run_number-1)/360., (2*TMath::Pi()*run_number)/360.);

	    px = mom_tot[0] * TMath::Sin( TMath::ACos(cost) ) * TMath::Cos( phi );
	    py = mom_tot[0] * TMath::Sin( TMath::ACos(cost) ) * TMath::Sin( phi );
	    pz = mom_tot[0] * cost; 

	    cx[0] = px/mom_tot[0]; cy[0] = py/mom_tot[0]; cz[0] = pz/mom_tot[0];
	    pid[0] = pid_e; mass[0] = mass_e; charge[0] = charge_e;
	  }
	//Get Info for the proton
        mom_tot[1] = gRandom->Uniform(.05*p_bin, .05*p_bin + .05);
        cost = gRandom->Uniform(-1+.02*cost_bin, -.98+.02*cost_bin);
        phi =  gRandom->Uniform(2*TMath::Pi()*(run_number-1)/360., (2*TMath::Pi()*run_number)/360.);

        px = mom_tot[1] * TMath::Sin( TMath::ACos(cost) ) * TMath::Cos( phi );
        py = mom_tot[1] * TMath::Sin( TMath::ACos(cost) ) * TMath::Sin( phi );
        pz = mom_tot[1] * cost; 

        cx[1] = px/mom_tot[1]; cy[1] = py/mom_tot[1]; cz[1] = pz/mom_tot[1];
        pid[1] = pid_p; mass[1] = mass_p; charge[1] = charge_p;

        //Get Info for the pi+
        mom_tot[2] = gRandom->Uniform(.05*p_bin, .05*p_bin + .05);
        cost = gRandom->Uniform(-1+.02*cost_bin, -.98+.02*cost_bin);
        phi =  gRandom->Uniform(2*TMath::Pi()*(run_number-1)/360., (2*TMath::Pi()*run_number)/360.);

        px = mom_tot[2] * TMath::Sin( TMath::ACos(cost) ) * TMath::Cos( phi );
        py = mom_tot[2] * TMath::Sin( TMath::ACos(cost) ) * TMath::Sin( phi );
        pz = mom_tot[2] * cost;

        cx[2] = px/mom_tot[2]; cy[2] = py/mom_tot[2]; cz[2] = pz/mom_tot[2];
        pid[2] = pid_pp; mass[2] = mass_pp; charge[2] = charge_pp;

        //Get Info for the pi-
        mom_tot[3] = gRandom->Uniform(.05*p_bin, .05*p_bin + .05);
        cost = gRandom->Uniform(-1+.02*cost_bin, -.98+.02*cost_bin);
        phi =  gRandom->Uniform(2*TMath::Pi()*(run_number-1)/360., (2*TMath::Pi()*run_number)/360.);

        px = mom_tot[3] * TMath::Sin( TMath::ACos(cost) ) * TMath::Cos( phi );
        py = mom_tot[3] * TMath::Sin( TMath::ACos(cost) ) * TMath::Sin( phi );
        pz = mom_tot[3] * cost;

        cx[3] = px/mom_tot[3]; cy[3] = py/mom_tot[3]; cz[3] = pz/mom_tot[3];
        pid[3] = pid_pm; mass[3] = mass_pm; charge[3] = charge_pm;
  
        //------------------------------------------------------------------------
	if (cost_bin >= 67)
	  {
	    for(Int_t k=0;k<top_num;k++){      
	      if(k==0) {outfile << top_num << endl;}
	      outfile << pid[k] <<" "<< cx[k] <<" "<< cy[k] <<" "<< cz[k] <<" "<< mom_tot[k] <<endl;
	      outfile << mass[k] <<" "<< charge[k] << endl;
	      outfile << x <<" "<< y <<" "<< z <<" "<< t_off <<" "<< flag <<endl;
	    }
	  }
	else
	  {
	    for(Int_t k=1;k<top_num;k++){      
	      if(k==1) {outfile << top_num-1 << endl;}
	      outfile << pid[k] <<" "<< cx[k] <<" "<< cy[k] <<" "<< cz[k] <<" "<< mom_tot[k] <<endl;
	      outfile << mass[k] <<" "<< charge[k] << endl;
	      outfile << x <<" "<< y <<" "<< z <<" "<< t_off <<" "<< flag <<endl;  
	    }
	  }
      }
    }
  }  
  //------------------------------------------------------------------------
    
  
  
  outfile.close();
  
  myapp->Run();
  return 0;
}
  