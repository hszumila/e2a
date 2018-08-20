#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TVector3.h"
#include "TGraph2D.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TLine.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

int which_sector(double phi);
int which_phi(double phi);
int which_theta(double theta);
int which_mom(double mom);

int main(int argc, char** argv)
{
  gSystem->Load("libTree");
  //Take in raw data and output a saved TF1
  if (argc < 3)
    {
      cerr << "Wrong number of arguments. Try instead:\n"
           << "\tfiducial_histmaker /path/to/input/file /path/to/output/file\n";
      return -1;
    }

  int numfiles = argc - 2;
  TFile * infile[numfiles];
  cout << "Starting program..." << endl;
  TChain intree("T");
  char *file;
  for (int i = 0;i<numfiles;i++)
    {
      file = argv[i+1];
      cout << file << endl;
      intree.Add(file); //Electron and proton Generated
    }

  TFile * outputFile = new TFile(argv[numfiles + 1], "RECREATE");

  const unsigned int nevents = intree.GetEntries();

  int maxPart = 50;
	int nRun;
  int nParticles;
	/*int nProtons, nNeutrons, nPiplus, nPiminus, nPi0;
	double Nu, Q2, Xb, Nu_unc, Q2_unc, Xb_unc, t0;
	double vtx_z_unc [maxPart], vtx_z_cor[maxPart], Mass[maxPart];

  double e_deltat  [maxPart];
	int    stat_sc   [maxPart], stat_ec  [maxPart], stat_dc[maxPart];
	double sc_time   [maxPart], sc_path  [maxPart];
	double ec_time   [maxPart], ec_path  [maxPart];
	double ec_in     [maxPart], ec_out   [maxPart], ec_tot [maxPart];
	double ec_x      [maxPart], ec_y     [maxPart], ec_z   [maxPart];
	double ec_u      [maxPart], ec_v     [maxPart], ec_w   [maxPart];
	double charge    [maxPart], beta     [maxPart];*/
	double mom_x     [maxPart], mom_y    [maxPart], mom_z  [maxPart];
  int Part_type    [maxPart];
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

	intree.SetBranchAddress("nRun"      , &nRun      );
  /*
	intree.SetBranchAddress("nProtons"  , &nProtons  );
	intree.SetBranchAddress("nNeutrons" , &nNeutrons );
	intree.SetBranchAddress("nPiplus"   , &nPiplus   );
	intree.SetBranchAddress("nPiminus"  , &nPiminus  );
	intree.SetBranchAddress("t0"        , &t0        );
	intree.SetBranchAddress("Nu"        , &Nu        );
	intree.SetBranchAddress("Q2"        , &Q2        );
	intree.SetBranchAddress("Xb"        , &Xb        );
	intree.SetBranchAddress("charge"    ,  charge    );
	intree.SetBranchAddress("beta"      ,  beta      );
	intree.SetBranchAddress("vtx_z_unc" ,  vtx_z_unc );
	intree.SetBranchAddress("e_deltat"  ,  e_deltat  );
	intree.SetBranchAddress("stat_sc"   ,  stat_sc   );
	intree.SetBranchAddress("stat_dc"   ,  stat_dc   );
	intree.SetBranchAddress("stat_ec"   ,  stat_ec   );
	intree.SetBranchAddress("sc_time"   ,  sc_time   );
	intree.SetBranchAddress("sc_path"   ,  sc_path   );
	intree.SetBranchAddress("ec_time"   ,  ec_time   );
	intree.SetBranchAddress("ec_path"   ,  ec_path   );
	intree.SetBranchAddress("ec_in"     ,  ec_in     );
	intree.SetBranchAddress("ec_out"    ,  ec_out    );
	intree.SetBranchAddress("ec_tot"    ,  ec_tot    );
  intree.SetBranchAddress("ec_x"      ,  ec_x      );
	intree.SetBranchAddress("ec_y"      ,  ec_y      );
	intree.SetBranchAddress("ec_z"      ,  ec_z      );
	intree.SetBranchAddress("ec_u"      ,  ec_u      );
	intree.SetBranchAddress("ec_v"      ,  ec_v      );
  intree.SetBranchAddress("vtx_z_cor" ,  vtx_z_cor );
  intree.SetBranchAddress("ec_w"      ,  ec_w      );
  intree.SetBranchAddress("Mass"      ,  Mass      );*/
	intree.SetBranchAddress("nParticles", &nParticles);
  intree.SetBranchAddress("mom_x"     ,  mom_x     );
	intree.SetBranchAddress("mom_y"     ,  mom_y     );
	intree.SetBranchAddress("mom_z"     ,  mom_z     );
  intree.SetBranchAddress("Part_type" ,  Part_type );

  // --------------------------------------------------------------------------------------------------
  //Create histograms for each subsection of theta and phi based on the central angles defined at the top of the program.
  TH2D * sector[sectors];
  TH2D * sector_p[sectors][p_bins];
  TH1D * sector_p_theta[sectors][p_bins][theta_bins];
  TH1D * sector_p_phi[sectors][p_bins][phi_bins];
  char name[100];
  char title[100];
  TH1D * vertex[numfiles];
  
  for(int l = 0;l<sectors;l++)
    {
      double lower = sector_middle[l]-30.;
      double upper = sector_middle[l]+30.;
      sprintf(title,"Sector_%d",l+1);
      sector[l] = new TH2D(title,title,theta_bins*2,0,theta_bins*2,60,lower,upper);
      for(int i=0;i<p_bins;i++)
        {
          sprintf(title,"Sector_%d_p_%f_to_%f",l+1,p_middle[i]-.05,p_middle[i]+.05);
          sector_p[l][i] = new TH2D(title, title,theta_bins*2,0,theta_bins*2,60,lower,upper);
          for(int k = 0;k<theta_bins;k++)
            {
              sprintf(title,"Sector_%d_p_%f_to_%f_theta_%d_to_%d",l+1,p_middle[i]-.05,p_middle[i]+.05,2*k,2*k+2);
              sector_p_theta[l][i][k] = new TH1D(title, title,phi_bins*2,lower,upper);
            }
          for (int j = 0;j<phi_bins;j++)
            {
              sprintf(title,"Sector_%d_p_%f_to_%f_phi_%f_to_%f",l+1,p_middle[i]-.05,p_middle[i]+.05,phi_middle[l][j]-1,phi_middle[l][j]+1);
              sector_p_phi[l][i][j] = new TH1D(title, title,theta_bins*2,0,theta_bins*2);
            }
        }
    }
  cout << "Finished naming histograms...\n"
       << "Starting loop through events..." << endl;

  double p_p_mag = 0;
  double p_p_theta = 0;
  double p_p_phi = 0;
  TVector3 p_p;
  //Loop through events. Meant to now fill the histograms.
  for (int event = 0; event < nevents; event++)
    {
      if (event%100000==0)
        cout << "Working on event " << event << endl;

      intree.GetEvent(event);

      for(int part = 0; part<nParticles; part++)
        {
          if (Part_type[part] != -11)
            continue;

          TVector3 p_p(mom_x[part], mom_y[part], mom_z[part]);
          p_p_mag = p_p.Mag();
          p_p_theta = p_p.Theta()*(180/3.141592);
          p_p_phi = p_p.Phi()*(180/3.141592);
          if(p_p_phi <-150)
            p_p_phi +=360;

          int prot_sector = which_sector(p_p_phi);
          if (prot_sector < 0)
            continue;

          sector[prot_sector]->Fill(p_p_theta,p_p_phi);

          int prot_mom = which_mom(p_p_mag);
          if (prot_mom<0)
            continue;

          sector_p[prot_sector][prot_mom]->Fill(p_p_theta,p_p_phi);

          int prot_theta = which_theta(p_p_theta);
          if (prot_theta<0)
            continue;

          sector_p_theta[prot_sector][prot_mom][prot_theta]->Fill(p_p_phi);

          int prot_phi = which_phi(p_p_phi);
          if (prot_phi<0)
            continue;

          sector_p_phi[prot_sector][prot_mom][prot_phi]->Fill(p_p_theta);
        }
    }
  for (int sect = 0; sect<sectors;sect++)
    {
      sector[sect]->Write();
      for (int mom = 0; mom<p_bins;mom++)
        {
          sector_p[sect][mom]->Write();
          for (int phi = 0; phi<phi_bins;phi++)
            {
              sector_p_phi[sect][mom][phi]->Write();
            }

          for (int theta = 0; theta<theta_bins;theta++)
            {
              sector_p_theta[sect][mom][theta]->Write();
            }
        }
    }
  return 0;
}

//This function takes in the phi of a specific event. It then loops through the list of central phi and figures out which one it belongs to. Remember that we are looking to place every phi within two degrees of the central one.
int which_phi(double phi)
{
  for (int i=0 ; i<phi_bins ; i++)
    {
      int prot_sect = which_sector(phi);
      if (fabs(phi-phi_middle[prot_sect][i])<1.)
        {
          return i;
        }
    }
  return -1;
}

//Same as which_phi, but for theta!
int which_theta(double theta)
{
  
  for (int i=0 ; i<theta_bins ; i++)
    {
      if (fabs(theta-theta_middle[i])<1)
        {
          return i;
        }
    }
  return -1;
}

int which_sector(double phi)
{
  for (int i=0 ; i<sectors ; i++)
    {
      if (fabs(phi-sector_middle[i])<30.)
        {
          return i;
        }
    }
  return -1;
}

int which_mom(double mom)
{
  for (int i=0;i<p_bins;i++)
    {
      if (fabs(mom-p_middle[i])<.15)
        {
          return i;
        }
    }
  return -1;
}

double gaus2(double *x, double *p)
{
  double A = p[0];
  double m = p[1];
  double s1 = p[2];
  double s2 = p[3];
  if (*x <= m)
    return A*exp(-((*x-m)*(*x-m))/(2*s1*s1));

  if (*x > m)
    return A*exp(-((*x-m)*(*x-m))/(2*s2*s2));
}
