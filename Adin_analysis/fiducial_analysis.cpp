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
double gaus2(double *x, double *p);
double fitguess(double*x, double*p);

int main(int argc, char** argv)
{
  gSystem->Load("libTree");
  //Take in raw data and output a saved TF1
  if (argc != 3)
    {
      cerr << "Wrong number of arguments. Try instead:\n"
           << "\tfiducial /path/to/input/file /path/to/output/file\n";
      return -1;
    }
  TFile * infile = new TFile(argv[1]);
  cout << "Starting program..." << endl;
  TFile * outfile = new TFile(argv[2], "RECREATE");

  //Create histograms for each subsection of theta and phi based on the central angles defined at the top of the program.
  TH2D * sector[sectors];
  TH2D * sector_p[sectors][p_bins];
  TH1D * sector_p_theta[sectors][p_bins][theta_bins];
  TH1D * sector_p_phi[sectors][p_bins][phi_bins];

  char title[256];
  for(int l = 0;l<sectors;l++)
    {
      sprintf(title,"Sector_%d",l+1);
      sector[l] = (TH2D*)infile->Get(title);
      for(int i=0;i<p_bins;i++)
        {
          sprintf(title,"Sector_%d_p_%f_to_%f",l+1,p_middle[i]-.05,p_middle[i]+.05);
          sector_p[l][i] = (TH2D*)infile->Get(title);
          for(int k = 0;k<theta_bins;k++)
            {
              sprintf(title,"Sector_%d_p_%f_to_%f_theta_%d_to_%d",l+1,p_middle[i]-.05,p_middle[i]+.05,2*k,2*k+2);
              sector_p_theta[l][i][k] = (TH1D*)infile->Get(title);
            }
          for (int j = 0;j<phi_bins;j++)
            {
              sprintf(title,"Sector_%d_p_%f_to_%f_phi_%f_to_%f",l+1,p_middle[i]-.05,p_middle[i]+.05,phi_middle[l][j]-1,phi_middle[l][j]+1);
              sector_p_phi[l][i][j] = (TH1D*)infile->Get(title);
            }
        }
    }
  cout << "Finished naming histograms...\n"
       << "Starting fitting process..." << endl;

  double p_p_mag = 0;
  double p_p_theta = 0;
  double p_p_phi = 0;
  TVector3 p_p;
  //Loop through events. Meant to now fill the histograms.
  TF1 * prot_boxfit;
  int theta_bin[p_bins];
  int theta_oi[p_bins];
  double phi_lower[p_bins];
  double phi_upper[p_bins];

  for (int sec = 0;sec<sectors;sec++)
    {
      cout << "Working on Sector " << sec+1 << endl;
      for (int p = 0; p<p_bins; p++)
        {
          theta_bin[p] = sector_p[sec][p]->GetXaxis()->FindBin(30-p*.42);
          theta_oi[p] = theta_bin[p]/2;
          for (int theta = theta_bins-1;theta>=theta_oi[p];theta--)
            {
              if (sector_p_theta[sec][p][theta]->Integral(0,2*theta_bins)<500)
                continue;
              prot_boxfit = new TF1("prot_box", fitguess, sector_middle[sec]-30,sector_middle[sec]+30,4);
              double max = sector_p_theta[sec][p][theta]->GetMaximum();
              double left_bin = sector_p_theta[sec][p][theta]->FindFirstBinAbove(max/3);
              double left = sector_p_theta[sec][p][theta]->GetBinCenter(left_bin);
              double right_bin = sector_p_theta[sec][p][theta]->FindLastBinAbove(max/3);
              double right = sector_p_theta[sec][p][theta]->GetBinCenter(right_bin);
              prot_boxfit->SetParameters(left,5,right,max/1.5);
              sector_p_theta[sec][p][theta]->Fit("prot_box","q");
              sector_p_theta[sec][p][theta]->Write();
              if (theta==theta_oi[p])
                {
                  phi_lower[p] = sector_p_theta[sec][p][theta_oi[p]]->GetFunction("prot_box")->GetParameter(0);
                  phi_upper[p] = sector_p_theta[sec][p][theta_oi[p]]->GetFunction("prot_box")->GetParameter(2);
                }
              delete sector_p_theta[sec][p][theta];
            }
          int phi_bin_lower = sector_p[sec][p]->GetYaxis()->FindBin(phi_lower[p]);
          int phi_bin_upper = sector_p[sec][p]->GetYaxis()->FindBin(phi_upper[p]);
          cout << p << " " << phi_lower[p] << " " << phi_bin_lower << endl;
          cout << p << " " << phi_upper[p] << " " << phi_bin_upper << endl;

          for (int phi = 0;phi<phi_bins;phi++)
            {
              if (sector_p_phi[sec][p][phi]->Integral(0,phi_bins*2)<500)
                continue;

              TF1 * fit2 = new TF1("fit",gaus2,10,40,4);
              fit2->SetParameters(sector_p_phi[sec][p][phi]->GetMaximum(),sector_p_phi[sec][p][phi]->GetMean(),2,5);
              sector_p_phi[sec][p][phi]->Fit("fit","q");
              sector_p_phi[sec][p][phi]->Write();
              delete sector_p_phi[sec][p][phi];
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

double fitguess(double*x, double*p)
{
  //Lower bound
  double a = p[0];
  double width = p[1];
  double b = a+width;
  //Upper bound
  double d = p[2];
  double c = d -width;
  double height = p[3];

  if ((*x < a) || (*x > d))
    return 0;

  if ((*x > a) && (*x < b))
    return height*(*x-a)/(b-a);

  //if ((*x > b) && (*x < d))
  //return height;
  if ((*x > b) && (*x < c))
    return height;
  
  if ((*x > c) && (*x < d))
    return height*(d-*x)/(d-c);
}
