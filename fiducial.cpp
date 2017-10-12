#include "fiducial.h"
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>

Fiducial::Fiducial(std::string E_beam, int torus, int mini)
{
  // Initialize the key run settings
  E1=E_beam;
  torus_current = torus;
  mini_current = mini;

  std::string homedir(getenv("HOME"));

  // Read in the Fiducial Cut Parameters
  char param_file_name[256];
  sprintf(param_file_name,"%s/.e2a/FCP_%s_%d.dat",homedir.c_str(),E1.c_str(),torus);
  std::ifstream param_file(param_file_name);
  int param_type, sector;
  double data[6];
  while ( param_file >> param_type )
    {
      param_file >> sector >> data[0] >> data[1] >> data[2] >> data[3] >> data[4] >> data[5];

      // Test the type of parameter and assign it to the proper data array
      switch (param_type)
	{
	case  0:
	  for(int k=0; k<2; k++) fgPar_Efid_t0_p[sector-1][k] = data[k];
	  break;
	case  1:
	  for(int k=0; k<6; k++) fgPar_Efid_t1_p[sector-1][k] = data[k];
	  break;
	case 10:
	  for(int k=0; k<6; k++) fgPar_Efid_b_p[sector-1][0][k] = data[k];
	  break;
	case 11:
	  for(int k=0; k<6; k++) fgPar_Efid_b_p[sector-1][1][k] = data[k];
	  break;
	case 20:
	  for(int k=0; k<6; k++) fgPar_Efid_a_p[sector-1][0][k] = data[k];
	  break;
	case 21:
	  for(int k=0; k<6; k++) fgPar_Efid_a_p[sector-1][1][k] = data[k];
	  break;
	default:
	  printf("Error in Efid parameter file!\nReceived parameter type %d, which is not found.\nAborting!\n\n\n",param_type);
	  exit(-1);
	}
    } // Done reading in Fiducial Region Parameters
  param_file.close();

  // Read in the Momentum Correction Parameters
  sprintf(param_file_name,"%s/.e2a/EMCP_%s_%d.dat",homedir.c_str(),E1.c_str(),torus);
  param_file.open(param_file_name);
  int cj;
  while (param_file >> param_type)
    {
      param_file >> sector >> cj >> data[0];

      // Assign the data to the correct arrays
      switch (param_type)
	{
	case 0:
	  fgPar_Phi[sector-1][cj] = data[0];
	  break;
	case 1:
	  fgPar_Theta[sector-1][cj] = data[0];
	  break;
	default: 
	  printf("Error in EMCP parameter file!\nReceived parameter type %d, which is not found.\nAborting!\n\n\n",param_type);
	  exit(-2);
	}
    } // Done reading in momentum correction parameters
  param_file.close();
  
}

Fiducial::~Fiducial()
{

}

void Fiducial::getElectronPhiLimits(double mom,double theta, int sector, double &phiMin, double &phiMax)
{
  if ((sector < 0) || (sector > 5))
    {
      std::cerr << "Sector " << sector << " passed to getElectronPhiLimits and is out of range. Check it and fix it!\n";
      exit(-3);
    }

  // Sanitize theta
  double theta_deg = theta * 180./M_PI;
  if (theta_deg < 15.)
    {
      std::cerr << "Theta " << theta_deg << " passed to getElectronPhiLimits and is out of range. Check it and fix it!\n";
      exit(-3);
    }
   
  // Sanitize momentum
  if (mom > 3.7) mom=3.7;
  if (mom < 0.9)
    {
      std::cerr << "Momentum " << mom << " passed to getElectronPhiLimits and is out of range. Check it and fix it!\n";
      exit(-3);
    }

  // Assemble the polynomials
  double t0 = fgPar_Efid_t0_p[sector][0]/pow(mom, fgPar_Efid_t0_p[sector][1]);
  double t1 = 0.; 
  double b[2]={0.,0.};
  double a[2]={0.,0.};
  for(int k=0; k<6; k++)
    {
      double mom_to_the_k = pow(mom,k);
      t1 += fgPar_Efid_t1_p[sector][k]*mom_to_the_k;
      for(int l=0; l<2; l++)
	{
	  a[l] += fgPar_Efid_a_p[sector][l][k]*mom_to_the_k;
	  b[l] += fgPar_Efid_b_p[sector][l][k]*mom_to_the_k;
	}
    }

  // Calculate the limits
  phiMin=sector*M_PI/3.; // Default is the center line of each sector
  phiMax=sector*M_PI/3.;
  if(t1 < 45.) t1 = 45.;
  if((t0 < theta_deg) && (theta_deg < t1))
    {
      phiMin -= M_PI/180.*b[0]*(1. - 1/((theta_deg - t0)/(b[0]/a[0]) + 1.));
      phiMax += M_PI/180.*b[1]*(1. - 1/((theta_deg - t0)/(b[1]/a[1]) + 1.));
    }
}

bool Fiducial::inFidRegion(TVector3 mom, int charge)
{
  // Establish the sector;
  double phi = mom.Phi();
  if (phi < -M_PI/6.)
    phi+= 2.*M_PI;
  int sector = (phi+M_PI/6.)/(M_PI/3.);

  // Get the boundaries of phi
  double phiMin, phiMax;
  getElectronPhiLimits(mom.Mag(),mom.Theta(),sector, phiMin,phiMax);

  return ((phi < phiMax) && (phi>phiMin));
}
