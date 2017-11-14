#include "Fiducial.h"
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>

#include "TF1.h"
#include "TFile.h"
// ===================================================================================================================================
Fiducial::Fiducial(int E_beam, int torus, int mini)
{
	// Initialize some variables
	el_Ep_ratio_mean=NULL;
	el_Ep_ratio_sig=NULL;
	prot_deltat_sig=NULL;
	prot_deltat_mean=NULL;

	// Initialize the key run settings
	E1=E_beam;
	torus_current = torus;
	mini_current = mini;

	homedir = std::string(getenv("HOME"));

	// Read in the various parameters
	bool all_ok=true;
	all_ok &= read_e_fid_params(); // Electron fiducial regions
	all_ok &= read_e_pcor_params(); // Electron momentum corrections
	all_ok &= read_e_pid_params(); // Electron E/p params
	all_ok &= read_p_pid_params(); // Proton delta_t vs mom params
	if (all_ok)
		std::cerr << "Successfully read in the various parameters...\n";
	else
	{
		std::cerr << "Failed to read in the parameter files. Exiting...\n";
		exit(-1);
	}  
}
// ===================================================================================================================================
Fiducial::~Fiducial()
{
	// Memory clean up
	if (prot_deltat_sig)
		delete prot_deltat_sig;
	if (prot_deltat_mean)
		delete prot_deltat_mean;
	if (el_Ep_ratio_sig)
		delete el_Ep_ratio_sig;
	if (el_Ep_ratio_mean)
		delete el_Ep_ratio_mean;
}
// ===================================================================================================================================
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
// ===================================================================================================================================
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
// ===================================================================================================================================
bool Fiducial::read_e_fid_params()
{
	char param_file_name[256];
	sprintf(param_file_name,"%s/.e2a/FCP_%d_%d.dat",homedir.c_str(),E1,torus_current);
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

	return true;
}
// ===================================================================================================================================
bool Fiducial::read_e_pcor_params()
{
	char param_file_name[256];
	sprintf(param_file_name,"%s/.e2a/EMCP_%d_%d.dat",homedir.c_str(),E1,torus_current);
	std::ifstream param_file(param_file_name);
	param_file.open(param_file_name);
	int param_type, sector;
	double data[6];
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

	return true;
}
// ===================================================================================================================================
bool Fiducial::read_e_pid_params()
{
	char param_file_name[256];
	sprintf(param_file_name,"%s/.e2a/el_Epratio_mom_%d_%d.root",homedir.c_str(),E1,torus_current);
	TFile * old_gfile = gFile;
	TFile * cal_file = new TFile(param_file_name);

	// If we previously set these, we should clean up their memory
	if (el_Ep_ratio_mean)
		delete el_Ep_ratio_mean;
	if (el_Ep_ratio_sig)
		delete el_Ep_ratio_sig;

	// Pull from file
	el_Ep_ratio_mean=(TF1*)cal_file->Get("f_mean")->Clone();
	el_Ep_ratio_sig=(TF1*)cal_file->Get("f_sig")->Clone();

	// Put the root global file pointer back to where it was. I hate ROOT. 
	cal_file->Close();
	gFile = old_gfile;

	// Test that the histograms were pulled successfully
	if (!el_Ep_ratio_mean)
		return false;
	if (!el_Ep_ratio_sig)
		return false;

	return true;
}
// ===================================================================================================================================
bool Fiducial::read_p_pid_params()
{
	char param_file_name[256];
	sprintf(param_file_name,"%s/.e2a/protdeltat_mom_%d_%d.root",homedir.c_str(),E1,torus_current);
	TFile * old_gfile = gFile;
	TFile * file_in1 = new TFile(param_file_name);

	// If we previously set these, we should clean up their memory
	if (prot_deltat_sig)
		delete prot_deltat_sig;
	if (prot_deltat_mean)
		delete prot_deltat_mean;

	// Pull from file
	prot_deltat_sig=(TF1*)file_in1->Get("sig_pol9")->Clone();
	prot_deltat_mean=(TF1*)file_in1->Get("mean_pol9")->Clone();

	// Put the root global file pointer back to where it was. I hate ROOT. 
	file_in1->Close();
	gFile = old_gfile;

	// Test that the histograms were pulled successfully
	if (!prot_deltat_sig)
		return false;
	if (!prot_deltat_mean)
		return false;

	return true;
}
// ===================================================================================================================================
bool Fiducial::in_p_deltaT(double delta_t, double mom, double cut_sigma)
{
	const double prot_mom_lim=2.15;
	if (mom > prot_mom_lim) mom = prot_mom_lim;

	double delta_t_up_limit = prot_deltat_mean->Eval(mom) + cut_sigma * prot_deltat_sig->Eval(mom);
	double delta_t_lo_limit = prot_deltat_mean->Eval(mom) - cut_sigma * prot_deltat_sig->Eval(mom);

	if ((delta_t > delta_t_lo_limit) && (delta_t < delta_t_up_limit))
		return true;
	else
		return false;
}
// ===================================================================================================================================
// ALREADY IMPLEMENTED FOR ELECTRONS IN SKIM_TREE
bool Fiducial::in_e_EoverP(double EoverP, double mom, double cut_sigma)
{
	const double min_el_mom = 1.5;
	const double max_el_mom = 3.8;

	if (mom < min_el_mom)
		return false;
	if (mom > max_el_mom)
		mom = max_el_mom;

	double min_EoverP = el_Ep_ratio_mean->Eval(mom) - cut_sigma * el_Ep_ratio_sig->Eval(mom);
	double max_EoverP = el_Ep_ratio_mean->Eval(mom) + cut_sigma * el_Ep_ratio_sig->Eval(mom);

	if ((EoverP > min_EoverP) && (EoverP < max_EoverP))
		return true;
	else
		return false;
}
