#include "Fiducial.h"
#include "e2a_constants.h"
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
	el_Ep_ratio_mean = NULL;
	el_Ep_ratio_sig  = NULL;
	prot_deltat_sig  = NULL;
	prot_deltat_mean = NULL;

	// Initialize the key run settings
	E1=E_beam;
	torus_current = torus;
	mini_current = mini;

	homedir = std::string(getenv("HOME"));

	// Read in the various parameters
	bool all_ok=true;
	all_ok &= read_e_fid_params (); // Electron fiducial regions
	all_ok &= read_e_pcor_params(); // Electron momentum corrections
	all_ok &= read_e_pid_params (); // Electron E/p params
	all_ok &= read_p_pid_params (); // Proton delta_t vs mom params
	//	all_ok &= read_vz_cor_params(); // vz corrections
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
bool Fiducial::read_vz_cor_params()
{

	char param_file_name[256];
	sprintf(param_file_name,"%s/.e2a/vz_He_%d.root",homedir.c_str(),E1);

	TFile * old_gfile = gFile;
	TFile * cal_file = new TFile(param_file_name);

	// If we previously set these, we should clean up their memory
	if (vz_corr_func)
		delete vz_corr_func;


	// Pull from file
	vz_corr_func=(TF1*)cal_file->Get("f_vz")->Clone();

	// Put the root global file pointer back to where it was. I hate ROOT. 
	cal_file->Close();
	gFile = old_gfile;

	// Test that the histograms were pulled successfully
	if (!vz_corr_func)
		return false;

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
// ===================================================================================================================================
TVector3 Fiducial::eMomentumCorrection(TVector3 V3el)
{
	// Electron Momentum correction, Pass the electron 4 vector, return corrected 4 Vector pointer.
	// Check out "http://nuclear.unh.edu/~maurik/Personal/E2Root/html/TE2AnaTool.html"

	TVector3      V3ecor(V3el);
	Float_t p   = V3el.Mag();
	Float_t cz  = V3el.CosTheta();
	Float_t phi = 180*V3el.Phi()/pi; 
	if (phi<-30.) phi += 360;
	Float_t theta =  57.29578*(V3el.Theta()); 
	Int_t sectInd = (Int_t)(phi+30)/60;
	if(sectInd>5) sectInd = 5;
	if(sectInd<0) sectInd = 0;
	// -----------------------------------------------------------------------
	// Correction for Ebeam = 2.2GeV and 2250A data
	if ( E1 > 2000 && E1 < 3000 && torus_current > 2240. && torus_current < 2260.){
		phi -= 60.*sectInd;
		p = p*(fgPar_Phi[sectInd][0] + fgPar_Phi[sectInd][1]*phi + fgPar_Phi[sectInd][2]*phi*phi);
		if     (cz>0.800 && cz<0.885){p=p*(fgPar_Theta[0][0]*sin(fgPar_Theta[0][1]*(cz+fgPar_Theta[0][2])) + fgPar_Theta[0][3]);}
		else if(cz>0.885 && cz<0.935){p=p*(fgPar_Theta[1][0] + fgPar_Theta[1][1]*cz + fgPar_Theta[1][2]*cz*cz);}
		else if(cz>0.935 && cz<0.970){p=p*(fgPar_Theta[2][0] + fgPar_Theta[2][1]*cz + fgPar_Theta[2][2]*cz*cz);}
		else {
			std::cerr << "eMomentumCorrection doesn't have correction parameters for the given input. Check it and fix it!\n";
			exit(-3);
		}
	}
	// -----------------------------------------------------------------------
	// Correction for Ebeam = 4.4GeV and 2250A data (corrections valid only for theta > 16)
	else if( E1 > 4000 && E1 < 5000 && torus_current > 2240. && torus_current < 2260.){
		p *= ((fgPar_Phi[sectInd][0] + fgPar_Phi[sectInd][1]*phi                    
					+ fgPar_Phi[sectInd][2]*phi*phi)*(fgPar_Theta[sectInd][3]
						+ (fgPar_Theta[sectInd][2]*theta*theta + fgPar_Theta[sectInd][1]*theta                    
							+ fgPar_Theta[sectInd][0])/exp(theta)));
	}
	// -----------------------------------------------------------------------
	else {                          
		std::cerr << "eMomentumCorrection doesn't have correction parameters for the given input. Check it and fix it!\n";
		exit(-3);
	}
	// -----------------------------------------------------------------------
	if(p!=0.) V3el.SetMag(p);
	return V3ecor;
}
// ===================================================================================================================================
Bool_t Fiducial::pFiducialCut(TVector3 momentum){ //Positive Hadron Fiducial Cut
	//Check out "http://www.jlab.org/Hall-B/secure/e2/bzh/pfiducialcut.html"

	// *********************************************************************************
	//Parameters for 4 GeV proton's Fiducial Cut Rustam Niyazov
	// <A HREF="http://www.physics.odu.edu/~rust/clas/fidp.html"</A> --Rustam Niyazov (ODU).

	const Float_t fgPar_4Gev_2250_Pfidft1l[6][6]={
		{26.2564,0.441269,-29.7632,94.5137,7.71903,2.10915},
		{29.7455,-0.826489,4.09596,91.8187,8.38108,1.5016},
		{29.5399,-0.878321,43.1909,64.9772,11.1844,0.825411},
		{28.5857,0.4061,98.6296,95.5022,13.7297,0.415071},
		{31.9803,0.341766,257.124,103.504,14.2357,0.43387},
		{29.2846,-0.257616,51.1709,84.3207,10.2963,1.69991}};
	const Float_t fgPar_4Gev_2250_Pfidft1r[6][6]={
		{34.7359,-1.45301,660.653,-79.1375,11.3239,1.05352},
		{30.6992,0.71858,442.087,4.20897,3.62722,3.35155},
		{19.1518,3.71404,-197.134,177.828,9.63173,1.35402},
		{23.9897,1.52101,23.9288,71.4476,8.89464,1.69512},
		{22.6619,2.4697,-54.5174,112.22,11.2561,0.687839},
		{20.9859,3.86504,-56.5229,230.635,13.6587,0.270987}};
	const Float_t fgPar_4Gev_2250_Pfidft2l[6][6]={
		{24.683,0.470268,124.501,-9.04329,8.60129,1.66063},
		{26.2736,-0.591497,182.954,-51.059,7.65701,2.29757},
		{24.8681,1.15526,111.322,22.2304,9.46319,1.6834},
		{29.3639,1.307,282.797,89.5863,11.7162,0.376266},
		{36.8099,-0.785452,655.368,46.4935,12.0443,0.500522},
		{25.8401,0.899645,141.723,27.6687,9.62103,1.7379}};
	const Float_t fgPar_4Gev_2250_Pfidft2r[6][6]={
		{32.9905,-0.580968,464.263,30.5379,11.7414,0.320415},
		{26.8867,0.748481,150.349,51.4182,8.70942,1.51013},
		{26.0729,0.357197,136.456,24.1839,6.70568,0.820883},
		{25.8339,1.018,149.648,38.7987,6.56928,0.527773},
		{27.997,0.0685368,268.87,-45.3343,5.26386,3.08026},
		{30.3568,1.60206,359.39,197.047,11.1523,0.451219}};
	const Float_t fgPar_4Gev_2250_Pfidbt1l[6][6]= {
		{-24.4118,4.20154,-0.0480933,-0.0800641,0.000311929,0.000511191},
		{-34.5523,8.81812,0.221281,-0.203846,-0.00115322,0.00119883},
		{-29.4962,6.57417,0.0830637,-0.142094,-0.000271087,0.000801481},
		{-29.5177,6.23458,0.183415,-0.160458,-0.00121912,0.0010282},
		{-19.8091,4.37431,-0.046672,-0.124147,-7.21454e-05,0.000931229},
		{-38.1865,10.6462,0.363126,-0.267793,-0.00212252,0.00162732}};
	const Float_t fgPar_4Gev_2250_Pfidbt1r[6][6]={
		{-15.6987,3.34818,-0.155291,-0.102923,0.000736214,0.000775517},
		{-15.9442,1.75807,-0.196246,-0.0524198,0.00118102,0.000398854},
		{-14.4453,1.65733,-0.269699,-0.0423913,0.00187485,0.000274252},
		{-18.5972,1.41622,-0.144491,-0.0369631,0.000874762,0.000326006},
		{-17.1008,0.577868,-0.173353,-0.021315,0.00108238,0.000189545},
		{2.21904,-3.38706,-0.636698,0.0953525,0.0038789,-0.000559086}};
	const Float_t fgPar_4Gev_2250_Pfidbt2l[6][6]={
		{-13.7253,-1.53789,-0.296133,0.0648705,0.00269427,-0.000928492},
		{-12.356,-2.62192,-0.366191,0.115155,0.0033624,-0.00137599},
		{-2.52638,-9.6591,-0.743505,0.380195,0.0067055,-0.00369404},
		{-34.5804,15.3815,0.417723,-0.489802,-0.00337546,0.00370894},
		{1.87747,-7.70598,-0.919924,0.376373,0.00776553,-0.00354661},
		{-12.3968,-2.37408,-0.367352,0.114661,0.00352523,-0.00148841}};
	const Float_t fgPar_4Gev_2250_Pfidbt2r[6][6]={
		{-29.5895,10.9088,0.248994,-0.326966,-0.00154954,0.00202508},
		{-7.20087,-6.19132,-0.568426,0.257971,0.00476513,-0.00236084},
		{-10.0076,-3.66545,-0.468027,0.163446,0.00421363,-0.00175242},
		{-9.03582,-5.14009,-0.515592,0.221044,0.00482855,-0.00237549},
		{-8.55955,-5.27785,-0.504058,0.201472,0.00404296,-0.00175892},
		{-21.122,5.19264,-0.0761427,-0.0826774,0.0018747,-0.000390706}};
	const Float_t fgPar_4Gev_2250_Pfidbl[6][6]={
		{131.839,-6.64199,-22.8623,4.91185,126.5,20},
		{132.055,-5.2283,2.20945,-1.57951,128.429,11.4286},
		{137.945,-7.90553,-12.8716,3.94534,119.857,22.8571},
		{124.743,-3.54503,-22.8263,5.62231,130.429,11.4286},
		{136.455,-7.59559,-18.6847,4.52149,123.5,20},
		{126.556,-4.02284,-22.2328,5.23298,124.857,22.8571}};
	const Float_t fgPar_4Gev_2250_Pfidbr[6][6]={
		{97.3917,2.99764,26.7715,-5.95695,126.5,20},
		{132.154,-6.60261,0.000146616,1.53542,128.429,11.4286},
		{113.746,-1.24667,32.0728,-9.35241,119.857,22.8571},
		{118.596,-2.44983,22.2973,-5.40976,130.429,11.4286},
		{125.129,-3.96273,21.6178,-5.86908,123.5,20},
		{111.201,-0.178015,25.1267,-6.55928,124.857,22.8571}};
	// *********************************************************************************

	Bool_t status = kTRUE;

	// ----------------------------------------------------------------------------------------------------------------
	if (E1>4000 && E1<5000 && torus_current>2240. && torus_current<2260.){

		Float_t theta = momentum.Theta()*180/pi;
		Float_t phi   = momentum.Phi()  *180/pi;
		if(phi<-30) phi+=360;
		Int_t sector = Int_t ((phi+30)/60);
		if(sector<0) sector=0;
		if(sector>5) sector=5;
		phi -= sector*60;
		Float_t p = momentum.Mag();

		Float_t parfidl [3];    for(Int_t i=0; i<3; i++){parfidl [i]=0;}
		Float_t parfidr [3];    for(Int_t i=0; i<3; i++){parfidr [i]=0;}
		Float_t parfidbl[2];    for(Int_t i=0; i<2; i++){parfidbl[i]=0;}
		Float_t parfidbr[2];    for(Int_t i=0; i<2; i++){parfidbr[i]=0;}
		Float_t cphil =0;       Float_t cphir =0;
		Float_t phi45l=0;       Float_t phi45r=0;
		Float_t phi60l=0;       Float_t phi60r=0;

		Float_t theta_min = 11;
		Float_t theta_max =140;

		bool Forward=kFALSE;            // defines if particle in Forward (Forward=kTRUE) or Backward (Forward=kFALSE) region.
		Int_t   thetab    =45;          // this variable defines the edge point for Forward<->Backward regions
		Float_t p1        =0.575;       // last bin momentum for region p<0.6 GeV/c
		if(p<0.2) p=0.2;                // momentum less than 0.2 GeV/c, use 0.2 GeV/c
		if(p>4.4) p=4.4;                // momentum greater than 4.4 GeV/c, use 4.4 GeV/c
		// ----------------------------------------------------------
		//Get parametrized values of theta_max for p<0.6 GeV/c region
		if(p<0.6){theta_max=fgPar_4Gev_2250_Pfidbl[sector][4]+fgPar_4Gev_2250_Pfidbl[sector][5]*p;}
		//Get parametrized values of theta_max for p>0.6 GeV/c region
		else{theta_max=fgPar_4Gev_2250_Pfidbl[sector][4]+fgPar_4Gev_2250_Pfidbl[sector][5]*p1;}

		//Get the momentum dependent parameters for Forward Region (theta <45 deg)   
		Forward=kTRUE;
		if(p<0.6){//forward1 defines  regions of momenta p<0.6 GeV/c
			//parameters for hyperbolic function
			for (Int_t i=0; i<3; i++){
				Int_t j=2*i;
				parfidl[i]=fgPar_4Gev_2250_Pfidft1l[sector][j]+fgPar_4Gev_2250_Pfidft1l[sector][j+1]/p;
				parfidr[i]=fgPar_4Gev_2250_Pfidft1r[sector][j]+fgPar_4Gev_2250_Pfidft1r[sector][j+1]/p;
			}
		}
		else{//forward2 defines  regions of momenta and p>0.6 GeV/c
			for (Int_t i=0; i<3; i++){
				Int_t j=2*i;
				parfidl[i]=fgPar_4Gev_2250_Pfidft2l[sector][j]+fgPar_4Gev_2250_Pfidft2l[sector][j+1]/p;
				parfidr[i]=fgPar_4Gev_2250_Pfidft2r[sector][j]+fgPar_4Gev_2250_Pfidft2r[sector][j+1]/p;
			}
		}
		// ----------------------------------------------------------
		phi45l=parfidl [0]*(parfidl[2]-45)/(45-parfidl[2]+(parfidl[1]/parfidl[0])); //parametrized value of phi at theta=45 deg.
		phi45r=-parfidr[0]*(parfidr[2]-45)/(45-parfidr[2]+(parfidr[1]/parfidr[0]));
		if(theta>thetab){//backward region defined by theta >45 deg. 
			if(theta>140) theta =140; //theta greater than 140 degrees, use 140 degrees
			if(p>1)p=1.; //momentum greater than 1.0 GeV/c, use 1.0 GeV/c

			//Get the momentum dependent parameters for Backward Region

			Forward=kFALSE;
			if(p<0.6){//backward1 defines  regions of momenta p<0.6 GeV/c
				//parameters for quadratic function
				for (Int_t i=0; i<3; i++){
					Int_t j=2*i;
					parfidl[i]=fgPar_4Gev_2250_Pfidbt1l[sector][j]+fgPar_4Gev_2250_Pfidbt1l[sector][j+1]/p;
					parfidr[i]=fgPar_4Gev_2250_Pfidbt1r[sector][j]+fgPar_4Gev_2250_Pfidbt1r[sector][j+1]/p;
				}
				//these parameters determine theta_flat and phi_edge at p<0.6 GeV/c
				for (Int_t i=0; i<2; i++){
					Int_t j=2*i;
					parfidbl[i]=fgPar_4Gev_2250_Pfidbl[sector][j]+fgPar_4Gev_2250_Pfidbl[sector][j+1]/p;
					parfidbr[i]=fgPar_4Gev_2250_Pfidbr[sector][j]+fgPar_4Gev_2250_Pfidbr[sector][j+1]/p;
				}
			}
			else{//backward2 defines  regions of momenta p>0.6 GeV/c
				//parameters for quadratic function
				for (Int_t i=0; i<3; i++){
					Int_t j=2*i;
					parfidl[i]=fgPar_4Gev_2250_Pfidbt2l[sector][j]+fgPar_4Gev_2250_Pfidbt2l[sector][j+1]/p;
					parfidr[i]=fgPar_4Gev_2250_Pfidbt2r[sector][j]+fgPar_4Gev_2250_Pfidbt2r[sector][j+1]/p;
				}
				//these parameters determine theta_flat and phi_edge at p=0.575 GeV/c momentum
				for (Int_t i=0; i<2; i++){
					Int_t j=2*i;
					parfidbl[i]=fgPar_4Gev_2250_Pfidbl[sector][j]+fgPar_4Gev_2250_Pfidbl[sector][j+1]/p1;
					parfidbr[i]=fgPar_4Gev_2250_Pfidbr[sector][j]+fgPar_4Gev_2250_Pfidbr[sector][j+1]/p1;
				}
			}
		}
		// -------------------------------------------------------------------------------------------------------------------
		if(Forward){//Forward region
			if(p<0.6) theta_min=14; else theta_min=11;//for p<0.6 GeV/c Region theta starts from 14 deg., otherwise 11 deg.   
			cphil=parfidl[0]*(parfidl[2]-theta)/(theta-parfidl[2]+(parfidl[1]/parfidl[0]));//hyperbolic function
			cphir=-parfidr[0]*(parfidr[2]-theta)/(theta-parfidr[2]+(parfidr[1]/parfidr[0]));
		}
		// -------------------------------------------------------------------------------------------------------------------
		else{//Backward region
			phi60l=parfidl[0]+ parfidl[1]*60.+ parfidl[2]*3600.;//parametrized value of phi at theta=60 deg.
			phi60r=-(parfidr[0]+ parfidr[1]*60.+ parfidr[2]*3600.);

			if(theta<60){
				cphil=parfidl[0]+ parfidl[1]*theta+ parfidl[2]*theta*theta; //quadratic function
				cphir=-(parfidr[0]+ parfidr[1]*theta+ parfidr[2]*theta*theta);
			}
			Float_t dl,el,dr,er; //dl and el are theta_flat and phi_edge parameters for phi<0; 
			//dr and er are theta_flat and phi_edge parameters for phi>0;  
			dl=parfidbl[0];el=parfidbl[1];
			dr=parfidbr[0];er=parfidbr[1];

			if(theta>45&&theta<60){ //BackwardA region
				//try to match parametrized values from Forward region to Backward region parameters
				if(cphil>phi45l)cphil=phi45l;
				if(cphir<phi45r)cphir=phi45r;
			}
			//BackwardB region & phi<0
			else if(theta>=60&&theta<=dl){cphil=phi60l;} //phi=constant 
			else if(theta>dl&&theta<=theta_max){
				cphil=(140-theta)*(phi60l-el)/(140-dl) +el;}//phi=stright line 
			else if(theta>theta_max){cphil=0;} //cut out if theta>theta_max
			//BackwardB region & phi>0
			if(theta>=60&&theta<=dr){cphir=phi60r;} //phi=constant 
			else if(theta>dr&&theta<=theta_max){
				cphir=(140-theta)*(phi60r-er)/(140-dr) +er;}//phi=stright line 
			else if(theta>theta_max){cphir=0;} //cut out if theta>theta_max
		}//Backward Region
		// -------------------------------------------------------------------------------------------------------------------
		if(phi<0) status=(phi>cphil); //check the constrains 
		else if(phi>=0) {status=(phi<cphir);
		}

		if(theta<theta_min) status=kFALSE; //Cutting out events below theta_min
		if(Forward && p<0.6 && theta<20.6-11.4*p)status=kFALSE; //function defines cut of the edge at low theta for p<0.6 GeV/c

		//p>0.6 GeV/c. Cut of the edge at low theta  for some sectors and for 
		//some range of momentum, where edge does not look good.
		bool s1s4 =(theta<11.7&&(sector==0||sector==3));
		bool s5   =(theta<12.2&& sector==4);
		bool s6   =(theta<11.4&& sector==5);
		if( p>=0.6 && p<1.5 && (s1s4||s5||s6) ) status=kFALSE;
	}
	// ----------------------------------------------------------------------------------------------------------------
	else{
		std::cerr << "pFiducialCut doesn't have correction parameters for the given input. Check it and fix it!\n";
		exit(-3);
	}
	// ----------------------------------------------------------------------------------------------------------------
	return status;

}

// ===================================================================================================================================
/*
   double vz_corr(double phi,double theta)            //correction function for vertex , takes the arguments in deg.
   {
   return ((-vz_corr_func->GetParameter(1)))*cos((phi-(vz_corr_func->GetParameter(2)))*TMath::DegToRad())/tan(theta*TMath::DegToRad()); 
// vertex correction function obtained for the empty runs 18522, works fine for 3He runs at 4.461[GeV/c] beam energy 
}
 */
// ===================================================================================================================================
TVector3 Fiducial::FindUVW(TVector3 xyz)
{       
	// get the U V W distance to EC edge for the purpose of geometry cut
	// ported from Stepan's function ec_xyz_duvw. the input is lab coordinates of the EC hit.

	Float_t x = xyz.X();
	Float_t y = xyz.Y();
	Float_t z = xyz.Z();
	Float_t xi,yi,zi,u,v,w;
	Float_t ec_the = 0.4363323;
	Float_t ylow   = -182.974 ;
	Float_t yhi    =  189.956 ;
	Float_t tgrho  =  1.95325 ;
	Float_t sinrho = 0.8901256;
	Float_t cosrho = 0.4557150;

	Float_t phi=xyz.Phi()*180./pi;
	if(phi<-30) phi+=360;

	Int_t ec_sect = (phi+30)/60.;
	if(ec_sect<0)ec_sect=0;
	if(ec_sect>5)ec_sect=5;

	Float_t ec_phi = ec_sect*pi/3.;
	xi = -x*sin(ec_phi) + y*cos(ec_phi);
	yi = x*cos(ec_the)*cos(ec_phi) + y*cos(ec_the)*sin(ec_phi) - z*sin(ec_the);
	zi = x*sin(ec_the)*cos(ec_phi) + y*sin(ec_the)*sin(ec_phi) + z*cos(ec_the);
	zi -= 510.32; 

	u =  ( yi-ylow)/sinrho;
	v =  (yhi-ylow)/tgrho - xi + (yhi-yi)/tgrho;
	w = ((yhi-ylow)/tgrho + xi + (yhi-yi)/tgrho)/2./cosrho;

	TVector3 uvw(u,v,w);
	return uvw;
}       
// ===================================================================================================================================
Bool_t Fiducial::CutUVW(TVector3 ecxyz)
{       
	// Cut the edges of EC according to UVW distance threshold defined by par_EcUVW array.
	// If it passes the cut, return true, if not return false

	//parameters for EC edge cuts
	const Float_t par_EcUVW[6][3] = {{60, 360, 400}, {55, 360, 400}, {50, 363, 400}, {52, 365, 396}, {60, 360, 398}, {50, 362, 398}};

	TVector3 ecuvw = FindUVW(ecxyz);
	Float_t phi=ecxyz.Phi()*180/pi;
	if(phi<-30) phi+=360;
	Int_t sector = (phi+30)/60;
	if(sector<0)sector=0;
	if(sector>5) sector=5;
	return (ecuvw.X()>par_EcUVW[sector][0] && ecuvw.Y()<par_EcUVW[sector][1] && ecuvw.Z()<par_EcUVW[sector][2]);
} 
