#ifndef __FIDUCIAL_H__
#define __FIDUCIAL_H__

#include <string>
#include "TVector3.h"

class TF1;

class Fiducial
{
	public:
		Fiducial(int E_beam, int torus_current, int mini_current);
		~Fiducial();

		// Functions to test fiducial and pid bounds
		bool inFidRegion(TVector3 mom, int charge);
		bool in_e_EoverP(double EoverP, double mom, double cut_sigma);
		bool in_p_deltaT(double delta_t, double mom, double cut_sigma);
		TVector3 eMomentumCorrection(TVector3 V3el);
		Bool_t pFiducialCut(TVector3 momentum);
	private:
		int E1;
		int torus_current;
		int mini_current;
		std::string homedir;

		// Helper functions
		bool read_e_fid_params();
		bool read_e_pcor_params();
		bool read_e_pid_params();
		bool read_p_pid_params();

		void getElectronPhiLimits(double mom, double theta, int sector, double &phiMin, double &phiMax);

		// Fiducial Cut Data
		double fgPar_Efid_t0_p[6][2];
		double fgPar_Efid_t1_p[6][6];
		double fgPar_Efid_b_p[6][2][6];
		double fgPar_Efid_a_p[6][2][6];

		// Momentum Correction Data
		double fgPar_Phi[6][3];
		double fgPar_Theta[6][4];

		// Electron PID data
		TF1 *el_Ep_ratio_mean;
		TF1 *el_Ep_ratio_sig;

		// Proton PID data
		TF1 *prot_deltat_sig;
		TF1 *prot_deltat_mean;
};

#endif
