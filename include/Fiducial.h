#ifndef __FIDUCIAL_H__
#define __FIDUCIAL_H__

#include <string>
#include "TVector3.h"

class TF1;

class Fiducial
{
	public:
		Fiducial(int E_beam, int torus_current, int mini_current, std::string target);
		~Fiducial();

		// Functions to test fiducial and pid bounds
		bool inFidRegion(TVector3 mom, int charge);
		bool in_e_EoverP(double EoverP, double mom, double cut_sigma);
		bool in_p_deltaT(double delta_t, double mom, double cut_sigma);
		bool pFiducialCut(TVector3 momentum);
		bool CutUVW(TVector3 ecxyz);
		double vz_corr(TVector3 T3_mom);
		TVector3 eMomentumCorrection(TVector3 V3el);
		TVector3 FindUVW(TVector3 xyz);
	private:
		int E1;
		int torus_current;
		int mini_current;
		std::string homedir;
		std::string tar;

		// Helper functions
		bool read_e_fid_params ();
		bool read_e_pcor_params();
		bool read_e_pid_params ();
		bool read_p_fid_params ();
		bool read_p_pid_params ();
		bool read_vz_cor_params();	

		void getElectronPhiLimits(double mom, double theta, int sector, double &phiMin, double &phiMax);

		// Fiducial Cut Data
		double fgPar_Efid_t0_p[6][2];
		double fgPar_Efid_t1_p[6][6];
		double fgPar_Efid_b_p[6][2][6];
		double fgPar_Efid_a_p[6][2][6];

		// Momentum Correction Data
		double fgPar_Phi[6][3];
		double fgPar_Theta[6][4];

		// Hadron Correction Data
		double fgPar_Pfidft1l[6][6];
		double fgPar_Pfidft1r[6][6];
		double fgPar_Pfidft2l[6][6];
		double fgPar_Pfidft2r[6][6];
		double fgPar_Pfidbt1l[6][6];
		double fgPar_Pfidbt1r[6][6];
		double fgPar_Pfidbt2l[6][6];
		double fgPar_Pfidbt2r[6][6];
		double fgPar_Pfidbl  [6][6];
		double fgPar_Pfidbr  [6][6];

		// Vertex z correction data
		TF1 *vz_corr_func;		

		// Electron PID data
		TF1 *el_Ep_ratio_mean;
		TF1 *el_Ep_ratio_sig;

		// Proton PID data
		TF1 *prot_deltat_sig;
		TF1 *prot_deltat_mean;
};

#endif
