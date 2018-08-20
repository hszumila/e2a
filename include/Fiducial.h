#ifndef __FIDUCIAL_H__
#define __FIDUCIAL_H__

#include <string>
#include "TVector3.h"

class TF1;

class Fiducial
{
	public:
		Fiducial(int E_beam, int torus_current, int mini_current, std::string target, bool data);
		~Fiducial();

		// Functions to test fiducial and pid bounds
		bool e_inFidRegion(TVector3 mom);
		bool in_e_EoverP(double EoverP, double mom, double cut_sigma);
		bool in_p_deltaT(double delta_t, double mom, double cut_sigma);
		bool in_pip_deltaT(double delta_t, double mom, double cut_sigma);
		bool in_pim_deltaT(double delta_t, double mom, double cut_sigma);
		bool pFiducialCut(TVector3 momentum);
		bool CutUVW_e(TVector3 ecxyz);
		bool CutUVW  (TVector3 ecuvw, double dist);
		double vz_corr(TVector3 T3_mom);
		double corrected_path_length( double uncorrected_path_length , double E_in , double E_out  );
		TVector3 eMomentumCorrection(TVector3 V3el);
		TVector3 FindUVW(TVector3 xyz);
		double EC_in_cut();
		double el_EC_cut();		  

	private:
		int E1;
		int torus_current;
		int mini_current;
		std::string homedir;
		std::string tar;
		bool is_data;

		// Helper functions
		bool read_e_fid_params     ();
		bool read_e_pcor_params    ();
		bool read_e_pid_params     ();
		bool read_p_fid_params     ();
		bool read_p_pid_params     ();
		bool read_pip_pid_params   ();
		bool read_pim_pid_params   ();
		bool read_vz_cor_params    ();
		bool read_n_pathlength_corr();

		// Fiducial Cut Data
		double fgPar_Efid_t0_p [6][2];
		double fgPar_Efid_t1_p [6][6];
		double fgPar_Efid_b_p  [6][2][6];
		double fgPar_Efid_a_p  [6][2][6];
		double fgPar_Efid      [6][6][9];
		double fgPar_Efid_Theta_S3[4][8];
		double fgPar_Efid_Theta_S4[2][8];
		double fgPar_Efid_Theta_S5[8][8];
		// Momentum Correction Data
		double fgPar_Phi  [6][3];
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
		double fgPar_Pfid_For[6][4][7];
		double fgPar_Pfid_Bak[6][4][7];
		double fgPar_Pfid_ScpdS2[2][6];
		double fgPar_Pfid_ScpdS3[8][6];
		double fgPar_Pfid_ScpdS4[4][6];
		double fgPar_Pfid_ScpdS5[8][6];
		// Vertex z correction data
		TF1 *vz_corr_func;		

		// Electron PID data
		TF1 *el_Ep_ratio_mean;
		TF1 *el_Ep_ratio_sig;

		// Proton PID data
		TF1 *prot_deltat_sig;
		TF1 *prot_deltat_mean;
		
		// Pi+ PID data
		TF1 *pip_deltat_sig;
                TF1 *pip_deltat_mean;

		// Pi- PID data
                TF1 *pim_deltat_sig;
                TF1 *pim_deltat_mean;

		// Neutron pathlength correction
		double pl_corr_in  ;
		double pl_corr_out ;
		double pl_corr_both;
};

#endif
