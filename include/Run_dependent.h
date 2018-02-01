#ifndef __RUN_DEPENDENT_H__
#define __RUN_DEPENDENT_H__

#include <string>
#include "TVector3.h"

class TF1;

class Run_dependent
{
	public:
		Run_dependent(int run_number);
		~Run_dependent();

		float ProtonMomCorrection_He3_4Cell(TVector3 V3Pr, float vertex_p );

	private:
		int run;
	
		std::string homedir;
		
		// Helper functions	
		bool read_p_pcor_params();

		// Proton Momentum Correction Data
		double up_parm[6];
		double down_parm[6];


};

#endif
