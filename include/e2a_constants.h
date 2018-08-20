#ifndef __E2A_CONSTANTS_H__
#define __E2A_CONSTANTS_H__

const double c_m_s   = 2.99792E+10;
const double me = 0.000510998;
const double mP = 0.938272;
const double mN = 0.93957;
const double mpc= 0.139570;
const double mp0= 0.134977;
const double pi = 3.14159;
const double ns_to_s = 1.0E-09;

const double sector_middle[6] = {-120,-60,0,60,120,180};
//const double p_middle[100];
//const double phi_middle[6][100];
//const double theta_middle[100];

const int p_bins = 13;
const int sectors = 6;
const int theta_bins = 30;
const int phi_bins = 30;

/*for (int i = 0; i<p_bins ;i++)
{
	p_middle[i] = 1.05+.2*i
}
for (int i = 0;i<theta_bins;i++)
{
	theta_middle[i] = 1.0+2.0*i;
}
for (int i = 0;i<sectors;i++)
{
	for (int j =-15;j<15+phi_bins;j++)
	{
	phi_middle[i][j+15] = sector_middle[i]+j*2;
	}
}*/

#endif
