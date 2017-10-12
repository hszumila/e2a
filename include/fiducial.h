#ifndef __FIDUCIAL_H__
#define __FIDUCIAL_H__

#include <string>
#include "TVector3.h"

class Fiducial
{
 public:
  Fiducial(std::string E1, int torus_current, int mini_current);
  ~Fiducial();
  bool inFidRegion(TVector3 mom, int charge);

 private:
  std::string E1;
  int torus_current;
  int mini_current;

  // Fiducial Cut Data
  double fgPar_Efid_t0_p[6][2];
  double fgPar_Efid_t1_p[6][6];
  double fgPar_Efid_b_p[6][2][6];
  double fgPar_Efid_a_p[6][2][6];

  // Momentum Correction Data
  double fgPar_Phi[6][3];
  double fgPar_Theta[6][4];

  void getElectronPhiLimits(double mom, double theta, int sector, double &phiMin, double &phiMax);

};


#endif
