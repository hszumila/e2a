#include <iostream>
#include <cmath>

#include "TVector3.h"

#include "Fiducial.h"

using namespace std;

int main(int argc, char ** argv)
{
  if (argc != 2)
    {
      cerr << "Wrong number of arguments. Instead use: \n"
	   << "\tfiducial_test [E_beam (MeV)]\n\n";
      return -1;
    }

  int E_beam = atoi(argv[1]);
  if (!((E_beam == 4461)||(E_beam == 2261)))
    {
      cerr << "Invalid beam energy provided. Try 4461 or 2261!\n\n";
      return -2;
    }

  Fiducial my_Fid(E_beam,2250,5996,"3He",true);

  // Loop over cin and get p theta phi points
  double p, theta_deg, phi_deg;
  while (cin >> p)
    {
      cin >> theta_deg >> phi_deg;

      double theta = theta_deg*M_PI/180.;
      double phi = phi_deg*M_PI/180.;

      TVector3 mom(p*sin(theta)*cos(phi),p*sin(theta)*sin(phi),p*cos(theta));
      int passedFid = my_Fid.e_inFidRegion(mom)? 1:0;

      cout << p << " " << theta_deg << " " << phi_deg << " " << passedFid << "\n";
    }

  return 0;
}
