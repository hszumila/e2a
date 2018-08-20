#include "TH3.h"
#include "TH2.h"
#include "TH1.h"
#include "TFile.h"

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

int main(int argc, char ** argv)
{
  if (argc != 3)
    {
      cerr << "Wrong number of arguments. Instead use:\n"
	   << "   basic_test /path/to/map/file /path/to/output.root\n\n";
      return -1;
    }

  TFile * infile = new TFile(argv[1]);
  TFile * outfile = new TFile(argv[2],"RECREATE");

  TH3D * gen = (TH3D*) infile->Get("Generated Particles");
  TH3D * acc = (TH3D*) infile->Get("Accepted Particles");
  TH1D * badbins = new TH1D("bins_under_10","bins_under_10",12,-1.5,10.5);
  // Tests
  // Double check axis ranges and bins
  // All bins have 10 generated
  // 0 momentum electrons should have no acceptance
  // 5 GeV electrons should have great acceptance
 
  int nPbins, nCbins, nFbins;
  nPbins = gen->GetXaxis()->GetNbins();
  nCbins = gen->GetYaxis()->GetNbins();
  nFbins = gen->GetZaxis()->GetNbins();
  cout << "X axis:\n\t";
  cout << nPbins << " bins, from " << gen->GetXaxis()->GetBinLowEdge(1) << " to " << gen->GetXaxis()->GetBinUpEdge(nPbins) << "\n";
  cout << "Y axis:\n\t";
  cout << nCbins << " bins, from " << gen->GetYaxis()->GetBinLowEdge(1) << " to " << gen->GetYaxis()->GetBinUpEdge(nCbins) << "\n";
  cout << "Z axis:\n\t";
  cout << nFbins << " bins, from " << gen->GetZaxis()->GetBinLowEdge(1) << " to " << gen->GetZaxis()->GetBinUpEdge(nFbins) << "\n";

  std::vector<int> bad_bins;
  for (int pBin=1 ; pBin<= nPbins; pBin++)
    for (int cBin=1 ; cBin <=nCbins ; cBin++)
      for (int fBin=1 ; fBin<=nFbins ; fBin++)
	{
	  int bin = gen->GetBin(pBin,cBin,fBin);
	  if (gen->GetBinContent(bin) <10.)
      {
	    bad_bins.push_back(bin);
      badbins->Fill(gen->GetBinContent(bin));
      }
  }
  cout << "\n\t There are " << bad_bins.size() << " bins with fewer than 10 generated events\n";
  
  gen->GetXaxis()->SetRange(1,1);
  acc->GetXaxis()->SetRange(1,1);
  TH2D * gen_P0 = (TH2D*) gen->Project3D("zy");
  TH2D * acc_P0 = (TH2D*) acc->Project3D("zy");

  gen->GetXaxis()->SetRange(100,100);
  acc->GetXaxis()->SetRange(100,100);
  TH2D * gen_P5 = (TH2D*) gen->Project3D("zy");
  TH2D * acc_P5 = (TH2D*) acc->Project3D("zy");

  gen_P0->Write();
  acc_P0->Write();
  gen_P5->Write();
  acc_P5->Write();
  badbins->Write();
  outfile->Close();

  return 0;
}
