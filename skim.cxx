#include "Riostream.h"
#include "TApplication.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TRint.h"
#include "TLatex.h"
#include "TChain.h"
#include "TStyle.h"
#include "TF1.h"
#include "TMath.h"
#include "TVector3.h"
#include "TH3D.h"
#include "TLorentzVector.h"

#include <cstdio>
#include <cstdlib>

using namespace std;




Float_t fgPar_4Gev_2250_Efid_t0_p[2][6];
Float_t fgPar_4Gev_2250_Efid_t1_p[6][6],fgPar_4Gev_2250_Efid_b_p[6][2][6],fgPar_4Gev_2250_Efid_a_p[6][2][6];
Float_t fgPar_4Gev_2250_Phi[6][3],fgPar_4Gev_2250_Theta[6][4];

void SetMomCorrParameters();
void SetFiducialCutParameters();
Bool_t GetEPhiLimits(Float_t momentum, Float_t theta, Int_t sector,
		     Float_t *EPhiMin, Float_t *EPhiMax);
Bool_t EFiducialCut2(TVector3 momentum);
Bool_t PFiducialCut(TVector3 momentum, Float_t *, Float_t *);


const double en_beam = 4.461;
const double fTorusCurrent = 2250;




const double m_n = .9395;
const double m_d = 1.8756;
const double m_p = .9383;
const double E = 4.461;

const double alpha = .007297;


//use Mariana's cuts to get accurate reconstructed data. Use simulation_data.root (after using ./write_simulation on the outfiles, which gives reconstruction according to eg2) to *only* get the generated variables, phi_g and the like.


double rutherford(TVector3 p_e){
  //double E = 5.014;
  //double E = 4.461;

  int ZZ = 1;
  double C = TMath::Pi() * .5 * ZZ * alpha*alpha;
  double hbar_c = 1; //natural units

  return C * (hbar_c / E) * (hbar_c / E) * (1 / ( 1 - p_e.CosTheta() )) * ( 1 / ( 1 - p_e.CosTheta() ) );

}



double mott(TVector3 p_e){
  //double E = 5.014;
  //double E = 4.461;
  
  return alpha*alpha * .5 * (1 + p_e.CosTheta()) / (E*E* (1 - p_e.CosTheta()) * (1 - p_e.CosTheta())
						    * ( 1 + (E * (1 - p_e.CosTheta()) / m_p )      ));
}


double neutrino_CS(TVector3 p_e) {
  //double E = 5.014;
  double G = 1.16638e-5;//GeV ^ -2

  double C = G*G / (TMath::Pi() * TMath::Pi());
  return C * ( E * m_n + (.5 * ( m_n*m_n - m_p*m_p ))) / ( m_n + (E * ( 1 - p_e.CosTheta() ) ))
    * ( E * m_n + (.5 * ( m_n*m_n - m_p*m_p ))) / ( m_n + (E * ( 1 - p_e.CosTheta() ) ));
}



int main(int argc, char ** argv) {
  
  /*#ifdef WITHRINT
  TRint *myapp = new TRint("RootSession", &argc, argv, NULL, 0);
#else
  TApplication *myapp = new TApplication("myapp", 0,0);
  #endif*/

  //TFile *f = new TFile("output/new/simulation_data_e2a.root");
  //TFile *f = new TFile("output/jul25/simulation_data_e2a.root");
  //TFile * f = new TFile("output/simulation_data_eff.root");
  TFile * f = new TFile("output/simulation_data_p.root");
  TTree *t = (TTree*)f->Get("data");



  
  //TChain *t = new TChain("data");
  //for(int i=1; i<=5; i++) {
  //  t->AddFile(Form("output/simulation_data_e2a_%d.root", i));
  //}
  

  TH2D * electron_pid_cuts_p_theta = new TH2D("electron_pid_cuts_p_theta", "Electrons passing PID cuts, p vs #theta", 90, 0, 90, 100, 0, 5);
  electron_pid_cuts_p_theta->GetXaxis()->SetTitle("#theta");
  electron_pid_cuts_p_theta->GetYaxis()->SetTitle("p [GeV/c]");

  TH2D * electron_fid_cuts_p_theta = new TH2D("electron_fid_cuts_p_theta", "Electrons passing PID and fiducial cuts, p vs #theta", 90, 0, 90, 100, 0, 5);
  electron_fid_cuts_p_theta->GetXaxis()->SetTitle("#theta");
  electron_fid_cuts_p_theta->GetYaxis()->SetTitle("p [GeV/c]");

  TH2D * electron_acc_p_theta = new TH2D("electron_acc_p_theta", "Electrons passing PID and fiducial cuts, p vs #theta", 90, 0, 90, 100, 0, 5);
  electron_acc_p_theta->GetXaxis()->SetTitle("#theta");
  electron_acc_p_theta->GetYaxis()->SetTitle("p [GeV/c]");

  
  
  TH2D * electron_pid_cuts_theta_phi = new TH2D("electron_pid_cuts_theta_phi", "Electrons passing PID cuts, #theta vs #phi", 360, 0, 360, 90, 0, 90);
  electron_pid_cuts_theta_phi->GetXaxis()->SetTitle("#phi");
  electron_pid_cuts_theta_phi->GetYaxis()->SetTitle("#theta");
  
  
  TH2D * electron_fid_cuts_theta_phi = new TH2D("electron_fid_cuts_theta_phi", "Electrons passing PID and fiducial cuts, #theta vs #phi", 360, 0, 360, 90, 0, 90);
  electron_fid_cuts_theta_phi->GetXaxis()->SetTitle("#phi");
  electron_fid_cuts_theta_phi->GetYaxis()->SetTitle("#theta");

  TH2D * electron_acc_theta_phi = new TH2D("electron_acc_theta_phi", "Electrons passing PID and fiducial cuts, #theta vs #phi", 360, 0, 360, 90, 0, 90);
  electron_acc_theta_phi->GetXaxis()->SetTitle("#phi");
  electron_acc_theta_phi->GetYaxis()->SetTitle("#theta");

  
  TH2D * electron_pid_cuts_p_phi = new TH2D("electron_pid_cuts_p_phi", "Electrons passing PID cuts, p vs #phi", 360, 0, 360, 100, 0, 5);
  electron_pid_cuts_p_phi->GetXaxis()->SetTitle("#phi");
  electron_pid_cuts_p_phi->GetYaxis()->SetTitle("p [GeV/c]");

  TH2D * electron_fid_cuts_p_phi = new TH2D("electron_fid_cuts_p_phi", "Electrons passing PID and fiducial cuts, p vs #phi", 360, 0, 360, 100, 0, 5);
  electron_fid_cuts_p_phi->GetXaxis()->SetTitle("#phi");
  electron_fid_cuts_p_phi->GetYaxis()->SetTitle("p [GeV/c]");

  TH2D * electron_acc_p_phi = new TH2D("electron_acc_p_phi", "Electrons passing PID and fiducial cuts, p vs #phi", 360, 0, 360, 100, 0, 5);
  electron_acc_p_phi->GetXaxis()->SetTitle("#phi");
  electron_acc_p_phi->GetYaxis()->SetTitle("p [GeV/c]");


  TH3D * electron_fid_cuts_3D = new TH3D("sim_electron_fid_cuts_3D", "Electrons passing PID and fiducial cuts, Reconstructed #phi vs cos(#theta) vs p", 100, 0, 5, 100, 0, 1, 360, 0, 360);
  electron_fid_cuts_3D->GetXaxis()->SetTitle("p [GeV/c]");
  electron_fid_cuts_3D->GetYaxis()->SetTitle("cos(#theta)");
  electron_fid_cuts_3D->GetZaxis()->SetTitle("#phi");

  TH3D * electron_acc_3D = new TH3D("electron_acc_3D", "Electrons passing PID and fiducial cuts, Generated #phi vs cos(#theta) vs p", 100, 0, 5, 100, 0, 1, 360, 0, 360);
  electron_acc_3D->GetXaxis()->SetTitle("p [GeV/c]");
  electron_acc_3D->GetYaxis()->SetTitle("cos(#theta)");
  electron_acc_3D->GetZaxis()->SetTitle("#phi");


  
  TH2D * electron_gen_p_theta = new TH2D("electron_gen_p_theta", "Generated Electrons, p vs #theta", 90, 0, 90, 100, 0, 5);
  electron_gen_p_theta->GetXaxis()->SetTitle("#theta");
  electron_gen_p_theta->GetYaxis()->SetTitle("p [GeV/c]");

  TH2D * electron_gen_theta_phi = new TH2D("electron_gen_theta_phi", "Generated Electrons, #theta vs #phi", 360, 0, 360, 90, 0, 90);
  //Remember to add 30 degrees to phi to make sure the gen and accepted line up- range for accepted is 0 to 360.
  electron_gen_theta_phi->GetXaxis()->SetTitle("#phi");
  electron_gen_theta_phi->GetYaxis()->SetTitle("#theta");
  
  TH2D * electron_gen_p_phi = new TH2D("electron_gen_p_phi", "Generated Electrons, p vs #phi", 360, 0, 360, 100, 0, 5);
  //same note as above
  electron_gen_p_phi->GetXaxis()->SetTitle("#phi");
  electron_gen_p_phi->GetYaxis()->SetTitle("p [GeV/c]");

  TH3D * electron_gen_3D = new TH3D("electron_gen_3D", "Generated Electrons, #phi vs cos(#theta) vs p", 100, 0, 5, 100, 0, 1, 360, 0, 360);
  //same note as above
  electron_gen_3D->GetXaxis()->SetTitle("p [GeV/c]");
  electron_gen_3D->GetYaxis()->SetTitle("cos(#theta)");
  electron_gen_3D->GetZaxis()->SetTitle("#phi");

  
  TH2D * electron_tphi[6];
  for(int i=0; i<6; i++){
    electron_tphi[i] = new TH2D(Form("electron_sim_tphi_%d",i), Form("Slice %d", i+1), 360, 0, 360, 90, 0, 90);
  }

  

  TH2D * proton_fid_cuts_p_theta = new TH2D("proton_fid_cuts_p_theta", "Protons passing fiducial cuts, p vs #theta", 180, 0, 180, 100, 0, 5);
  proton_fid_cuts_p_theta->GetXaxis()->SetTitle("#theta");
  proton_fid_cuts_p_theta->GetYaxis()->SetTitle("p [GeV/c]");

  TH2D * proton_fid_cuts_theta_phi = new TH2D("proton_fid_cuts_theta_phi", "Protons passing fiducial cuts, #theta vs #phi", 360, 0, 360, 180, 0, 180);
  proton_fid_cuts_theta_phi->GetXaxis()->SetTitle("#phi");
  proton_fid_cuts_theta_phi->GetYaxis()->SetTitle("#theta");

  TH2D * proton_fid_cuts_p_phi = new TH2D("proton_fid_cuts_p_phi", "Protons passing fiducial cuts, p vs #phi", 360, 0, 360, 100, 0, 5);
  proton_fid_cuts_p_phi->GetXaxis()->SetTitle("#phi");
  proton_fid_cuts_p_phi->GetYaxis()->SetTitle("p [GeV/c]");

  
  TH2D * proton_acc_p_theta = new TH2D("proton_acc_p_theta", "Protons passing fiducial cuts, p vs #theta", 180, 0, 180, 100, 0, 5);
  proton_acc_p_theta->GetXaxis()->SetTitle("#theta");
  proton_acc_p_theta->GetYaxis()->SetTitle("p [GeV/c]");

  TH2D * proton_acc_theta_phi = new TH2D("proton_acc_theta_phi", "Protons passing fiducial cuts, #theta vs #phi", 360, 0, 360, 180, 0, 180);
  proton_acc_theta_phi->GetXaxis()->SetTitle("#phi");
  proton_acc_theta_phi->GetYaxis()->SetTitle("#theta");

  TH2D * proton_acc_p_phi = new TH2D("proton_acc_p_phi", "Protons passing fiducial cuts, p vs #phi", 360, 0, 360, 100, 0, 5);
  proton_acc_p_phi->GetXaxis()->SetTitle("#phi");
  proton_acc_p_phi->GetYaxis()->SetTitle("p [GeV/c]");




  
  TH3D * proton_fid_cuts_3D = new TH3D("sim_proton_fid_cuts_3D", "Protons passing fiducial cuts, Reconstructed #phi vs cos(#theta) vs p", 100, 0, 5, 200, -1, 1, 360, 0, 360);
  proton_fid_cuts_3D->GetXaxis()->SetTitle("p [GeV/c]");
  proton_fid_cuts_3D->GetYaxis()->SetTitle("cos(#theta)");
  proton_fid_cuts_3D->GetZaxis()->SetTitle("#phi");

  TH3D * proton_acc_3D = new TH3D("proton_acc_3D", "Protons passing fiducial cuts, Generated #phi vs cos(#theta) vs p", 100, 0, 5, 200, -1, 1, 360, 0, 360);
  proton_acc_3D->GetXaxis()->SetTitle("p [GeV/c]");
  proton_acc_3D->GetYaxis()->SetTitle("cos(#theta)");
  proton_acc_3D->GetZaxis()->SetTitle("#phi");
  

  
  TH2D * proton_gen_p_theta = new TH2D("proton_gen_p_theta", "Generated Protons, p vs #theta", 180, 0, 180, 100, 0, 5);
  proton_gen_p_theta->GetXaxis()->SetTitle("#theta");
  proton_gen_p_theta->GetYaxis()->SetTitle("p [GeV/c]");

  TH2D * proton_gen_theta_phi = new TH2D("proton_gen_theta_phi", "Generated Protons, #theta vs #phi", 360, 0, 360, 180, 0, 180);
  
  proton_gen_theta_phi->GetXaxis()->SetTitle("#phi");
  proton_gen_theta_phi->GetYaxis()->SetTitle("#theta");

  TH2D * proton_gen_p_phi = new TH2D("proton_gen_p_phi", "Generated Protons, p vs #phi", 360, 0, 360, 100, 0, 5);
  //same note as above
  proton_gen_p_phi->GetXaxis()->SetTitle("#phi");
  proton_gen_p_phi->GetYaxis()->SetTitle("p [GeV/c]");

  TH3D * proton_gen_3D = new TH3D("proton_gen_3D", "Generated Protons, #phi vs cos(#theta) vs p", 100, 0, 5, 200, -1, 1, 360, 0, 360);
  //same note as above
  proton_gen_3D->GetXaxis()->SetTitle("p [GeV/c]");
  proton_gen_3D->GetYaxis()->SetTitle("cos(#theta)");
  proton_gen_3D->GetZaxis()->SetTitle("#phi");



  /*
  TH2D * electron_fid_cuts_theta_phi_pcut = new TH2D("electron_fid_cuts_theta_phi_pcut", "Electrons passing PID and fiducial cuts, p>1.5 GeV/c, #theta vs #phi", 360, 0, 360, 90, 0, 90);
  electron_fid_cuts_theta_phi_pcut->GetXaxis()->SetTitle("#phi");
  electron_fid_cuts_theta_phi_pcut->GetYaxis()->SetTitle("#theta");
  
  TH2D * electron_gen_theta_phi_pcut = new TH2D("electron_gen_theta_phi_pcut", "Generated Electrons, p>1.5 GeV/c, #theta vs #phi", 360, 0, 360, 90, 0, 90);
  electron_gen_theta_phi_pcut->GetXaxis()->SetTitle("#phi");
  electron_gen_theta_phi_pcut->GetYaxis()->SetTitle("#theta");
  */


  TH2F *h2_prot_px_py_p = new TH2F("h2_prot_px_py_p","",100,-1,1,100,-1,1);
  TH2F *h2_prot_px_py_p_fidcut = new TH2F("h2_prot_px_py_p_fidcut","",100,-1,1,100,-1,1);//copied both from Mariana's code


  
  
  TH1D * electron_zvert = new TH1D("electron_zvert", "Reconstructed Electron Vertex for Events Passing Cuts", 100, -10, 10);
  //TH2D * electron_xyvert = new TH2D("electron_xyvert", "Reconstructed Electron Vertex (x,y) for Events Passing Cuts", 100, -4, 4, 100, -4, 4);

  //TH1D * electron_zvert_gen = new TH1D("electron_zvert_gen", "Generated Electron Vertex", 100, -10, 10);
  //TH2D * electron_xyvert_gen = new TH2D("electron_xyvert_gen", "Generated Electron Vertex", 100, -4, 4, 100, -4, 4);

  TH1D * num_prot = new TH1D("num_prot", "Number of Protons", 4, 0, 4);

  TH1D * proton_deltt_hist = new TH1D("proton_deltt_hist", "#Delta t, in ns", 100, -10, 10);
  int num_passing_cuts = 0;
  int num_passing_fid_cuts = 0;

  
  /*
  //Int_t           dc[40];   //[gpart]
  Int_t           cc[40];   //[gpart]
  Int_t           sc[40];   //[gpart]
  Int_t           ec[40];   //[gpart]

  Float_t         etot[40];   //[ec_part]
  Float_t         ec_ei[40];   //[ec_part]
  Float_t         ec_eo[40];   //[ec_part]
  Float_t         p[40];


  Float_t cx[40];
  Float_t cy[40];
  Float_t cz[40];
  //Int_t cc_segm[40];
  //Int_t cc_sect[40];
  //Float_t cc_t[40];
  //Float_t cc_r[40];
  Float_t sc_r[40];
  Float_t sc_t[40];

  Int_t q[40];
  Float_t vx[40];
  Float_t vy[40];
  Float_t vz[40];

  Int_t gpart;
  Int_t stat[40];
  Int_t id[40];*/
  

  int num_g[40];
  int id_g[40];
  Float_t p_g[40];
  Float_t px_g[40];
  Float_t py_g[40];
  Float_t pz_g[40];

  Float_t z_g[40];
  Float_t theta_g[40];
  Float_t phi_g[40];
  
  Int_t gpart;
  Int_t cc[40];
  Int_t dc[40];
  Int_t ec[40];
  Int_t sc[40];
  Float_t stat[40];
  
  Float_t EC_in[40];
  Float_t EC_out[40];
  Float_t EC_tot[40];
  Float_t q[40];
  Int_t id[40];
  Int_t num[40];
  Float_t p[40], px[40], py[40], pz[40];
  Float_t theta[40];
  Float_t phi[40];
  Float_t z[40];

  Float_t SC_Path[40];
  Float_t SC_Time[40];



  
  double delt_uplim, delt_lowlim;
  Float_t cphil = 0;
  Float_t cphir= 0;

  Float_t tr_time;
  //int el_segment;
  //double el_sccc_timediff;
  /*
  t->SetBranchAddress("cc", cc);
  t->SetBranchAddress("ec", ec);
  t->SetBranchAddress("sc", sc);
  t->SetBranchAddress("sc", sc);
  t->SetBranchAddress("etot", etot);
  t->SetBranchAddress("ec_ei", ec_ei);
  t->SetBranchAddress("ec_eo", ec_eo);
  t->SetBranchAddress("p", p);
  //t->SetBranchAddress("cc_segm", cc_segm);
  //t->SetBranchAddress("cc_sect", cc_sect);
  //t->SetBranchAddress("cc_t", cc_t);
  //t->SetBranchAddress("cc_r", cc_t);  
  t->SetBranchAddress("sc_r", sc_r);
  t->SetBranchAddress("sc_t", sc_t);
  t->SetBranchAddress("q", q);
  t->SetBranchAddress("vx", vx);
  t->SetBranchAddress("vy", vy);
  t->SetBranchAddress("vz", vz);
  t->SetBranchAddress("cx", cx);
  t->SetBranchAddress("cy", cy);
  t->SetBranchAddress("cz", cz);
  t->SetBranchAddress("id", id);
  t->SetBranchAddress("tr_time", &tr_time);*/
  
  t->SetBranchAddress("Number_g", num_g);
  t->SetBranchAddress("particle_g", id_g);
  t->SetBranchAddress("Momentum_g", p_g);
  t->SetBranchAddress("Momentumx_g", px_g);
  t->SetBranchAddress("Momentumy_g", py_g);
  t->SetBranchAddress("Momentumz_g", pz_g);
  t->SetBranchAddress("TargetZ_g", z_g);
  t->SetBranchAddress("Theta_g", theta_g);
  t->SetBranchAddress("Phi_g", phi_g);
  
  t->SetBranchAddress("StatCC", cc);
  t->SetBranchAddress("StatDC", dc);
  t->SetBranchAddress("StatEC", ec);
  t->SetBranchAddress("StatSC", sc);
  t->SetBranchAddress("EC_in", EC_in);
  t->SetBranchAddress("EC_out", EC_out);
  t->SetBranchAddress("EC_tot", EC_tot);
  t->SetBranchAddress("Charge", q);
  t->SetBranchAddress("particle", id);
  t->SetBranchAddress("Momentum", p);
  t->SetBranchAddress("Momentumx", px);
  t->SetBranchAddress("Momentumy", py);
  t->SetBranchAddress("Momentumz", pz);
  t->SetBranchAddress("Theta", theta);
  t->SetBranchAddress("Phi", phi);
  t->SetBranchAddress("TargetZ", z);
  t->SetBranchAddress("gPart", &gpart);
  t->SetBranchAddress("Stat", stat);
  t->SetBranchAddress("SC_Path", SC_Path);
  t->SetBranchAddress("SC_Time", SC_Time);










  TFile *file_in=new TFile("files/el_Epratio_mom4461.root");
  
    
  TF1 *el_Epratio_mean=(TF1*)file_in->Get("f_mean");
  TF1 *el_Epratio_sig=(TF1*)file_in->Get("f_sig");
  
  //proton stuff
  TFile *file_in1 = new TFile("files/protdeltat_mom.root");
  TF1 *prot_deltat_sig=(TF1*)file_in1->Get("sig_pol9");
  TF1 *prot_deltat_mean=(TF1*)file_in1->Get("mean_pol9"); 
  

  double prot_delt_cutrange = 3.;
  double min_good_mom=1.5;//,mom_max_lim=1.6;
  double max_mom=3.8;
    
  double epratio_sig_cutrange=3.;
  double ece; 
  double sum_val=0;
  double sub_val=0;
  double bett, deltt;


  //proton stuff
  double prot_mom_lim=2.15;
  double m_prot=0.938;
  const Double_t c = 2.99792E+10;

  double ns_to_s = 1E-9;
  //int n_elec = 0;
  const int ind_em=0;


  SetFiducialCutParameters();
  SetMomCorrParameters();

  int nEntries = t->GetEntries();
  cout<<"Number of entries: "<<nEntries<<endl;
  cout<<"number with gPart>0="<<t->GetEntries("gPart>0")<<endl;
  for(int i=0; i<nEntries; i++){
    
    t->GetEntry(i);
    if(i%100000==0) cout<<"Processing Event "<<i<<endl;

    
    //generated stuff
    if(p_g[0] > (m_d / (1 + (m_d/E) - TMath::Cos(theta_g[0] * TMath::DegToRad()) ) ) || p_g[0] < 1.5) continue;//elastic maximum for momentum
    //if(p_g[0] < 2 || p_g[0] > 2.1)continue;
    if(phi_g[0]<0) phi_g[0]+=360;
    if(p_g[0] > 1.5){
      electron_gen_p_theta->Fill(theta_g[0], p_g[0]);
      electron_gen_theta_phi->Fill(phi_g[0], theta_g[0]);
      electron_gen_p_phi->Fill(phi_g[0], p_g[0]);

      electron_gen_3D->Fill(p_g[0], TMath::Cos( theta_g[0] * TMath::DegToRad() ), phi_g[0]);
    }

    //if( p_g[0]>1.5 ) electron_gen_theta_phi_pcut->Fill(phi_g[0], theta_g[0]);
    //end generated stuff
    
    //reconstructed stuff
    //ece = TMath::Max( ec_ei[ec[ind_em] - 1] + ec_eo[ec[ind_em] - 1], etot[ec[ind_em] - 1]);
    ece = TMath::Max(EC_in[0] + EC_out[0], EC_tot[0]);
  
    //el_segment=int((cc_segm[cc[ind_em]-1]-int(cc_segm[cc[ind_em]-1]/1000)*1000)/10);
    //el_cc_sector=cc_sect[cc[ind_em]-1];
    //el_sccc_timediff=sc_t[cc[ind_em]-1]-cc_t[cc[ind_em]-1]-(sc_r[cc[ind_em]-1]-cc_r[cc[ind_em]-1])/(c*ns_to_s);

    if(p[ind_em]>min_good_mom && p[ind_em]<max_mom) sum_val=el_Epratio_mean->Eval(p[ind_em])+epratio_sig_cutrange*el_Epratio_sig->Eval(p[ind_em]);
    if(p[ind_em]>max_mom) sum_val= el_Epratio_mean->Eval(max_mom)+epratio_sig_cutrange*el_Epratio_sig->Eval(max_mom);
      
    if(p[ind_em]>min_good_mom && p[ind_em]<max_mom) sub_val= el_Epratio_mean->Eval(p[ind_em])-epratio_sig_cutrange*el_Epratio_sig->Eval(p[ind_em]);
    if(p[ind_em]>max_mom) sub_val= el_Epratio_mean->Eval(max_mom)-epratio_sig_cutrange*el_Epratio_sig->Eval(max_mom);
   


    if( gpart > 0 && ec[ind_em] > 0/*0.5*/ && cc[ind_em] > 0/* 0.5*/ &&  sc[ind_em] > 0/*0.5*/ && q[ind_em] < 0 &&


	//ec_ei[ec[ind_em] - 1] > 0.055
	EC_in[0] > 0.055


	&& ece>0.33

	&& ece/p[ind_em]>sub_val


	&& ece/p[ind_em] <sum_val
	

	&& p[ind_em]>min_good_mom

	&& z[0] > 4.1 && z[0] < 6
	&& p_g[0] > 1.5
	) //electron pid cuts (no need for CC cut at 4.46 GeV)
      {
	num_passing_cuts++;
	//fill a histogram for electrons passing pid cuts
	Float_t cx = px[0] / p[0];
	Float_t cy = py[0] / p[0];
	Float_t cz = pz[0] / p[0];
	TVector3 el_mom1(px[ind_em], py[ind_em], pz[ind_em]);

	double m = mott(el_mom1);
	//double n = neutrino_CS(el_mom1);
	
	//TVector3 el_mom1(p[ind_em]*cx[ind_em], p[ind_em]*cy[ind_em], p[ind_em]*cz[ind_em]);
	
	
	double el_theta = TMath::ACos(cz)*TMath::RadToDeg();
	double el_phi_mod = TMath::ATan2(cy, cx)*TMath::RadToDeg();
	if(el_phi_mod<0)el_phi_mod+=360;
	float pmag = p[ind_em];
	
	electron_pid_cuts_p_theta->Fill(el_theta, pmag);
	electron_pid_cuts_theta_phi->Fill(el_phi_mod, el_theta);
	electron_pid_cuts_p_phi->Fill(el_phi_mod, pmag);
	
       
	
	if(pmag > 1.5 && pmag < 2){
	  electron_tphi[0]->Fill(el_phi_mod, el_theta);
	}else if(pmag > 2 && pmag < 2.5){
	  electron_tphi[1]->Fill(el_phi_mod, el_theta);
	}else if(pmag > 2.5 && pmag < 3){
	  electron_tphi[2]->Fill(el_phi_mod, el_theta);
	}else if(pmag > 3 && pmag < 3.5){
	  electron_tphi[3]->Fill(el_phi_mod, el_theta);
	}else if(pmag > 3.5 && pmag < 4){
	  electron_tphi[4]->Fill(el_phi_mod, el_theta);
	}else if(pmag > 4 && pmag < 4.5){
	  electron_tphi[5]->Fill(el_phi_mod, el_theta);
	}
	  
	
	
	
	
	
	if(EFiducialCut2(el_mom1)) {
	  electron_zvert->Fill(z[0]);
	  //electron_xyvert->Fill(vx[0], vy[0]);
	  num_passing_fid_cuts++;

	  
	  electron_fid_cuts_p_theta->Fill(el_theta, pmag);
	  electron_fid_cuts_theta_phi->Fill(el_phi_mod, el_theta);
	  electron_fid_cuts_p_phi->Fill(el_phi_mod, pmag);

	  //if(pmag>1.5) electron_fid_cuts_theta_phi_pcut->Fill(el_phi_mod, el_theta);
	  
	  electron_fid_cuts_3D->Fill(pmag, cz, el_phi_mod);//cz equivalent to cos(theta_e)


	  
	  electron_acc_p_theta->Fill(theta_g[0], p_g[0]);
	  electron_acc_theta_phi->Fill(phi_g[0], theta_g[0]);
	  electron_acc_p_phi->Fill(phi_g[0], p_g[0]);

	  
	  electron_acc_3D->Fill(p_g[0], TMath::Cos( theta_g[0] * TMath::DegToRad() ), phi_g[0]); 



	  //proton stuff

	  //generated stuff
	  if(phi_g[1] < 0) phi_g[1]+=360;
	  
	  
	  proton_gen_p_theta->Fill(theta_g[1], p_g[1]);
	  proton_gen_theta_phi->Fill(phi_g[1], theta_g[1]);
	  proton_gen_p_phi->Fill(phi_g[1], p_g[1]);

	  //if( p_g[1]>1.5 ) electron_gen_theta_phi_pcut->Fill(phi_g[1], theta_g[1]);
    
	  proton_gen_3D->Fill(p_g[1], TMath::Cos( theta_g[1] * TMath::DegToRad() ), phi_g[1]);


	  //end generated stuff









	  
	  int index_p[20],ind_p,index_pi[20],ind_pi;
	  int num_p=0;
	  //int num_pi=0,num_pimi=0,num_pipl=0;
	  // int index_n[20]={},ec_index_n[20]={},index_pipl[20]={},index_pimi[20]={};
	  //int num_n = 0,ec_num_n = 0;
	  //double pi_phi,pi_phi_mod, pi_theta,pimi_phi,pimi_phi_mod,pimi_theta,pipl_phi,pipl_phi_mod, pipl_theta;

	  for( int i = 1; i < TMath::Min(gpart, 20); i++ )
	    {
	      if( sc[i] > 0 && stat[i] > 0 &&  id[i] == 2212 )
		{
		  ind_p=i;

		  bett=p[ind_p]/TMath::Sqrt(p[ind_p]*p[ind_p]+m_prot*m_prot);

		  Float_t tr_time = SC_Time[0]-SC_Path[0]/(c*ns_to_s);//beta=1 for electron
		    //deltt=sc_t[sc[ind_p]-1]-sc_r[sc[ind_p]-1]/(bett*c*ns_to_s) - tr_time;               
		  deltt = SC_Time[ind_p] - SC_Path[ind_p]/(bett*c*ns_to_s) - tr_time;

		  proton_deltt_hist->Fill(deltt);

		  
                
		  if(p[ind_p]<prot_mom_lim) delt_uplim=prot_deltat_mean->Eval(p[ind_p])+prot_delt_cutrange*prot_deltat_sig->Eval(p[ind_p]);
		  if(p[ind_p]>prot_mom_lim) delt_uplim=prot_deltat_mean->Eval(prot_mom_lim)+prot_delt_cutrange*prot_deltat_sig->Eval(prot_mom_lim);
                
                
		  if(p[ind_p]<prot_mom_lim) delt_lowlim= prot_deltat_mean->Eval(p[ind_p])-prot_delt_cutrange*prot_deltat_sig->Eval(p[ind_p]);
		  if(p[ind_p]>prot_mom_lim) delt_lowlim= prot_deltat_mean->Eval(prot_mom_lim)-prot_delt_cutrange*prot_deltat_sig->Eval(prot_mom_lim);
		  //prot_delt_p->Fill(p[ind_p], deltt);
		  if(deltt<delt_uplim || deltt>delt_lowlim){  //PID cut
		    //prot_delt_p_pidcut->Fill(p[ind_p], deltt);
		    TLorentzVector V4_uncorrprot(px[ind_p], py[ind_p], pz[ind_p], /*p[ind_p]*cx,p[ind_p]*cy[ind_p],p[ind_p]*cz[ind_p],*/TMath::Sqrt(m_prot*m_prot+p[ind_p]*p[ind_p]));
		    h2_prot_px_py_p->Fill(px[ind_p] / p[ind_p], py[ind_p] / p[ind_p]);


		    if(PFiducialCut(V4_uncorrprot.Vect(), &cphil, &cphir)){ //proton fiducial cuts
		      h2_prot_px_py_p_fidcut->Fill(px[ind_p] / p[ind_p], py[ind_p] / p[ind_p]);//(cx[ind_p],cy[ind_p]);
         

		      num_p = num_p + 1;
		      
		      index_p[num_p-1]=i;
		    }
		  }


		}
	    }
	  num_prot->Fill(num_p);//it's only ever 0 or 1
	  if(num_p==1) {

	    ind_p = index_p[0];

	    //int ind_p=1;
	    double cost = pz[ind_p] / p[ind_p];
	    double prot_theta = TMath::ACos(cost)*TMath::RadToDeg();
	    //double el_theta = TMath::ACos(cz[ind_p])*TMath::RadToDeg();
	    //double el_phi_mod = TMath::ATan2(cy[ind_p], cx[ind_p])*TMath::RadToDeg();
	    double prot_phi_mod = TMath::ATan2(py[ind_p] / p[ind_p],  px[ind_p] / p[ind_p])*TMath::RadToDeg();
	    if(prot_phi_mod<0)prot_phi_mod+=360;
	    //float pmag = sqrt( p[0]*p[0] + p[1]*p[1] + p[2]*p[2] );
	    float pmag = p[ind_p];

	    
	    proton_fid_cuts_p_theta->Fill(prot_theta, pmag);
	    proton_fid_cuts_theta_phi->Fill(prot_phi_mod, prot_theta);
	    proton_fid_cuts_p_phi->Fill(prot_phi_mod, pmag);
	    proton_fid_cuts_3D->Fill(pmag, cost, prot_phi_mod);

	    proton_acc_p_theta->Fill(theta_g[1], p_g[1]);
	    proton_acc_theta_phi->Fill(phi_g[1], theta_g[1]);
	    proton_acc_p_phi->Fill(phi_g[1], p_g[1]);
	    proton_acc_3D->Fill(p_g[1], TMath::Cos( theta_g[1] * TMath::DegToRad() ), phi_g[1]);
	    //cout<<"Accepted proton: "<<endl<<"p: "<<p_g[1]<<endl<<"cos(theta): "<<TMath::Cos( theta_g[1] * TMath::DegToRad())<<endl<<"phi: "<<phi_g[1]<<endl<<endl;

	  }
	  
	  //end proton stuff

	  
	  

	  
	  
	  	  
	}

	

      }
 
    
  }



  TH2D * electron_rec_ratio_p_theta = (TH2D*)electron_fid_cuts_p_theta->Clone("electron_rec_ratio_p_theta");
  electron_rec_ratio_p_theta->SetTitle("Electron Acceptance Rec_Ratio, p vs #theta");
  electron_rec_ratio_p_theta->Divide(electron_gen_p_theta);

  TH2D * electron_rec_ratio_theta_phi = (TH2D*)electron_fid_cuts_theta_phi->Clone("electron_rec_ratio_theta_phi");
  electron_rec_ratio_theta_phi->SetTitle("Electron Acceptance Rec_Ratio, #theta vs #phi");
  electron_rec_ratio_theta_phi->Divide(electron_gen_theta_phi);

  TH2D * electron_rec_ratio_p_phi = (TH2D*)electron_fid_cuts_p_phi->Clone("electron_rec_ratio_p_phi");
  electron_rec_ratio_p_phi->SetTitle("Electron Acceptance Rec_Ratio, p vs #phi");
  electron_rec_ratio_p_phi->Divide(electron_gen_p_phi);

  
  TH2D * electron_acc_ratio_p_theta = (TH2D*)electron_acc_p_theta->Clone("electron_acc_ratio_p_theta");
  electron_acc_ratio_p_theta->SetTitle("Electron Acceptance Acc_Ratio, p vs #theta");
  electron_acc_ratio_p_theta->Divide(electron_gen_p_theta);

  TH2D * electron_acc_ratio_theta_phi = (TH2D*)electron_acc_theta_phi->Clone("electron_acc_ratio_theta_phi");
  electron_acc_ratio_theta_phi->SetTitle("Electron Acceptance Acc_Ratio, #theta vs #phi");
  electron_acc_ratio_theta_phi->Divide(electron_gen_theta_phi);

  TH2D * electron_acc_ratio_p_phi = (TH2D*)electron_acc_p_phi->Clone("electron_acc_ratio_p_phi");
  electron_acc_ratio_p_phi->SetTitle("Electron Acceptance Acc_Ratio, p vs #phi");
  electron_acc_ratio_p_phi->Divide(electron_gen_p_phi);

  //TH2D * electron_ratio_theta_phi_pcut = (TH2D*)electron_fid_cuts_theta_phi_pcut->Clone("electron_ratio_theta_phi_pcut");
  //electron_ratio_theta_phi_pcut->SetTitle("Electron Acceptance Ratio, #theta vs #phi [p > 1.5 GeV/c]");
  //electron_ratio_theta_phi_pcut->Divide(electron_gen_theta_phi_pcut);

  //TH3D * electron_ratio_3D = (TH3D*)electron_fid_cuts_3D->Clone("electron_ratio_3D");
  TH3D * electron_ratio_3D = (TH3D*)electron_acc_3D->Clone("electron_ratio_3D");
  electron_ratio_3D->SetTitle("Electron Acceptance Ratio, #phi vs cos(#theta) vs p");
  electron_ratio_3D->Divide(electron_gen_3D);



  TH2D * proton_rec_ratio_p_theta = (TH2D*)proton_fid_cuts_p_theta->Clone("proton_rec_ratio_p_theta");
  proton_rec_ratio_p_theta->SetTitle("Proton Acceptance Rec_Ratio, p vs #theta");
  proton_rec_ratio_p_theta->Divide(proton_gen_p_theta);

  TH2D * proton_rec_ratio_theta_phi = (TH2D*)proton_fid_cuts_theta_phi->Clone("proton_rec_ratio_theta_phi");
  proton_rec_ratio_theta_phi->SetTitle("Proton Acceptance Rec_Ratio, #theta vs #phi");
  proton_rec_ratio_theta_phi->Divide(proton_gen_theta_phi);

  TH2D * proton_rec_ratio_p_phi = (TH2D*)proton_fid_cuts_p_phi->Clone("proton_rec_ratio_p_phi");
  proton_rec_ratio_p_phi->SetTitle("Proton Acceptance Rec_Ratio, p vs #phi");
  proton_rec_ratio_p_phi->Divide(proton_gen_p_phi);


  TH2D * proton_acc_ratio_p_theta = (TH2D*)proton_acc_p_theta->Clone("proton_acc_ratio_p_theta");
  proton_acc_ratio_p_theta->SetTitle("Proton Acceptance Acc_Ratio, p vs #theta");
  proton_acc_ratio_p_theta->Divide(proton_gen_p_theta);

  TH2D * proton_acc_ratio_theta_phi = (TH2D*)proton_acc_theta_phi->Clone("proton_acc_ratio_theta_phi");
  proton_acc_ratio_theta_phi->SetTitle("Proton Acceptance Acc_Ratio, #theta vs #phi");
  proton_acc_ratio_theta_phi->Divide(proton_gen_theta_phi);

  TH2D * proton_acc_ratio_p_phi = (TH2D*)proton_acc_p_phi->Clone("proton_acc_ratio_p_phi");
  proton_acc_ratio_p_phi->SetTitle("Proton Acceptance Acc_Ratio, p vs #phi");
  proton_acc_ratio_p_phi->Divide(proton_gen_p_phi);

  
  //TH3D * proton_ratio_3D = (TH3D*)proton_fid_cuts_3D->Clone("proton_ratio_3D");
  TH3D * proton_ratio_3D = (TH3D*)proton_acc_3D->Clone("proton_ratio_3D");
  proton_ratio_3D->SetTitle("Proton Acceptance Ratio, #phi vs cos(#theta) vs p");
  proton_ratio_3D->Divide(proton_gen_3D);




  
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1");
  electron_gen_p_theta->Draw("colz");

  TCanvas *c2 = new TCanvas("c2");
  electron_fid_cuts_p_theta->Draw("colz");
  TCanvas *c2a = new TCanvas("c2a");
  electron_pid_cuts_p_theta->Draw("colz");

  TCanvas *c3 = new TCanvas("c3");
  electron_rec_ratio_p_theta->Draw("colz");

  TCanvas *c3a = new TCanvas("c3a");
  electron_acc_ratio_p_theta->Draw("colz");



  TCanvas *c4 = new TCanvas("c4");
  electron_gen_theta_phi->Draw("colz");

  TCanvas *c5 = new TCanvas("c5");
  electron_fid_cuts_theta_phi->Draw("colz");
  TCanvas *c5a = new TCanvas("c5a");
  electron_pid_cuts_theta_phi->Draw("colz");

  TCanvas *c6 = new TCanvas("c6");
  electron_rec_ratio_theta_phi->Draw("colz");

  TCanvas *c6a = new TCanvas("c6a");
  electron_acc_ratio_theta_phi->Draw("colz");

  //TCanvas *c4a = new TCanvas("c4a");
  //electron_gen_theta_phi_pcut->Draw("colz");

  //TCanvas *c5a = new TCanvas("c5a");
  //electron_fid_cuts_theta_phi_pcut->Draw("colz");

  //TCanvas *c6a = new TCanvas("c6a");
  //electron_ratio_theta_phi_pcut->Draw("colz"); //originally end comment here......


  TCanvas *c7 = new TCanvas("c7");
  electron_gen_p_phi->Draw("colz");

  TCanvas *c8 = new TCanvas("c8");
  electron_fid_cuts_p_phi->Draw("colz");
  TCanvas *c8a = new TCanvas("c8a");
  electron_pid_cuts_p_phi->Draw("colz");

  TCanvas *c9 = new TCanvas("c9");
  electron_rec_ratio_p_phi->Draw("colz");

  TCanvas * c9a = new TCanvas("c9a");
  electron_acc_ratio_p_phi->Draw("colz");
  

  c1->Print("plots/new_ratios.pdf[");
  c1->Print("plots/new_ratios.pdf");
  c2->Print("plots/new_ratios.pdf");
  c2a->Print("plots/new_ratios.pdf");
  c3->Print("plots/new_ratios.pdf");
  c3a->Print("plots/new_ratios.pdf");
  c4->Print("plots/new_ratios.pdf");
  c5->Print("plots/new_ratios.pdf");
  c5a->Print("plots/new_ratios.pdf");
  c6->Print("plots/new_ratios.pdf");
  c6a->Print("plots/new_ratios.pdf");
  
  //c4a->Print("plots/new_ratios.pdf");
  //c5a->Print("plots/new_ratios.pdf");
  //c6a->Print("plots/new_ratios.pdf");
  //end comment
  c7->Print("plots/new_ratios.pdf");
  c8->Print("plots/new_ratios.pdf");
  c8a->Print("plots/new_ratios.pdf");
  c9->Print("plots/new_ratios.pdf");
  c9a->Print("plots/new_ratios.pdf");
  c9a->Print("plots/new_ratios.pdf]");
  
  TCanvas *c10 = new TCanvas("c10");
  num_prot->Draw();
  c10->Print("plots/new_num_prot.pdf");








  
  TCanvas *p1 = new TCanvas("p1");
  proton_gen_p_theta->Draw("colz");

  TCanvas *p2 = new TCanvas("p2");
  proton_fid_cuts_p_theta->Draw("colz");

  TCanvas *p3 = new TCanvas("p3");
  proton_rec_ratio_p_theta->Draw("colz");

  TCanvas * p3a = new TCanvas("p3a");
  proton_acc_ratio_p_theta->Draw("colz");

  
  TCanvas *p4 = new TCanvas("p4");
  proton_gen_theta_phi->Draw("colz");

  TCanvas *p5 = new TCanvas("p5");
  proton_fid_cuts_theta_phi->Draw("colz");

  TCanvas *p6 = new TCanvas("p6");
  proton_rec_ratio_theta_phi->Draw("colz");

  TCanvas * p6a = new TCanvas("p6a");
  proton_acc_ratio_theta_phi->Draw("colz");
  

  
  TCanvas *p7 = new TCanvas("p7");
  proton_gen_p_phi->Draw("colz");

  TCanvas *p8 = new TCanvas("p8");
  proton_fid_cuts_p_phi->Draw("colz");

  TCanvas *p9 = new TCanvas("p9");
  proton_rec_ratio_p_phi->Draw("colz");

  TCanvas *p9a = new TCanvas("p9a");
  proton_acc_ratio_p_phi->Draw("colz");

  p1->Print("plots/new_p_ratios.pdf[");
  p1->Print("plots/new_p_ratios.pdf");
  p2->Print("plots/new_p_ratios.pdf");
  p3->Print("plots/new_p_ratios.pdf");
  p3a->Print("plots/new_p_ratios.pdf");
  p4->Print("plots/new_p_ratios.pdf");
  p5->Print("plots/new_p_ratios.pdf");
  p6->Print("plots/new_p_ratios.pdf");
  p6a->Print("plots/new_p_ratios.pdf");
  p7->Print("plots/new_p_ratios.pdf");
  p8->Print("plots/new_p_ratios.pdf");
  p9->Print("plots/new_p_ratios.pdf");
  p9a->Print("plots/new_p_ratios.pdf");
  p9a->Print("plots/new_p_ratios.pdf]");
  

  TCanvas * p10 = new TCanvas("p10");
  proton_deltt_hist->Draw();
  p10->Print("plots/deltt_hist.pdf");

  
  
  //TFile *rootfile = new TFile("e2a_carbon_acc.root", "RECREATE");
  /*TFile *rootfile = new TFile("jul25.root", "RECREATE");
  electron_gen_3D->Write();
  electron_fid_cuts_3D->Write();
  electron_acc_3D->Write();
  
  proton_gen_3D->Write();
  proton_fid_cuts_3D->Write();
  proton_acc_3D->Write();
  rootfile->Close();*/

  TFile *rootfile2 = new TFile("sim.root", "RECREATE");
  electron_fid_cuts_3D->Write();
  proton_fid_cuts_3D->Write();
  rootfile2->Close();

  TFile *p_slices_file = new TFile("p_slices_sim.root", "RECREATE");
  for(int k=0; k<6; k++){
    electron_tphi[k]->Write();
  }
  p_slices_file->Close();


  
  //myapp->Run();
  f->Close();
  return 0;
}







void SetFiducialCutParameters(){
  // reads from a file the parameters of the fiducial cut functions
  // Please refer to <A HREF="http://einstein.unh.edu/protopop/FiducialCuts/fc4E2.html">Fiducial Cuts</A> -- D.Protopopescu(UNH)

  if(en_beam>4. && en_beam<5.){
    //
    // reads FC parameters for 4.4GeV
    //
    FILE *fiducial_par;
    Int_t ptype; // Parameter Type (see encoding in URL above)
    Int_t ci;    // sector number
    Float_t par[6];
    Char_t Filename[100], ParLocation[80];
    sprintf(ParLocation,"files"); // where do we put this parameter file ??
    //
    printf("Reading fiducial cut parameters for 4.4GeV/2250A ...\n");
    sprintf(Filename,"%s/FCP_4GeV.par", ParLocation);
    fiducial_par = fopen(Filename,"r");
    while(!feof(fiducial_par)){
      fscanf(fiducial_par, "%i   %i   %f   %f   %f   %f  %f   %f", &ptype, &ci, &par[0], &par[1], &par[2], &par[3], &par[4], &par[5]);
      switch (ptype){
      case  0: for(int k=0; k<2; k++) fgPar_4Gev_2250_Efid_t0_p[ci-1][k] = par[k];// power function for t0
	break;
      case  1: for(int k=0; k<6; k++) fgPar_4Gev_2250_Efid_t1_p[ci-1][k] = par[k];// poly5 for the others
	break;
      case 10: for(int k=0; k<6; k++) fgPar_4Gev_2250_Efid_b_p[ci-1][0][k] = par[k];
	break;
      case 11: for(int k=0; k<6; k++) fgPar_4Gev_2250_Efid_b_p[ci-1][1][k] = par[k];
	break;
      case 20: for(int k=0; k<6; k++) fgPar_4Gev_2250_Efid_a_p[ci-1][0][k] = par[k];
	break;
      case 21: for(int k=0; k<6; k++) fgPar_4Gev_2250_Efid_a_p[ci-1][1][k] = par[k];
	break;
      default: printf("error!\n");
      }
    }
    fclose(fiducial_par);
  }
  else printf("Don't know how to set fiducial cut parameters for %3.1f GeV!\n", en_beam);

}


Bool_t GetEPhiLimits(Float_t momentum, Float_t theta, Int_t sector,Float_t *EPhiMin, Float_t *EPhiMax){
  //Begin_Html
  /*</pre>
    Information for electron fiducial cut, 
    returns the minimum and maximum phi accepted for a given momentum, theta and sector
    momentum is in GeV/c, theta is in degrees, 0 <= sector <= 5
    EPhiMin and EPhiMax are in degrees
    Function returns False if inputs are out of bounds
    1.1 GeV not implemented yet
    tested against EFiducialCut to make sure the limits are identical
    2.2 GeV: tested for 10 < theta < 65, -30 < phi < 360, 0.1 < Ef < 2.261
    2 inconsistent events out of 10^6
    4.4 GeV: tested for 10 < theta < 65, -30 < phi < 360, 0.3 < Ef < 4.461
    0 inconsistent events out of 10^6
    Please refer to <A HREF="http://www.jlab.org/Hall-B/secure/e2/bzh/efiducialcut.html">Electron Fiducial Cuts</A> -- Bin Zhang (MIT).
    For 4.4GeV please refer to <A HREF="http://einstein.unh.edu/protopop/FiducialCuts/fc4E2.html">Fiducial Cuts</A> -- D.Protopopescu (UNH)
    <pre>
  */
  //End_Html
  if (sector < 0 || sector > 5) return kFALSE;    // bad input

  if(en_beam>4. && en_beam<5. && fTorusCurrent>2240. && fTorusCurrent<2260.){// 4.4GeV fiducial cuts by protopop@jlab.org
    if ((theta < 15.) || (momentum < 0.9)) return kFALSE;         // out of range
    Float_t t0, t1, b[2], a[2];

    if (momentum > 3.7) momentum = 3.7; // don't extrapolate past the data


    // uncomment this if you want 100MeV energy bins
    //Enrgy = 0.100*int(Enrgy/0.100);


    // calculates parameters of cut functions for this energy
    t0 = fgPar_4Gev_2250_Efid_t0_p[sector][0]/pow(momentum, fgPar_4Gev_2250_Efid_t0_p[sector][1]);
    t1 = 0.; for(int k=0; k<6; k++) t1 += (fgPar_4Gev_2250_Efid_t1_p[sector][k]*pow(momentum, k));
    for(int l=0; l<2; l++){
      b[l] = 0.; for(int k=0; k<6; k++) b[l] += (fgPar_4Gev_2250_Efid_b_p[sector][l][k]*pow(momentum, k));
      a[l] = 0.; for(int k=0; k<6; k++) a[l] += (fgPar_4Gev_2250_Efid_a_p[sector][l][k]*pow(momentum, k));
    }

    // adjust upper limit according to hardware
    if(t1 < 45.) t1 = 45.;
    if(t0 < theta && theta < t1){

      *EPhiMin = 60.*sector - b[0]*(1. - 1/((theta - t0)/(b[0]/a[0]) + 1.));
      *EPhiMax = 60.*sector + b[1]*(1. - 1/((theta - t0)/(b[1]/a[1]) + 1.));
      // if(momentum<1.65 && momentum>1.60)cout<<sector<<"  "<<a[0]<<"    "<<a[1]<<"    "<<a[2]<<endl;
    }
    else {
      *EPhiMin = 60.*sector;
      *EPhiMax = 60.*sector;
    }


  }   // 4.4 GeV e2a
  else {
    return kFALSE;     // wrong beam energy/torus
  }
  return kTRUE;
}




Bool_t EFiducialCut2(TVector3 momentum){
  //Begin_Html
  /*</pre>
    Electron fiducial cut, return kTRUE if the electron is in the fiducial volume
    modified 14 May 2001 lbw
    Now calls GetEPhiLimits for 2.2 and 4.4 GeV 
    tested against EFiducialCut for both 2.2 (with and without bad scintillator cuts) and 4.4 GeV
    discrepancy less than 2 in 10^6 events
    Please refer to <A HREF="http://www.jlab.org/Hall-B/secure/e2/bzh/efiducialcut.html">Electron Fiducial Cuts</A> -- Bin Zhang (MIT).
    For 4.4GeV please refer to <A HREF="http://einstein.unh.edu/protopop/FiducialCuts/fc4E2.html">Fiducial Cuts</A> -- D.Protopopescu (UNH)
    Please refer to <a href="http://www.jlab.org/Hall-B/secure/e2/stevenmc/FiducialCuts/index.html">1.1 GeV fiducial cuts</a> -- Steven McLauchlan (GU).
    <pre>
  */
  //End_Html

  Bool_t status = kTRUE;
  Float_t phiMin, phiMax;
  Float_t mom = momentum.Mag();
  Float_t phi = momentum.Phi()*180./TMath::Pi();
  if(phi<-30.) phi += 360.;
  Float_t theta = momentum.Theta()*180./TMath::Pi();
  Int_t  sector = (Int_t)((phi+30.)/60.);
  if(sector < 0) sector = 0;
  if(sector > 5) sector = 5; // to match array index
  // all the work is now done in GetEPhiLimits
  status = GetEPhiLimits(mom, theta, sector, &phiMin, &phiMax);

  if (status) {
    status = status && (phi > phiMin) && (phi < phiMax);
  }
  return status;
}
























Bool_t PFiducialCut(TVector3 momentum, Float_t *ptr_cphil, Float_t *ptr_cphir){
  //Positive Hadron Fiducial Cut
  //Please refer to <A HREF="http://www.jlab.org/Hall-B/secure/e2/bzh/pfiducialcut.html">Electron Fiducial Cuts</A> -- Bin Zhang (MIT).

  //copied from the header file
  const Float_t fgPar_2Gev_2250_Phi[6][3] = {{0.9903, -3.926E-4, 1.318E-5}, {0.9803, -3.177E-4, 1.706E-5}, {0.9922, 1.836E-4, 1.474E-5}, {0.9898,
																	  5.259E-5, 1.45E-5}, {0.9906, -1.604E-4, 1.66E-5}, {0.9925, 1.902E-4, 1.985E-5}};
  const Float_t fgPar_2Gev_2250_Theta[3][4] = {{8.56526E-4, 7.89140E+1, -8.41321E-1, 1.00082}, {6.10625E-1, 8.30600E-1, -4.40544E-1, 0}, {-5.02481, 1.29011E+1, -6.90397, 0}};






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
//end header material
 

  
  Bool_t status = kTRUE;
  if (en_beam>4. && en_beam<5. && fTorusCurrent>2240. && fTorusCurrent<2260.){//4 GeV Fiducial Cut Rustam Niyazov
    Float_t phi=momentum.Phi()*180/TMath::Pi(); if(phi<-30) phi+=360;
    Int_t sector = Int_t ((phi+30)/60); if(sector<0)sector=0; if(sector>5) sector=5;
    phi -= sector*60;
    Float_t theta = momentum.Theta()*180/TMath::Pi();
    Float_t p = momentum.Mag();

    Float_t parfidl[3];for(Int_t i=0; i<3; i++){parfidl[i]=0;} 
    Float_t parfidr[3];for(Int_t i=0; i<3; i++){parfidr[i]=0;}
    Float_t parfidbl[2];for(Int_t i=0; i<2; i++){parfidbl[i]=0;}
    Float_t parfidbr[2];for(Int_t i=0; i<2; i++){parfidbr[i]=0;}
    Float_t cphil=0; Float_t cphir=0;
    Float_t phi45l=0; Float_t phi45r=0;
    Float_t phi60l=0; Float_t phi60r=0;
    Float_t theta_min=11;



    //bool SCpdcut=true;
    bool Forward=kFALSE; //defines if particle in Forward (Forward=kTRUE) or Backward (Forward=kFALSE) region.
    Int_t thetab=45; //this variable defines the edge point for Forward<->Backward regions
    Float_t p1=0.575; //last bin momentum for region p<0.6 GeV/c

    Float_t theta_max=140;
    if(p<0.2)p=0.2; //momentum less than 0.2 GeV/c, use 0.2 GeV/c
    if(p>4.4)p=4.4; //momentum greater than 4.4 GeV/c, use 4.4 GeV/c

    //get parametrized values of theta_max for p<0.6 GeV/c region
    if(p<0.6){theta_max=fgPar_4Gev_2250_Pfidbl[sector][4]+fgPar_4Gev_2250_Pfidbl[sector][5]*p;}
    //get parametrized values of theta_max for p>0.6 GeV/c region
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
    phi45l=parfidl[0]*(parfidl[2]-45)/(45-parfidl[2]+(parfidl[1]/parfidl[0])); //parametrized value of phi at theta=45 deg.
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
      
    if(Forward){//Forward region
      if(p<0.6) theta_min=14; else theta_min=11;//for p<0.6 GeV/c Region theta starts from 14 deg., otherwise 11 deg.   
      cphil=parfidl[0]*(parfidl[2]-theta)/(theta-parfidl[2]+(parfidl[1]/parfidl[0]));//hyperbolic function
      cphir=-parfidr[0]*(parfidr[2]-theta)/(theta-parfidr[2]+(parfidr[1]/parfidr[0]));
    }
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

    if(phi<0) status=(phi>cphil); //check the constrains 
    else if(phi>=0) {status=(phi<cphir);
    }

    if(theta<theta_min) status=kFALSE; //Cutting out events below theta_min

    if(Forward && p<0.6 && theta<20.6-11.4*p)status=kFALSE; //function defines cut of the edge at low theta for p<0.6 GeV/c

    //p>0.6 GeV/c. Cut of the edge at low theta  for some sectors and for 
    //some range of momentum, where edge does not look good.
    bool s1s4=(theta<11.7&&(sector==0||sector==3));
    bool s5=(theta<12.2&&sector==4);
    bool s6=(theta<11.4&&sector==5);
    if(p>=0.6&&p<1.5&&(s1s4||s5||s6)) status=kFALSE; 
    //cout<<"cphil"<<cphir<<endl;

    *ptr_cphil = cphil;
    *ptr_cphir = cphir;


  }
  //  cout<<"cphil"<<cphir<<endl;
  //
  return status;
}


void SetMomCorrParameters(){
//reads from a file and calculates fiducial cut parameters 
//Please refer to 
// <A HREF="http://einstein.unh.edu/protopop/MomentumCorrections/emc4E2.html">Momentum Corrections</A> -- D.Protopopescu(UNH)

 if( en_beam>4. &&  en_beam<5.){
   //
   // reads momcorr parameters for 4.4GeV 
   //
   FILE *momcorr_par; 
   Int_t ptype; // Parameter Type = 0 for Phi function, 1 for Theta function
   Int_t ci, cj;
   Float_t par_val;
   Char_t Filename[100], ParLocation[80]; 
   sprintf(ParLocation,"files");// where do we put the parameter file ??
   //
   //   printf("Reading momentum corrections for 4.4GeV/2250A ...\n");
   sprintf(Filename,"%s/EMCP_4GeV.par", ParLocation);
   momcorr_par = fopen(Filename,"r");
   while(!feof(momcorr_par)){
     fscanf(momcorr_par, "%i    %i    %i    %e \n", &ptype, &ci, &cj, &par_val);
     switch (ptype){
      case 0: fgPar_4Gev_2250_Phi[ci-1][cj] = par_val;
        //printf("Sector %i b[%i] = %e \n", ci, cj, fgPar_4Gev_2250_Phi[ci-1][cj]);
        break;
      case 1: fgPar_4Gev_2250_Theta[ci-1][cj] = par_val; 
        //printf("Sector %i d[%i] = %e \n", ci, cj, fgPar_4Gev_2250_Theta[ci-1][cj]);
        break;
      default: printf("error!\n"); 
     }
   }
   fclose(momcorr_par);
 }
 else printf("Don't know how to set momentum corrections parameters for %3.1f GeV!\n",  en_beam); 
}

