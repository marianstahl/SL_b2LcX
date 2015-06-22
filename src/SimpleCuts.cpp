/*

    2015-05-25 M. Stahl

    Determine cut and background efficiencies for a given cut string

*/

#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <exception>

#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TObject.h"
#include "TStopwatch.h"
#include "Riostream.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TProfile.h"

#include "TArrow.h"
#include "TLatex.h"

#include "RooRealVar.h"
#include "RooDataHist.h"

#include "../include/Fit_Charm_Baryon.h"
#include "../include/configuration.h"
#include "../include/MyStyle.h"

using namespace std;
using namespace RooFit;

TString temp;

struct cH_mult_cand{
  double mass;
  float ip;
  int index;
};

struct Xiccp_mult_cand{
  double mass;
  float h1ip;
  float h2ip;
  int index;
};

struct Xiccpp_mult_cand{
  double mass;
  float h1ip;
  float h2ip;
  float h3ip;
  int index;
};

bool compareByIP(const cH_mult_cand &a, const cH_mult_cand &b){return a.ip < b.ip;}
bool compareByH1IP(const Xiccp_mult_cand &a, const Xiccp_mult_cand &b){return a.h1ip < b.h1ip;}
bool compareByH2IP(const Xiccp_mult_cand &a, const Xiccp_mult_cand &b){return a.h2ip < b.h2ip;}
bool compareXByH1IP(const Xiccpp_mult_cand &a, const Xiccpp_mult_cand &b){return a.h1ip < b.h1ip;}
bool compareXByH2IP(const Xiccpp_mult_cand &a, const Xiccpp_mult_cand &b){return a.h2ip < b.h2ip;}
bool compareXByH3IP(const Xiccpp_mult_cand &a, const Xiccpp_mult_cand &b){return a.h3ip < b.h3ip;}

template <class cH> void fill_profile(vector<cH> multiple_cH_candidate, TProfile *Xc_IP);
void SimpleCuts();
template<class T> void make_1D_Plot(T *hist, configuration* myconfig);
void make_2D_Plot(TH2D* hist, configuration* myconfig);
void make_overlay_Plot(TH1D *SS_hist, TH1D *OS_hist, configuration *myconfig);

int main(int argc, char **argv)
{
  SimpleCuts();
  return 0;
}

void SimpleCuts(){

  TStopwatch *clock = new TStopwatch();
  clock->Start(1);

  bool Lc_sel = false;
  /*bool Xiccp_sel = true;
  bool Xiccpp_sel = true;
  bool make_IP_profile_plots = false;*/
  bool clean_vertex = false;
  bool purge_vertex = false;
  bool fit_charm_hadron = true;

  const double protonmass = 938.272013; //MeV
  const double pionmass = 139.57018; //MeV
  const double kaonmass = 493.677; //MeV
  const double lcmass = 2286.46;
  const double lcPDGmass = 2286.46;

  MyStyle();

  configuration* myconfig = new configuration();
  myconfig->set_version(1);
  myconfig->fill_cs("Basic selection");
  myconfig->set_current_cs("IPCut");
  myconfig->set_particle("Lc");

  TH1D::AddDirectory(0);
  TProfile::AddDirectory(0);
  TH2D::AddDirectory(0);

  TH1D *Lc_hist = new TH1D("Lc",";;",myconfig->get_NBinsForFit(),myconfig->get_XcMlo(),myconfig->get_XcMhi());
  TH2D *Lc_beta_hist = new TH2D("Lc_beta",";#beta_{p};M_{inv}(pK#pi) (GeV)",150,-0.75,0.75,2*myconfig->get_NBinsForFit(),myconfig->get_XcMlo(),myconfig->get_XcMhi());
  TH2D *Lc_beta_cut_hist = new TH2D("Lc_beta_cut",";#beta_{p};M_{inv}(pK#pi) (GeV)",150,-0.75,0.75,2*myconfig->get_NBinsForFit(),myconfig->get_XcMlo(),myconfig->get_XcMhi());
  TH2D *Lc_2D_dump = new TH2D("Lc_p_CosTheta",";M_{inv}(pK#pi) (GeV);#Theta^{CM}_{p}",2*myconfig->get_NBinsForFit(),myconfig->get_XcMlo(),myconfig->get_XcMhi(),100,-1,1);

  int n_CorrM_bins = 84;
  double CorrM_lo = 3000, CorrM_hi = 7200;
  TH1D *Lb0_hist = new TH1D("Lambdab0",";M_{corr}(#Lambda_{c}^{+}_{c}#mu^{-} + c.c.) (GeV);Events/50 MeV",n_CorrM_bins,CorrM_lo,CorrM_hi);
  TH1D *Lb0_SB_hist = new TH1D("Lambdab0_SB",";M_{corr}(#Lambda_{c}^{+}_{c}#mu^{-}+ c.c.) (GeV);Events/50 MeV",n_CorrM_bins,CorrM_lo,CorrM_hi);

  /*int n_profilebins = 20;
  TProfile *Lc_IP = new TProfile("Lc_IP",";IP (mm);<nTracks>",n_profilebins,0,1," ");*/

  double Lcpilo = 2430, Lcpihi = 2830; int nbinsLcpi = 200;
  temp.Form(";M_{inv}(#Lambda_{c}^{+}#pi)-M_{inv}(#Lambda_{c}^{+})+M_{PDG}(#Lambda_{c}^{+}) (GeV); Events/%g MeV",(Lcpihi-Lcpilo)/nbinsLcpi);
  TH1D *Lcpi_hist = new TH1D("Lcpi_hist",temp,nbinsLcpi,Lcpilo,Lcpihi);
  TH1D *LcpiWS_hist = new TH1D("LcpiWS_hist",temp,nbinsLcpi,Lcpilo,Lcpihi);
  double LcKlo = 2780, LcKhi = 3380; int nbinsLcK = 150;
  temp.Form(";M_{inv}(#Lambda_{c}^{+}K)-M_{inv}(#Lambda_{c}^{+})+M_{PDG}(#Lambda_{c}^{+}) (GeV); Events/%g MeV",(LcKhi-LcKlo)/nbinsLcK);
  TH1D *LcK_hist = new TH1D("LcK_hist",temp,nbinsLcK,LcKlo,LcKhi);
  TH1D *LcKWS_hist = new TH1D("LcKWS_hist",temp,nbinsLcK,LcKlo,LcKhi);
  vector<TH1D*> Lc_hists;Lc_hists.push_back(Lcpi_hist);Lc_hists.push_back(LcpiWS_hist);Lc_hists.push_back(LcK_hist);Lc_hists.push_back(LcKWS_hist);

  double ScKlo = 2950, ScKhi = 3950; int nbinsScK = 200;
  temp.Form(";M_{inv}(#Sigma_{c}(2455)^{++}K)-M_{inv}(#Lambda_{c}^{+}#pi^{+})+M_{PDG}(#Sigma_{c}(2455)^{++} + c.c.) (GeV); Events/%g MeV",(ScKhi-ScKlo)/nbinsScK);
  TH1D *Xiccp_hist = new TH1D("Xiccp_hist",temp,nbinsScK,ScKlo,ScKhi);
  TH1D *XiccpWC_hist = new TH1D("XiccpWC_hist",temp,nbinsScK,ScKlo,ScKhi);
  temp.Form(";M_{inv}(#Sigma_{c}(2455)^{++}K#pi^{+})-M_{inv}(#Lambda_{c}^{+}#pi^{+})+M_{PDG}(#Sigma_{c}(2455)^{++} + c.c.) (GeV); Events/%g MeV",(ScKhi-ScKlo)/nbinsScK);
  TH1D *Xiccpp_hist = new TH1D("Xiccpp_hist",temp,nbinsScK,ScKlo+140,ScKhi+140);
  TH1D *XiccppWC_hist = new TH1D("XiccppWC_hist",temp,nbinsScK,ScKlo+140,ScKhi+140);

  /*TH1D *Xiccpp_mult = new TH1D("Xiccpp_mult",";multiplicity;# candidates",20,0,20);
  TH1D *Xiccp_mult = new TH1D("Xiccp_mult",";multiplicity;# candidates",20,0,20);*/

  temp = myconfig->get_tupledir() + "SLLcStrp21.root";
  gErrorIgnoreLevel = kError;
  TFile *fSLBS = new TFile(temp,"read");
  TTree *Lc_tree = (TTree*)gDirectory->Get("Lc/DecayTree");
  gErrorIgnoreLevel = kPrint;


  myconfig->set_particle("Lc");
  fSLBS->cd();
  double p_PT, p_ETA, p_PHI, p_CosTheta;
  double K_PT, K_ETA, K_PHI;
  double pi_PT, pi_ETA, pi_PHI;
  double Xb_OWNPV_X, Xb_OWNPV_Y, Xb_OWNPV_Z;
  double Xb_ENDVERTEX_X, Xb_ENDVERTEX_Y, Xb_ENDVERTEX_Z;
  double Xb_PT, Xb_ETA, Xb_PHI, Xb_M;
  double Xc_PT, Xc_ETA, Xc_PHI, Xc_M;
  double p_ProbNNp, K_ProbNNk, pi_ProbNNpi;
  float Added_H_PT[200], Added_H_ETA[200], Added_H_PHI[200], Added_CharmH_M[200], Added_CharmH_M_kaon[200], Added_CharmH_PT[200], Added_CharmH_VERTEXCHI2_NEW[200];
  int Added_n_Particles, Xc_ID;
  float Added_H_PROBNNPID[200], Added_H_ProbNNpi[200], Added_H_ProbNNk[200];
  float b_IPCHI2[200], PV_IPCHI2[200];

  if(Lc_sel){Lc_tree->SetBranchStatus("*",0); //disable all branches
  //now switch on the ones we need (saves a lot of time)
  Lc_tree->SetBranchStatus("Lb_M",1);
  Lc_tree->SetBranchStatus("Lb_PT",1);
  Lc_tree->SetBranchStatus("Lb_ETA",1);
  Lc_tree->SetBranchStatus("Lb_PHI",1);
  Lc_tree->SetBranchStatus("Lb_OWNPV_X",1);
  Lc_tree->SetBranchStatus("Lb_OWNPV_Y",1);
  Lc_tree->SetBranchStatus("Lb_OWNPV_Z",1);
  Lc_tree->SetBranchStatus("Lb_ENDVERTEX_X",1);
  Lc_tree->SetBranchStatus("Lb_ENDVERTEX_Y",1);
  Lc_tree->SetBranchStatus("Lb_ENDVERTEX_Z",1);

  Lc_tree->SetBranchStatus("Lc_M",1);
  Lc_tree->SetBranchStatus("Lc_ID",1);
  Lc_tree->SetBranchStatus("Lc_PT",1);
  Lc_tree->SetBranchStatus("Lc_ETA",1);
  Lc_tree->SetBranchStatus("Lc_PHI",1);

  Lc_tree->SetBranchStatus("Added_n_Particles",1);
  Lc_tree->SetBranchStatus("Added_H_PT",1);
  Lc_tree->SetBranchStatus("Added_H_ETA",1);
  Lc_tree->SetBranchStatus("Added_H_PHI",1);
  Lc_tree->SetBranchStatus("Added_H_PROBNNPID",1);
  Lc_tree->SetBranchStatus("Added_H_ProbNNpi",1);
  Lc_tree->SetBranchStatus("Added_H_ProbNNk",1);

  Lc_tree->SetBranchStatus("Added_CharmH_M",1);
  Lc_tree->SetBranchStatus("Added_CharmH_M_kaon",1);
  Lc_tree->SetBranchStatus("Added_CharmH_PT",1);

  Lc_tree->SetBranchStatus("Added_CharmH_LOGIPCHI2_NEW",1);
  Lc_tree->SetBranchStatus("Added_CharmH_LOGMINIPCHI2",1);
  Lc_tree->SetBranchStatus("Added_CharmH_VERTEXCHI2_NEW",1);

  Lc_tree->SetBranchStatus("p_ProbNNp",1);
  Lc_tree->SetBranchStatus("K_ProbNNk",1);
  Lc_tree->SetBranchStatus("pi_ProbNNpi",1);
  Lc_tree->SetBranchStatus("p_PT",1);
  Lc_tree->SetBranchStatus("p_CosTheta",1);
  Lc_tree->SetBranchStatus("p_ETA",1);
  Lc_tree->SetBranchStatus("p_PHI",1);
  Lc_tree->SetBranchStatus("K_PT",1);
  Lc_tree->SetBranchStatus("K_ETA",1);
  Lc_tree->SetBranchStatus("K_PHI",1);
  Lc_tree->SetBranchStatus("pi_PT",1);
  Lc_tree->SetBranchStatus("pi_ETA",1);
  Lc_tree->SetBranchStatus("pi_PHI",1);

  //set the branch addresses
  Lc_tree->SetBranchAddress("Lb_M",&Xb_M);
  Lc_tree->SetBranchAddress("Lb_PT",&Xb_PT);
  Lc_tree->SetBranchAddress("Lb_ETA",&Xb_ETA);
  Lc_tree->SetBranchAddress("Lb_PHI",&Xb_PHI);
  Lc_tree->SetBranchAddress("Lb_OWNPV_X",&Xb_OWNPV_X);
  Lc_tree->SetBranchAddress("Lb_OWNPV_Y",&Xb_OWNPV_Y);
  Lc_tree->SetBranchAddress("Lb_OWNPV_Z",&Xb_OWNPV_Z);
  Lc_tree->SetBranchAddress("Lb_ENDVERTEX_X",&Xb_ENDVERTEX_X);
  Lc_tree->SetBranchAddress("Lb_ENDVERTEX_Y",&Xb_ENDVERTEX_Y);
  Lc_tree->SetBranchAddress("Lb_ENDVERTEX_Z",&Xb_ENDVERTEX_Z);

  Lc_tree->SetBranchAddress("Lc_M",&Xc_M);
  Lc_tree->SetBranchAddress("Lc_ID",&Xc_ID);
  Lc_tree->SetBranchAddress("Lc_PT",&Xc_PT);
  Lc_tree->SetBranchAddress("Lc_ETA",&Xc_ETA);
  Lc_tree->SetBranchAddress("Lc_PHI",&Xc_PHI);

  Lc_tree->SetBranchAddress("Added_n_Particles",&Added_n_Particles);
  Lc_tree->SetBranchAddress("Added_H_PT",&Added_H_PT);
  Lc_tree->SetBranchAddress("Added_H_ETA",&Added_H_ETA);
  Lc_tree->SetBranchAddress("Added_H_PHI",&Added_H_PHI);
  Lc_tree->SetBranchAddress("Added_H_PROBNNPID",&Added_H_PROBNNPID);
  Lc_tree->SetBranchAddress("Added_H_ProbNNpi",&Added_H_ProbNNpi);
  Lc_tree->SetBranchAddress("Added_H_ProbNNk",&Added_H_ProbNNk);

  Lc_tree->SetBranchAddress("Added_CharmH_M",&Added_CharmH_M);
  Lc_tree->SetBranchAddress("Added_CharmH_M_kaon",&Added_CharmH_M_kaon);
  Lc_tree->SetBranchAddress("Added_CharmH_PT",&Added_CharmH_PT);
  Lc_tree->SetBranchAddress("Added_CharmH_LOGIPCHI2_NEW",&b_IPCHI2);
  Lc_tree->SetBranchAddress("Added_CharmH_LOGMINIPCHI2",&PV_IPCHI2);
  Lc_tree->SetBranchAddress("Added_CharmH_VERTEXCHI2_NEW",&Added_CharmH_VERTEXCHI2_NEW);

  Lc_tree->SetBranchAddress("p_ProbNNp",&p_ProbNNp);
  Lc_tree->SetBranchAddress("K_ProbNNk",&K_ProbNNk);
  Lc_tree->SetBranchAddress("pi_ProbNNpi",&pi_ProbNNpi);

  Lc_tree->SetBranchAddress("p_PT",&p_PT);
  Lc_tree->SetBranchAddress("p_ETA",&p_ETA);
  Lc_tree->SetBranchAddress("p_PHI",&p_PHI);
  Lc_tree->SetBranchAddress("p_CosTheta",&p_CosTheta);
  Lc_tree->SetBranchAddress("K_PT",&K_PT);
  Lc_tree->SetBranchAddress("K_ETA",&K_ETA);
  Lc_tree->SetBranchAddress("K_PHI",&K_PHI);
  Lc_tree->SetBranchAddress("pi_PT",&pi_PT);
  Lc_tree->SetBranchAddress("pi_ETA",&pi_ETA);
  Lc_tree->SetBranchAddress("pi_PHI",&pi_PHI);

  double p_beta;
  double p_as_piKpi_M, p_as_KKpi_M;
  RooRealVar Charm_mass("Xc_M","M_{inv}(pK^{-}pi^{+} + c.c.)",myconfig->get_XcMlo(),myconfig->get_XcMhi(),"MeV") ;

  UInt_t Lc_nevents = Lc_tree->GetEntries();
  cout << "Entries in " + myconfig->get_particle() + " tree: " << Lc_nevents << endl;
  for (UInt_t evt = 0; evt < Lc_nevents;evt++) {
    Lc_tree->GetEntry(evt);

    TLorentzVector proton;
    proton.SetPtEtaPhiM(p_PT,p_ETA,p_PHI,protonmass);
    TLorentzVector kaon;
    kaon.SetPtEtaPhiM(K_PT,K_ETA,K_PHI,kaonmass);
    TLorentzVector pion;
    pion.SetPtEtaPhiM(pi_PT,pi_ETA,pi_PHI,pionmass);

    p_beta  = (-proton.P()+kaon.P()+pion.P())/(proton.P()+kaon.P()+pion.P());

    TLorentzVector p_as_pi;
    p_as_pi.SetVectM(proton.Vect(),pionmass);
    TLorentzVector p_as_K;
    p_as_K.SetVectM(proton.Vect(),kaonmass);

    p_as_piKpi_M = (p_as_pi + kaon + pion).M();
    p_as_KKpi_M = (p_as_K + kaon + pion).M();

    TLorentzVector Xc;
    Xc.SetPtEtaPhiM(Xc_PT,Xc_ETA,Xc_PHI,Xc_M);

    //some cuts before we begin
    bool vetos = !(1860 < p_as_piKpi_M && p_as_piKpi_M < 1880) && !(1860 < p_as_KKpi_M && p_as_KKpi_M < 1880) && !(1955 < p_as_KKpi_M && p_as_KKpi_M < 1985);//&& !(2005 < MisID_D2piKpi_M && MisID_D2piKpi_M < 2025)
    bool PID_cuts = p_ProbNNp > 0.05 && K_ProbNNk > 0.05 && pi_ProbNNpi > 0.05;
    //bool CorrM_cut = 4000 < Xb_M && Xb_M < 7200;
    bool Lc_Mass_cut = fabs(Xc_M - lcmass) < 20;

    Lc_beta_hist->Fill(p_beta,Xc_M);
    if(vetos && PID_cuts){
      Lc_hist->Fill(Xc_M);
      if(Lc_Mass_cut)Lb0_hist->Fill(Xb_M);
      else Lb0_SB_hist->Fill(Xb_M);
    }

    if(!(Lc_Mass_cut && vetos && PID_cuts))continue;
    Lc_2D_dump->Fill(Xc_M,p_CosTheta);
    Lc_beta_cut_hist->Fill(p_beta,Xc_M);

    vector<cH_mult_cand> multiple_cpi_SS;
    vector<cH_mult_cand> multiple_cpi_OS;
    vector<cH_mult_cand> multiple_cK_SS;
    vector<cH_mult_cand> multiple_cK_OS;
    vector < vector<cH_mult_cand> > multiples;

    for(int ap = 0; ap < Added_n_Particles; ap++){
      bool OS = ((Added_H_PROBNNPID[ap] < 0 && Xc_ID > 0) || (Added_H_PROBNNPID[ap] > 0 && Xc_ID < 0));
      bool added_H_cuts = Added_H_PT[ap] > 180 && Added_CharmH_VERTEXCHI2_NEW[ap] < 6 &&  Added_CharmH_PT[ap] > 1500 && PV_IPCHI2[ap] > 0.6;
      bool added_K_cuts = Added_H_ProbNNk[ap] > 0.1 && Added_CharmH_M_kaon[ap]-Xc_M+lcPDGmass < LcKhi && TMath::Abs(Added_H_PROBNNPID[ap]) == 321;
      bool added_pi_cuts = Added_H_ProbNNpi[ap] > 0.1 && Added_CharmH_M[ap]-Xc_M+lcPDGmass < Lcpihi && TMath::Abs(Added_H_PROBNNPID[ap]) == 221;
      bool piSS_selection = !OS && added_H_cuts && added_pi_cuts;//basic_selection &&
      bool piOS_selection =  OS && added_H_cuts && added_pi_cuts;
      bool KSS_selection  = !OS && added_H_cuts && added_K_cuts;
      bool KOS_selection  =  OS && added_H_cuts && added_K_cuts;
      if(piSS_selection)multiple_cpi_SS.emplace_back(cH_mult_cand{Added_CharmH_M[ap]-Xc_M+lcPDGmass,b_IPCHI2[ap],ap});
      if(piOS_selection)multiple_cpi_OS.emplace_back(cH_mult_cand{Added_CharmH_M[ap]-Xc_M+lcPDGmass,b_IPCHI2[ap],ap});
      if(KSS_selection)multiple_cK_SS.emplace_back(cH_mult_cand{Added_CharmH_M_kaon[ap]-Xc_M+lcPDGmass,b_IPCHI2[ap],ap});
      if(KOS_selection)multiple_cK_OS.emplace_back(cH_mult_cand{Added_CharmH_M_kaon[ap]-Xc_M+lcPDGmass,b_IPCHI2[ap],ap});
    }
    multiples.push_back(multiple_cpi_SS);
    multiples.push_back(multiple_cpi_OS);
    multiples.push_back(multiple_cK_SS);
    multiples.push_back(multiple_cK_OS);

    for(unsigned int i = 0; i < multiples.size(); i++){
      if(multiples.at(i).size() > 1){
        if(myconfig->get_verbosity() > normal)cout << "we have a multiple candidate" << endl;
        sort(multiples.at(i).begin(), multiples.at(i).end(), compareByIP);
        //if(i == 1 && make_IP_profile_plots)fill_profile(multiples.at(i),Xic_IP);
        for(unsigned int j = 0; j < multiples.at(i).size(); j++){
          if(myconfig->get_verbosity() > normal)cout << "IP of sorted candidate # " << j << ": " << multiples.at(i).at(j).ip << endl;
        }
      }
      if(clean_vertex){
        if(multiples.at(i).size() > 1){
          int good_IP_particles = 1;
          int ih1 = 0, ih2 = 0;
          for(unsigned int j = 1; j < multiples.at(i).size(); j++){
            //eliminate Xc*->Xc+nH (where n > 1)
            //Assuming that if n candidates have similar IP, they come from the same vertex
            if(multiples.at(i).at(1).ip/multiples.at(i).at(0).ip < 1.25){
              good_IP_particles++;
              ih1 = multiples.at(i).at(0).index; ih2 = multiples.at(i).at(1).index;
              multiples.at(i).erase(multiples.at(i).begin());//it will always be the first element I hope ?!
              if(purge_vertex)multiples.at(i).clear();
            }
          }
        }
        if(!multiples.at(i).empty())Lc_hists.at(i)->Fill(multiples.at(i).at(0).mass);
      }
      else{
        if(!multiples.at(i).empty())Lc_hists.at(i)->Fill(multiples.at(i).at(0).mass);
        /*if(multiples.at(i).size() > 1){
          for(unsigned int j = 1; j < multiples.at(i).size(); j++){
            Xic_IPfails.at(i)->Fill(multiples.at(i).at(j).mass);
          }*/
      }
    }
  }

  make_overlay_Plot(Lb0_SB_hist,Lb0_hist,myconfig);
  //if(make_IP_profile_plots)make_1D_Plot(Lc_IP,myconfig);
  make_2D_Plot(Lc_beta_hist,myconfig);
  make_2D_Plot(Lc_beta_cut_hist,myconfig);
  make_2D_Plot(Lc_2D_dump,myconfig);
  make_overlay_Plot(Lcpi_hist,LcpiWS_hist,myconfig);
  make_overlay_Plot(LcK_hist,LcKWS_hist,myconfig);


  /*
  make_1D_Plot(Xiccpp_mult,myconfig);
  make_1D_Plot(Xiccp_mult,myconfig);*/

  /*make_overlay_Plot(LcpiSS_IPfail,LcpiOS_IPfail,myconfig);
  make_overlay_Plot(LcKSS_IPfail,LcKOS_IPfail,myconfig);*/
  if(fit_charm_hadron){
    RooDataHist *Lc = new RooDataHist("Xc","Xc",Charm_mass,Import(*Lc_hist)) ;
    Fit_Charm(Lc,myconfig);
  }

  clock->Stop();clock->Print();
  //Lc_tree->SetDirectory(0);
  //fSLBS->Close();
  cout << "SUCCESS" << endl;

  }

  gErrorIgnoreLevel = kError;
  temp = myconfig->get_tupledir() + "SLLcStrp21_friend.root";
  TFile *fXicc = new TFile(temp,"read");
  TTree *Xicc_tree = (TTree*)gDirectory->Get("Lc");
  TTree *XiccWS_tree = (TTree*)gDirectory->Get("LcWS");
  gErrorIgnoreLevel = kPrint;

  float Lcpi1_M[100], Lcpi2_M[100], LcKpi_M[100], LcKpipi_M[100];
  float K2_MINIPCHI2[100], pi2_MINIPCHI2[100], pi3_MINIPCHI2[100], K2_BIPCHI2[100], pi2_BIPCHI2[100], pi3_BIPCHI2[100];
  float Lc3_M[100];
  float K2_ProbNNk[100], pi2_ProbNNpi[100], pi3_ProbNNpi[100];
  int Xicc_cand;

  float Lcpi1_MWS[100], Lcpi2_MWS[100], LcKpi_MWS[100], LcKpipi_MWS[100];
  float K2_MINIPCHI2WS[100], pi2_MINIPCHI2WS[100], pi3_MINIPCHI2WS[100], K2_BIPCHI2WS[100], pi2_BIPCHI2WS[100], pi3_BIPCHI2WS[100];
  float Lc3_MWS[100];
  float K2_ProbNNkWS[100], pi2_ProbNNpiWS[100], pi3_ProbNNpiWS[100];
  double XcWS_M;
  int Xicc_WScand;

  Xicc_tree->SetBranchStatus("Xicc_cand",1);
  Xicc_tree->SetBranchStatus("Lc_M",1);
  Xicc_tree->SetBranchStatus("Lcpi1_M",1);
  Xicc_tree->SetBranchStatus("Lcpi2_M",1);
  Xicc_tree->SetBranchStatus("LcKpi_M",1);
  Xicc_tree->SetBranchStatus("LcKpipi_M",1);
  Xicc_tree->SetBranchStatus("K2_MINIPCHI2",1);
  Xicc_tree->SetBranchStatus("pi2_MINIPCHI2",1);
  Xicc_tree->SetBranchStatus("pi3_MINIPCHI2",1);
  Xicc_tree->SetBranchStatus("K2_BIPCHI2",1);
  Xicc_tree->SetBranchStatus("pi2_BIPCHI2",1);
  Xicc_tree->SetBranchStatus("pi3_BIPCHI2",1);
  Xicc_tree->SetBranchStatus("Lc3_M",1);
  Xicc_tree->SetBranchStatus("K2_ProbNNk",1);
  Xicc_tree->SetBranchStatus("pi2_ProbNNpi",1);
  Xicc_tree->SetBranchStatus("pi3_ProbNNpi",1);
  //set the branch addresses
  Xicc_tree->SetBranchAddress("Xicc_cand",&Xicc_cand);
  Xicc_tree->SetBranchAddress("Lc_M",&Xc_M);
  Xicc_tree->SetBranchAddress("Lcpi1_M",&Lcpi1_M);
  Xicc_tree->SetBranchAddress("Lcpi2_M",&Lcpi2_M);
  Xicc_tree->SetBranchAddress("LcKpi_M",&LcKpi_M);
  Xicc_tree->SetBranchAddress("LcKpipi_M",&LcKpipi_M);
  Xicc_tree->SetBranchAddress("K2_MINIPCHI2",&K2_MINIPCHI2);
  Xicc_tree->SetBranchAddress("pi2_MINIPCHI2",&pi2_MINIPCHI2);
  Xicc_tree->SetBranchAddress("pi3_MINIPCHI2",&pi3_MINIPCHI2);
  Xicc_tree->SetBranchAddress("K2_BIPCHI2",&K2_BIPCHI2);
  Xicc_tree->SetBranchAddress("pi2_BIPCHI2",&pi2_BIPCHI2);
  Xicc_tree->SetBranchAddress("pi3_BIPCHI2",&pi3_BIPCHI2);
  Xicc_tree->SetBranchAddress("Lc3_M",&Lc3_M);
  Xicc_tree->SetBranchAddress("K2_ProbNNk",&K2_ProbNNk);
  Xicc_tree->SetBranchAddress("pi2_ProbNNpi",&pi2_ProbNNpi);
  Xicc_tree->SetBranchAddress("pi3_ProbNNpi",&pi3_ProbNNpi);

  XiccWS_tree->SetBranchStatus("Xicc_cand",1);
  XiccWS_tree->SetBranchStatus("Lc_M",1);
  XiccWS_tree->SetBranchStatus("Lcpi1_M",1);
  XiccWS_tree->SetBranchStatus("Lcpi2_M",1);
  XiccWS_tree->SetBranchStatus("LcKpi_M",1);
  XiccWS_tree->SetBranchStatus("LcKpipi_M",1);
  XiccWS_tree->SetBranchStatus("K2_MINIPCHI2",1);
  XiccWS_tree->SetBranchStatus("pi2_MINIPCHI2",1);
  XiccWS_tree->SetBranchStatus("pi3_MINIPCHI2",1);
  XiccWS_tree->SetBranchStatus("K2_BIPCHI2",1);
  XiccWS_tree->SetBranchStatus("pi2_BIPCHI2",1);
  XiccWS_tree->SetBranchStatus("pi3_BIPCHI2",1);
  XiccWS_tree->SetBranchStatus("Lc3_M",1);
  XiccWS_tree->SetBranchStatus("K2_ProbNNk",1);
  XiccWS_tree->SetBranchStatus("pi2_ProbNNpi",1);
  XiccWS_tree->SetBranchStatus("pi3_ProbNNpi",1);
  //set the branch addresses
  XiccWS_tree->SetBranchAddress("Xicc_cand",&Xicc_WScand);
  XiccWS_tree->SetBranchAddress("Lc_M",&XcWS_M);
  XiccWS_tree->SetBranchAddress("Lcpi1_M",&Lcpi1_MWS);
  XiccWS_tree->SetBranchAddress("Lcpi2_M",&Lcpi2_MWS);
  XiccWS_tree->SetBranchAddress("LcKpi_M",&LcKpi_MWS);
  XiccWS_tree->SetBranchAddress("LcKpipi_M",&LcKpipi_MWS);
  XiccWS_tree->SetBranchAddress("K2_MINIPCHI2",&K2_MINIPCHI2WS);
  XiccWS_tree->SetBranchAddress("pi2_MINIPCHI2",&pi2_MINIPCHI2WS);
  XiccWS_tree->SetBranchAddress("pi3_MINIPCHI2",&pi3_MINIPCHI2WS);
  XiccWS_tree->SetBranchAddress("K2_BIPCHI2",&K2_BIPCHI2WS);
  XiccWS_tree->SetBranchAddress("pi2_BIPCHI2",&pi2_BIPCHI2WS);
  XiccWS_tree->SetBranchAddress("pi3_BIPCHI2",&pi3_BIPCHI2WS);
  XiccWS_tree->SetBranchAddress("Lc3_M",&Lc3_MWS);
  XiccWS_tree->SetBranchAddress("K2_ProbNNk",&K2_ProbNNkWS);
  XiccWS_tree->SetBranchAddress("pi2_ProbNNpi",&pi2_ProbNNpiWS);
  XiccWS_tree->SetBranchAddress("pi3_ProbNNpi",&pi3_ProbNNpiWS);

  for (UInt_t evt = 0; evt < Xicc_tree->GetEntries();evt++) {
    Xicc_tree->GetEntry(evt);

    //if(!(fabs(Xc_M - lcmass) < 20))continue;

    for(int i = 0; i < Xicc_cand; i++){
      bool vetos = !(fabs(Lc3_M[i] - lcmass) < 15);
      bool Sc_mass_cut = fabs(Lcpi1_M[i] - 2454) < 8;
      bool IPcuts = K2_MINIPCHI2[i] > 1.5 && K2_BIPCHI2[i] < 0.6;
      bool PIDcuts = K2_ProbNNk[i] > 0.25;

      if(!(vetos && IPcuts && PIDcuts ))continue;
      if(Sc_mass_cut){
        Xiccp_hist->Fill(LcKpi_M[i]);
      }
      bool Xiccpp_Sc = fabs(Lcpi1_M[i] - 2454) < 8 || fabs(Lcpi2_M[i] - 2454) < 8;
      if(Xiccpp_Sc && LcKpipi_M[i] > 0)
        Xiccpp_hist->Fill(LcKpipi_M[i]);
    }
  }

  for (UInt_t evt = 0; evt < XiccWS_tree->GetEntries();evt++) {
    XiccWS_tree->GetEntry(evt);

    //if(!(fabs(XcWS_M - lcmass) < 20))continue;

    for(int i = 0; i < Xicc_WScand; i++){
      bool vetos = !(fabs(Lc3_MWS[i] - lcmass) < 15);
      bool Sc_mass_cut = fabs(Lcpi1_MWS[i] - 2454) < 8;
      bool IPcuts = K2_MINIPCHI2WS[i] > 1.5 && K2_BIPCHI2WS[i] < 0.6;
      bool PIDcuts = K2_ProbNNkWS[i] > 0.25;

      if(!(vetos && IPcuts && PIDcuts ))continue;
      if(Sc_mass_cut){
        XiccpWC_hist->Fill(LcKpi_MWS[i]);
      }
      bool Xiccpp_Sc = fabs(Lcpi1_MWS[i] - 2454) < 8 || fabs(Lcpi2_MWS[i] - 2454) < 8;
      if(Xiccpp_Sc && LcKpipi_MWS[i] > 0)
        XiccppWC_hist->Fill(LcKpipi_MWS[i]);
    }
  }

  make_overlay_Plot(XiccppWC_hist,Xiccpp_hist,myconfig);
  make_overlay_Plot(XiccpWC_hist,Xiccp_hist,myconfig);


  return;
}

template <class cH>
void fill_profile(vector<cH> multiple_cH_candidate, TProfile *Xc_IP) {
  int n_profilebins = Xc_IP->GetNbinsX();
  double binwidth = Xc_IP->GetXaxis()->GetBinWidth(1);
  vector<int> multiplicity;
  for (int i = 1; i < n_profilebins + 1; i++){
    multiplicity.push_back(0);
    for (int j = 0; j < multiple_cH_candidate.size(); j++){
      if(i < n_profilebins){
        if(multiple_cH_candidate.at(j).ip < (float)(i*binwidth))
          multiplicity[i-1] += 1;
      }
      else multiplicity[i-1] += 1;//get all
    }
  }
  for (int i = 1; i < n_profilebins + 1; i++){
    //cout << "multiplicity below " << i*binwidth << " mm: " << multiplicity.at(i-1) << endl;
    Xc_IP->Fill(Xc_IP->GetXaxis()->GetBinCenter(i),(double)multiplicity.at(i-1));
  }
  return;
}

template<class T>
void make_1D_Plot(T *hist, configuration *myconfig){
  TCanvas *c1 = new TCanvas("c1","C1",10,10,1280,960);
  gPad->SetTopMargin(0.08);
  gPad->SetRightMargin(0.04);
  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.14);
  gPad->SetTickx();
  gPad->SetTicky();
  //gPad->SetLogy();
  //hist->GetYaxis()->SetTitleOffset(0.8);
  hist->SetLineWidth(2);
  if(static_cast<TString>(hist->GetName()).Contains("_IP")){
    hist->SetLineColor(kBlack);
    hist->SetMarkerColor(kBlack);
  }
  else hist->SetLineColor(kBlue+2);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetRangeUser(0,1.2*hist->GetMaximum());
  if(static_cast<TString>(hist->GetName()).Contains("_IP"))hist->Draw("e1p");
  else hist->Draw("hist");

  TPaveText *blank_box = new TPaveText(1.001-gPad->GetRightMargin(),gPad->GetBottomMargin(),1.0,gPad->GetBottomMargin()+0.05,"BRNDC");
  blank_box->SetBorderSize(0);blank_box->SetFillColor(kWhite);blank_box->SetTextAlign(12);blank_box->SetFillStyle(1001);
  blank_box->AddText(" ");
  blank_box->Draw();
  hist->Draw("axis same");

  if(static_cast<TString>(hist->GetName()).Contains("_IP")){
    TLegend *leg;
    leg = new TLegend(0.5,0.12+gPad->GetBottomMargin(),0.97-gPad->GetRightMargin(),0.17+gPad->GetBottomMargin());
    leg->SetBorderSize(0);leg->SetFillColor(kWhite);leg->SetFillStyle(1001);leg->SetTextAlign(12);leg->SetTextSize(0.05);leg->SetTextFont(42);
    temp = "#Xi_{c}^{+}+n negative Tracks";
    if(myconfig->get_particle().Contains("Xic0"))temp = "#Xi_{c}^{0} + n positive Tracks";
    if(myconfig->get_particle().Contains("Omega"))temp = "#Omega_{c}^{0} + n positive Tracks^{+}";
    leg->AddEntry(hist,temp,"lp");
    leg->Draw();
  }

  temp = myconfig->get_dumpdir()+"/Plots/";
  if(!gSystem->OpenDirectory(temp))gSystem->mkdir(temp);

  temp += static_cast<TString>(hist->GetName()) + "_" + myconfig->get_current_cs() + ".pdf";
  if(!myconfig->is_debug()) c1->SaveAs(temp);
  delete c1;
  return;
}

void make_2D_Plot(TH2D *hist, configuration *myconfig){
  TCanvas *c1 = new TCanvas("c1","C1",10,10,1280,960);
  gPad->SetTopMargin(0.08);
  gPad->SetRightMargin(0.14);
  gPad->SetLeftMargin(0.14);
  gPad->SetBottomMargin(0.14);
  gPad->SetTickx();
  gPad->SetTicky();
  //gPad->SetLogy();
  hist->GetYaxis()->SetTitleOffset(0.95);
  hist->Draw("colz");

  TPaveText *blank_box = new TPaveText(1.001-gPad->GetRightMargin(),gPad->GetBottomMargin(),1.0,gPad->GetBottomMargin()+0.05,"BRNDC");
  blank_box->SetBorderSize(0);blank_box->SetFillColor(kWhite);blank_box->SetTextAlign(12);blank_box->SetFillStyle(1001);
  blank_box->AddText(" ");
  blank_box->Draw();

  TPaveText *sbb_box = new TPaveText(0.001+gPad->GetLeftMargin(),1.001-gPad->GetTopMargin(),0.1+gPad->GetLeftMargin(),1.0,"BRNDC");
  sbb_box->SetBorderSize(0);sbb_box->SetFillColor(kWhite);sbb_box->SetTextAlign(12);sbb_box->SetFillStyle(1001);
  sbb_box->AddText(" ");
  sbb_box->Draw();

  hist->Draw("z axis same");

  temp = myconfig->get_dumpdir()+"/Plots/";
  if(!gSystem->OpenDirectory(temp))gSystem->mkdir(temp);

  temp += static_cast<TString>(hist->GetName()) + "_" + myconfig->get_current_cs() + ".pdf";
  if(!myconfig->is_debug()) c1->SaveAs(temp);
  delete c1;
  return;
}

void make_overlay_Plot(TH1D *SS_hist, TH1D *OS_hist, configuration *myconfig){
  TCanvas *c1 = new TCanvas("c1","C1",10,10,1280,960);
  gPad->SetTopMargin(0.08);
  gPad->SetRightMargin(0.04);
  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.14);
  gPad->SetTickx();
  gPad->SetTicky();
  //gPad->SetLogy();
  //hist->GetYaxis()->SetTitleOffset(0.8);

  double ymax = TMath::Max(SS_hist->GetMaximum(),OS_hist->GetMaximum());

  //SS_hist->SetLineWidth(2);
  SS_hist->SetLineColor(kOrange+1);

  //OS_hist->SetLineWidth(2);
  OS_hist->SetLineColor(kBlue+2);

  OS_hist->GetXaxis()->SetLabelSize(0.05);
  OS_hist->GetYaxis()->SetLabelSize(0.05);
  OS_hist->GetXaxis()->SetTitleSize(0.05);
  OS_hist->GetYaxis()->SetTitleSize(0.05);
  OS_hist->GetYaxis()->SetRangeUser(0,1.2*ymax);

  //SS_hist->SetFillColorAlpha(kOrange+1, 0.35);
  //SS_hist->SetFillColor(kOrange+1);
  OS_hist->Draw("hist");
  SS_hist->Draw("histsame");

  TPaveText *blank_box = new TPaveText(1.001-gPad->GetRightMargin(),gPad->GetBottomMargin(),1.0,gPad->GetBottomMargin()+0.05,"BRNDC");
  blank_box->SetBorderSize(0);blank_box->SetFillColor(kWhite);blank_box->SetTextAlign(12);blank_box->SetFillStyle(1001);
  blank_box->AddText(" ");
  blank_box->Draw();
  OS_hist->Draw("axis same");

  if(myconfig->get_particle().Contains("Lc") && static_cast<TString>(OS_hist->GetName()).Contains("Lcpi")){

    c1->Update();
    double Sc_mass = 2454;
    double Sc_2520_mass = 2518;
    double Xic0_mass = 2471;//107.5 : Xic' - Xic mass diffrrence

    TArrow Lcpi_arrows;
    Lcpi_arrows.SetAngle(40);
    Lcpi_arrows.SetLineWidth(2);
    Lcpi_arrows.DrawArrow(Sc_mass,0.85*gPad->GetUymax(),Sc_mass,0.88*gPad->GetUymax(),0.03,"<|");
    Lcpi_arrows.DrawArrow(Sc_2520_mass,0.6*gPad->GetUymax(),Sc_2520_mass,0.7*gPad->GetUymax(),0.03,"<|");
    Lcpi_arrows.DrawArrow(Xic0_mass,0.6*gPad->GetUymax(),Xic0_mass,0.7*gPad->GetUymax(),0.03,"<|");

    TLatex Lcpi_labels;Lcpi_labels.SetTextAlign(12);Lcpi_labels.SetTextFont(42);Lcpi_labels.SetTextSize(0.05);
    Lcpi_labels.DrawLatex(Sc_mass-10,0.93*gPad->GetUymax(),"#Sigma^{0}_{c}(2455)");
    Lcpi_labels.DrawLatex(Sc_2520_mass-10,0.75*gPad->GetUymax(),"#Sigma^{0}_{c}(2520)");
    //Lcpi_labels.SetTextAlign(22);
    Lcpi_labels.DrawLatex(Xic0_mass-10,0.75*gPad->GetUymax(),"#Xi^{0}_{c}");

  }

  if(myconfig->get_particle().Contains("Lc") && static_cast<TString>(OS_hist->GetName()).Contains("Xiccpp_hist")){

    c1->Update();
    double LQCD_Xicc_mass = 3610;

    TArrow Xiccpp_arrows;
    Xiccpp_arrows.SetAngle(40);
    Xiccpp_arrows.SetLineWidth(2);
    Xiccpp_arrows.DrawArrow(LQCD_Xicc_mass,0.3*gPad->GetUymax(),LQCD_Xicc_mass,0.45*gPad->GetUymax(),0.03,"<|");

    TLatex Xiccpp_labels;Xiccpp_labels.SetTextAlign(12);Xiccpp_labels.SetTextFont(42);Xiccpp_labels.SetTextSize(0.05);
    Xiccpp_labels.DrawLatex(LQCD_Xicc_mass-10,0.5*gPad->GetUymax(),"LQCD #Xi_{cc}");
  }

  if(myconfig->get_particle().Contains("Lc") && static_cast<TString>(OS_hist->GetName()).Contains("Xiccp_hist")){

    c1->Update();
    double Xic_2980_mass = 2971;
    double Xic_3080_mass = 3077;
    double LQCD_Xicc_mass = 3610;

    TArrow Xiccp_arrows;
    Xiccp_arrows.SetAngle(40);
    Xiccp_arrows.SetLineWidth(2);
    Xiccp_arrows.DrawArrow(Xic_2980_mass,0.6*gPad->GetUymax(),Xic_2980_mass,0.7*gPad->GetUymax(),0.03,"<|");
    Xiccp_arrows.DrawArrow(Xic_3080_mass,0.85*gPad->GetUymax(),Xic_3080_mass,0.88*gPad->GetUymax(),0.03,"<|");
    Xiccp_arrows.DrawArrow(LQCD_Xicc_mass,0.3*gPad->GetUymax(),LQCD_Xicc_mass,0.45*gPad->GetUymax(),0.03,"<|");

    TLatex Xiccp_labels;Xiccp_labels.SetTextAlign(12);Xiccp_labels.SetTextFont(42);Xiccp_labels.SetTextSize(0.05);
    Xiccp_labels.DrawLatex(Xic_2980_mass-10,0.75*gPad->GetUymax(),"#Xi^{+}_{c}(2980)");
    Xiccp_labels.DrawLatex(Xic_3080_mass-10,0.93*gPad->GetUymax(),"#Xi^{+}_{c}(3080)");
    Xiccp_labels.DrawLatex(LQCD_Xicc_mass-10,0.5*gPad->GetUymax(),"LQCD #Xi_{cc}");
  }

  if(myconfig->get_particle().Contains("Xic") && static_cast<TString>(OS_hist->GetName()).Contains("pi")){

    c1->Update();
    double Xic_2645_mass = 2646;
    double Xic_2790_mass = 2790 - 107.5;//107.5 : Xic' - Xic mass diffrrence
    double LQCD_Xicc_mass = 3610;

    TArrow Xic_arrows;
    Xic_arrows.SetAngle(40);
    Xic_arrows.SetLineWidth(2);
    Xic_arrows.DrawArrow(Xic_2645_mass,0.85*gPad->GetUymax(),Xic_2645_mass,0.88*gPad->GetUymax(),0.03,"<|");
    Xic_arrows.DrawArrow(Xic_2790_mass,0.6*gPad->GetUymax(),Xic_2790_mass,0.7*gPad->GetUymax(),0.03,"<|");
    Xic_arrows.DrawArrow(LQCD_Xicc_mass,0.10*gPad->GetUymax(),LQCD_Xicc_mass,0.2*gPad->GetUymax(),0.03,"<|");

    TLatex Xic_labels;Xic_labels.SetTextAlign(12);Xic_labels.SetTextFont(42);Xic_labels.SetTextSize(0.05);
    Xic_labels.DrawLatex(Xic_2645_mass-10,0.93*gPad->GetUymax(),"#Xi_{c}(2645)");
    Xic_labels.DrawLatex(Xic_2790_mass-10,0.75*gPad->GetUymax(),"#Xi_{c}(2790)#rightarrow#Xi_{c}'#pi");
    Xic_labels.SetTextAlign(22);
    Xic_labels.DrawLatex(LQCD_Xicc_mass,0.25*gPad->GetUymax(),"LQCD #Xi_{cc}");

  }

  if(myconfig->get_particle().CompareTo("Xic") == 0 && static_cast<TString>(OS_hist->GetName()).Contains("K")){

    c1->Update();

    double Marco_1_QVal = 40.0;
    double Marco_2_QVal = 90.0;
    double Marco_3_QVal = 105.0;
    double Marco_4_QVal = 128.0;
    double Marco_5_QVal = 160.0;

    double peak1_mass = Marco_1_QVal + 2469 + 494;
    double peak2_mass = Marco_2_QVal + 2469 + 494;
    double peak3_mass = Marco_3_QVal + 2469 + 494;
    double peak4_mass = Marco_4_QVal + 2469 + 494;
    double peak5_mass = Marco_5_QVal + 2469 + 494;

    TArrow OmegacStar_arrows;
    OmegacStar_arrows.SetAngle(40);
    OmegacStar_arrows.SetLineWidth(2);
    OmegacStar_arrows.DrawArrow(peak1_mass,0.83*gPad->GetUymax(),peak1_mass,0.88*gPad->GetUymax(),0.03,"<|");
    OmegacStar_arrows.DrawArrow(peak2_mass,0.83*gPad->GetUymax(),peak2_mass,0.88*gPad->GetUymax(),0.03,"<|");
    OmegacStar_arrows.DrawArrow(peak3_mass,0.83*gPad->GetUymax(),peak3_mass,0.88*gPad->GetUymax(),0.03,"<|");
    OmegacStar_arrows.DrawArrow(peak4_mass,0.83*gPad->GetUymax(),peak4_mass,0.88*gPad->GetUymax(),0.03,"<|");
    OmegacStar_arrows.DrawArrow(peak5_mass,0.63*gPad->GetUymax(),peak5_mass,0.88*gPad->GetUymax(),0.03,"<|");

    TLatex OmegacStar_labels;OmegacStar_labels.SetTextAlign(22);OmegacStar_labels.SetTextFont(42);OmegacStar_labels.SetTextSize(0.05);
    OmegacStar_labels.DrawLatex(peak1_mass,0.93*gPad->GetUymax(),"1");
    OmegacStar_labels.DrawLatex(peak2_mass,0.93*gPad->GetUymax(),"2");
    OmegacStar_labels.DrawLatex(peak3_mass,0.93*gPad->GetUymax(),"3");
    OmegacStar_labels.DrawLatex(peak4_mass,0.93*gPad->GetUymax(),"4");
    OmegacStar_labels.DrawLatex(peak5_mass,0.93*gPad->GetUymax(),"5");

  }

  TLegend *leg;
  if(static_cast<TString>(SS_hist->GetName()).Contains("SB"))leg = new TLegend(gPad->GetLeftMargin()+0.06,0.82-gPad->GetTopMargin(),0.35-gPad->GetLeftMargin(),0.97-gPad->GetTopMargin());
  else leg = new TLegend(0.6,0.72-gPad->GetTopMargin(),0.97-gPad->GetRightMargin(),0.97-gPad->GetTopMargin());
  leg->SetBorderSize(0);leg->SetFillColor(kWhite);leg->SetFillStyle(1001);leg->SetTextAlign(12);leg->SetTextSize(0.05);leg->SetTextFont(42);
  if(static_cast<TString>(OS_hist->GetName()).Contains("IPfail"))leg->SetHeader("Multiple candidates");

  if(static_cast<TString>(OS_hist->GetName()).Contains("pi")){
    temp = "#Xi_{c}^{+}#pi^{-} + c.c.";
    if(myconfig->get_particle().Contains("Lc"))temp = "#Lambda_{c}^{+}#pi^{-} + c.c.";
    if(myconfig->get_particle().Contains("Xic0"))temp = "#Xi_{c}^{0}#pi^{-} + c.c.";
    if(myconfig->get_particle().Contains("Omega"))temp = "#Omega_{c}^{0}#pi^{-} + c.c.";
  }
  else if(static_cast<TString>(OS_hist->GetName()).Contains("K")){
    temp = "#Xi_{c}^{+}K^{-} + c.c.";
    if(myconfig->get_particle().Contains("Lc"))temp = "#Lambda_{c}^{+}K^{-} + c.c.";
    if(myconfig->get_particle().Contains("Xic0"))temp = "#Xi_{c}^{0}K^{-} + c.c.";
    if(myconfig->get_particle().Contains("Omega"))temp = "#Omega_{c}^{0}K^{-} + c.c.";
  }
  else if(static_cast<TString>(OS_hist->GetName()).Contains("Xiccp_hist")){
    temp = "#Sigma_{c}^{++}K^{-} + c.c.";
  }
  else if(static_cast<TString>(OS_hist->GetName()).Contains("Xiccpp_hist")){
    temp = "#Sigma_{c}^{++}K^{-}#pi^{+} + c.c.";
  }
  else{
    temp = "from #Xi_{c}^{+} signal";
    if(myconfig->get_particle().Contains("Lc"))temp = "from #Lambda_{c}^{+} signal";
    if(myconfig->get_particle().Contains("Xic0"))temp = "from #Xi_{c}^{0} signal";
    if(myconfig->get_particle().Contains("Omega"))temp = "from #Omega_{c}^{0} signal";
  }
  leg->AddEntry(OS_hist,temp,"l");
  if(static_cast<TString>(OS_hist->GetName()).Contains("pi")){
    temp = "#Xi_{c}^{+}#pi^{+} + c.c.";
    if(myconfig->get_particle().Contains("Lc"))temp = "#Lambda_{c}^{+}#pi^{+} + c.c.";
    if(myconfig->get_particle().Contains("Xic0"))temp = "#Xi_{c}^{0}#pi^{+} + c.c.";
    if(myconfig->get_particle().Contains("Omega"))temp = "#Omega_{c}^{0}#pi^{+} + c.c.";
  }
  else if(static_cast<TString>(OS_hist->GetName()).Contains("K")){
    temp = "#Xi_{c}^{+}K^{+} + c.c.";
    if(myconfig->get_particle().Contains("Lc"))temp = "#Lambda_{c}^{+}K^{+} + c.c.";
    if(myconfig->get_particle().Contains("Xic0"))temp = "#Xi_{c}^{0}K^{+} + c.c.";
    if(myconfig->get_particle().Contains("Omega"))temp = "#Omega_{c}^{0}K^{+} + c.c.";
  }
  else if(static_cast<TString>(OS_hist->GetName()).Contains("Xiccp_hist")){
    temp = "#Sigma_{c}^{++}K^{+} + c.c.";
  }
  else if(static_cast<TString>(OS_hist->GetName()).Contains("Xiccpp_hist")){
    temp = "#Sigma_{c}^{++}K^{+}#pi^{+} + c.c.";
  }
  else temp = "from sidebands";
  leg->AddEntry(SS_hist,temp,"l");
  leg->Draw();

  temp = myconfig->get_dumpdir()+"/Plots/";
  if(!gSystem->OpenDirectory(temp))gSystem->mkdir(temp);

  temp += static_cast<TString>(OS_hist->GetName()) + "_" + myconfig->get_current_cs() + ".pdf";
  if(!myconfig->is_debug()) c1->SaveAs(temp);
  delete c1;
  return;
}
