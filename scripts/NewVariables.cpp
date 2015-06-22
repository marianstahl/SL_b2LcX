/*

    2015-05-28 M. Stahl

    Simple ACLIC compilable script to add new variables to the baryon spectroscopy tuples.

    Use e.g. root -l scripts/NewVariables.cpp++ to run.

*/

#include <iostream>
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
#include "TLorentzVector.h"
#include "TVector3.h"

#include "RooRealVar.h"
#include "RooDataHist.h"

TString temp;

void fill_friend(TString decay, TFile *origin, TFile *friendfile);
void NewVariables(){

  TStopwatch *clock = new TStopwatch();
  clock->Start(1);

  gErrorIgnoreLevel = kError;
  TFile *fSLBS = new TFile("/auto/data/mstahl/SLBaryonSpectroscopy/SLLcStrp21.root","read");
  TFile *f1 = new TFile("/auto/data/mstahl/SLBaryonSpectroscopy/SLLcStrp21_friend.root","RECREATE");
  fill_friend("Lc",fSLBS,f1);
  fSLBS->Print();
  //fill_friend("Lc_WSmu",fSLBS,f1);

  clock->Stop();clock->Print();delete clock;
  return;
}

void fill_friend(TString decay, TFile *origin, TFile *friendfile){

  const double protonmass = 938.272013; //MeV
  const double pionmass = 139.57018; //MeV
  const double kaonmass = 493.677; //MeV
  const double lcmass = 2286.46;
  //const double muonmass = 105.6583715; //MeV

  double p_PT, p_ETA, p_PHI;
  double K_PT, K_ETA, K_PHI;
  double pi_PT, pi_ETA, pi_PHI;
  double Xb_OWNPV_X, Xb_OWNPV_Y, Xb_OWNPV_Z;
  double Xb_ENDVERTEX_X, Xb_ENDVERTEX_Y, Xb_ENDVERTEX_Z;
  double Xb_PT, Xb_ETA, Xb_PHI, Xb_M;
  double Xc_PT, Xc_ETA, Xc_PHI, Xc_M;
  double p_ProbNNp, K_ProbNNk, pi_ProbNNpi;
  float Added_H_PT[200], Added_H_ETA[200], Added_H_PHI[200], Added_CharmH_M[200];
  int Added_n_Particles, Xc_ID;
  float Added_H_PROBNNPID[200], Added_H_ProbNNpi[200], Added_H_ProbNNk[200];
  float b_IPCHI2[200], PV_IPCHI2[200];

  origin->cd();
  temp = decay+ "/DecayTree";
  TTree *Lc_tree = (TTree*)gDirectory->Get(temp);

  Lc_tree->SetBranchStatus("*",0); //disable all branches
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

  Lc_tree->SetBranchStatus("Added_CharmH_LOGIPCHI2_NEW",1);
  Lc_tree->SetBranchStatus("Added_CharmH_LOGMINIPCHI2",1);

  Lc_tree->SetBranchStatus("p_ProbNNp",1);
  Lc_tree->SetBranchStatus("K_ProbNNk",1);
  Lc_tree->SetBranchStatus("pi_ProbNNpi",1);
  Lc_tree->SetBranchStatus("p_PT",1);
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
  Lc_tree->SetBranchAddress("Added_CharmH_LOGIPCHI2_NEW",&b_IPCHI2);
  Lc_tree->SetBranchAddress("Added_CharmH_LOGMINIPCHI2",&PV_IPCHI2);

  Lc_tree->SetBranchAddress("p_ProbNNp",&p_ProbNNp);
  Lc_tree->SetBranchAddress("K_ProbNNk",&K_ProbNNk);
  Lc_tree->SetBranchAddress("pi_ProbNNpi",&pi_ProbNNpi);

  Lc_tree->SetBranchAddress("p_PT",&p_PT);
  Lc_tree->SetBranchAddress("p_ETA",&p_ETA);
  Lc_tree->SetBranchAddress("p_PHI",&p_PHI);
  Lc_tree->SetBranchAddress("K_PT",&K_PT);
  Lc_tree->SetBranchAddress("K_ETA",&K_ETA);
  Lc_tree->SetBranchAddress("K_PHI",&K_PHI);
  Lc_tree->SetBranchAddress("pi_PT",&pi_PT);
  Lc_tree->SetBranchAddress("pi_ETA",&pi_ETA);
  Lc_tree->SetBranchAddress("pi_PHI",&pi_PHI);
  //SLBS_tree->AddBranchToCache("*");
  //SLBS_tree->LoadBaskets(1000000000);//Load baskets up to 1 GB to memory

  double p_beta, K_beta, pi_beta;
  float LcKpipi_M[100],LcKpi_M[100],Lcpi1_M[100],Lcpi2_M[100];
  float LcKpipi_PT[100],LcKpi_PT[100],Lcpi1_PT[100],Lcpi2_PT[100];
  float LcKpipi_ETA[100],LcKpi_ETA[100],Lcpi1_ETA[100],Lcpi2_ETA[100];
  float LcKpipi_PHI[100],LcKpi_PHI[100],Lcpi1_PHI[100],Lcpi2_PHI[100];
  float pi2_MINIPCHI2[100],K2_MINIPCHI2[100],pi3_MINIPCHI2[100],pi2_BIPCHI2[100],K2_BIPCHI2[100],pi3_BIPCHI2[100];
  double p_as_piKpi_M, p_as_KKpi_M, pK_as_pipi_M, pK_as_ppi_M, pKpi_as_K_M, pKpi_as_p_M;
  float K2_ProbNNk[100],pi2_ProbNNpi[100],pi3_ProbNNpi[100];
  int Xicc_cand;

  //backgrounds
  float Xic01_M[100],Xic02_M[100],Xic03_M[100];
  float Lc1_M[100],Lc2_M[100],Lc3_M[100],Lc4_M[100],Lc5_M[100];
  float Ds1_M[100],Ds2_M[100],Ds3_M[100],Ds4_M[100],Ds5_M[100];
  float MisID_D1_M[100],MisID_D2_M[100],MisID_D3_M[100],MisID_D4_M[100],MisID_D5_M[100];
  float D1_M[100],D2_M[100],D3_M[100],D4_M[100],D5_M[100],D6_M[100];
  float D01_M[100],D02_M[100],D03_M[100],D04_M[100],D05_M[100],D06_M[100];

  friendfile->cd();
  TTree added_Lc_tree(decay,decay);

  added_Lc_tree.Branch("Added_n_Particles", &Added_n_Particles, "Added_n_Particles/I");
  added_Lc_tree.Branch("Xicc_cand", &Xicc_cand, "Xicc_cand/I");
  added_Lc_tree.Branch("LcKpipi_M", &LcKpipi_M, "LcKpipi_M[Xicc_cand]/F");
  added_Lc_tree.Branch("LcKpipi_PT", &LcKpipi_PT, "LcKpipi_PT[Xicc_cand]/F");
  added_Lc_tree.Branch("LcKpipi_ETA", &LcKpipi_ETA, "LcKpipi_ETA[Xicc_cand]/F");
  added_Lc_tree.Branch("LcKpipi_PHI", &LcKpipi_PHI, "LcKpipi_PHI[Xicc_cand]/F");
  added_Lc_tree.Branch("LcKpi_M", &LcKpi_M, "LcKpi_M[Xicc_cand]/F");
  added_Lc_tree.Branch("LcKpi_PT", &LcKpi_PT, "LcKpi_PT[Xicc_cand]/F");
  added_Lc_tree.Branch("LcKpi_ETA", &LcKpi_ETA, "LcKpi_ETA[Xicc_cand]/F");
  added_Lc_tree.Branch("LcKpi_PHI", &LcKpi_PHI, "LcKpi_PHI[Xicc_cand]/F");
  added_Lc_tree.Branch("Lcpi1_M", &Lcpi1_M, "Lcpi1_M[Xicc_cand]/F");
  added_Lc_tree.Branch("Lcpi1_PT", &Lcpi1_PT, "Lcpi1_PT[Xicc_cand]/F");
  added_Lc_tree.Branch("Lcpi1_ETA", &Lcpi1_ETA, "Lcpi1_ETA[Xicc_cand]/F");
  added_Lc_tree.Branch("Lcpi1_PHI", &Lcpi1_PHI, "Lcpi1_PHI[Xicc_cand]/F");
  added_Lc_tree.Branch("Lcpi2_M", &Lcpi2_M, "Lcpi2_M[Xicc_cand]/F");
  added_Lc_tree.Branch("Lcpi2_PT", &Lcpi2_PT, "Lcpi2_PT[Xicc_cand]/F");
  added_Lc_tree.Branch("Lcpi2_ETA", &Lcpi2_ETA, "Lcpi2_ETA[Xicc_cand]/F");
  added_Lc_tree.Branch("Lcpi2_PHI", &Lcpi2_PHI, "Lcpi2_PHI[Xicc_cand]/F");
  added_Lc_tree.Branch("Lc_M", &Xc_M, "Lc_M/D");
  added_Lc_tree.Branch("Lb_M", &Xb_M, "Lb_M/D");
  added_Lc_tree.Branch("p_ProbNNp", &p_ProbNNp, "p_ProbNNp/D");
  added_Lc_tree.Branch("K1_ProbNNk", &K_ProbNNk, "K1_ProbNNk/D");
  added_Lc_tree.Branch("pi1_ProbNNpi", &pi_ProbNNpi, "pi1_ProbNNpi/D");
  added_Lc_tree.Branch("K2_ProbNNk", &K2_ProbNNk, "K2_ProbNNk[Xicc_cand]/F");
  added_Lc_tree.Branch("pi2_ProbNNpi", &pi2_ProbNNpi, "pi2_ProbNNpi[Xicc_cand]/F");
  added_Lc_tree.Branch("pi3_ProbNNpi", &pi3_ProbNNpi, "pi3_ProbNNpi[Xicc_cand]/F");
  added_Lc_tree.Branch("K2_MINIPCHI2", &K2_MINIPCHI2, "K2_MINIPCHI2[Xicc_cand]/F");
  added_Lc_tree.Branch("pi2_MINIPCHI2", &pi2_MINIPCHI2, "pi2_MINIPCHI2[Xicc_cand]/F");
  added_Lc_tree.Branch("pi3_MINIPCHI2", &pi3_MINIPCHI2, "pi3_MINIPCHI2[Xicc_cand]/F");
  added_Lc_tree.Branch("K2_BIPCHI2", &K2_BIPCHI2, "K2_BIPCHI2[Xicc_cand]/F");
  added_Lc_tree.Branch("pi2_BIPCHI2", &pi2_BIPCHI2, "pi2_BIPCHI2[Xicc_cand]/F");
  added_Lc_tree.Branch("pi3_BIPCHI2", &pi3_BIPCHI2, "pi3_BIPCHI2[Xicc_cand]/F");
  added_Lc_tree.Branch("p_beta", &p_beta, "p_beta/D");
  added_Lc_tree.Branch("K_beta", &K_beta, "K_beta/D");
  added_Lc_tree.Branch("pi_beta", &pi_beta, "pi_beta/D");
  added_Lc_tree.Branch("p_as_piKpi_M", &p_as_piKpi_M, "p_as_piKpi_M/D");
  added_Lc_tree.Branch("p_as_KKpi_M", &p_as_KKpi_M, "p_as_KKpi_M/D");
  added_Lc_tree.Branch("pK_as_pipi_M", &pK_as_pipi_M, "pK_as_pipi_M/D");
  added_Lc_tree.Branch("pK_as_ppi_M", &pK_as_ppi_M, "pK_as_ppi_M/D");
  added_Lc_tree.Branch("pKpi_as_K_M", &pKpi_as_K_M, "pKpi_as_K_M/D");
  added_Lc_tree.Branch("pKpi_as_p_M", &pKpi_as_p_M, "pKpi_as_p_M/D");
  added_Lc_tree.Branch("Xic01_M", &Xic01_M, "Xic01_M[Xicc_cand]/F");
  added_Lc_tree.Branch("Xic02_M", &Xic02_M, "Xic02_M[Xicc_cand]/F");
  added_Lc_tree.Branch("Xic03_M", &Xic03_M, "Xic03_M[Xicc_cand]/F");
  added_Lc_tree.Branch("Lc1_M", &Lc1_M, "Lc1_M[Xicc_cand]/F");
  added_Lc_tree.Branch("Lc2_M", &Lc2_M, "Lc2_M[Xicc_cand]/F");
  added_Lc_tree.Branch("Lc3_M", &Lc3_M, "Lc3_M[Xicc_cand]/F");
  added_Lc_tree.Branch("Lc4_M", &Lc4_M, "Lc4_M[Xicc_cand]/F");
  added_Lc_tree.Branch("Lc5_M", &Lc5_M, "Lc5_M[Xicc_cand]/F");
  added_Lc_tree.Branch("Ds1_M", &Ds1_M, "Ds1_M[Xicc_cand]/F");
  added_Lc_tree.Branch("Ds2_M", &Ds2_M, "Ds2_M[Xicc_cand]/F");
  added_Lc_tree.Branch("Ds3_M", &Ds3_M, "Ds3_M[Xicc_cand]/F");
  added_Lc_tree.Branch("Ds4_M", &Ds4_M, "Ds4_M[Xicc_cand]/F");
  added_Lc_tree.Branch("Ds5_M", &Ds5_M, "Ds5_M[Xicc_cand]/F");
  added_Lc_tree.Branch("MisID_D1_M", &MisID_D1_M, "MisID_D1_M[Xicc_cand]/F");
  added_Lc_tree.Branch("MisID_D2_M", &MisID_D2_M, "MisID_D2_M[Xicc_cand]/F");
  added_Lc_tree.Branch("MisID_D3_M", &MisID_D3_M, "MisID_D3_M[Xicc_cand]/F");
  added_Lc_tree.Branch("MisID_D4_M", &MisID_D4_M, "MisID_D4_M[Xicc_cand]/F");
  added_Lc_tree.Branch("MisID_D5_M", &MisID_D5_M, "MisID_D5_M[Xicc_cand]/F");
  added_Lc_tree.Branch("D1_M", &D1_M, "D1_M[Xicc_cand]/F");
  added_Lc_tree.Branch("D2_M", &D2_M, "D2_M[Xicc_cand]/F");
  added_Lc_tree.Branch("D3_M", &D3_M, "D3_M[Xicc_cand]/F");
  added_Lc_tree.Branch("D4_M", &D4_M, "D4_M[Xicc_cand]/F");
  added_Lc_tree.Branch("D5_M", &D5_M, "D5_M[Xicc_cand]/F");
  added_Lc_tree.Branch("D6_M", &D6_M, "D6_M[Xicc_cand]/F");
  added_Lc_tree.Branch("D01_M", &D01_M, "D01_M[Xicc_cand]/F");
  added_Lc_tree.Branch("D02_M", &D02_M, "D02_M[Xicc_cand]/F");
  added_Lc_tree.Branch("D03_M", &D03_M, "D03_M[Xicc_cand]/F");
  added_Lc_tree.Branch("D04_M", &D04_M, "D04_M[Xicc_cand]/F");
  added_Lc_tree.Branch("D05_M", &D05_M, "D05_M[Xicc_cand]/F");
  added_Lc_tree.Branch("D06_M", &D06_M, "D06_M[Xicc_cand]/F");

  UInt_t Lc_nevents = Lc_tree->GetEntries();
  cout << "Entries in " + decay + " tree: " << Lc_nevents << endl;

  for (UInt_t evt = 0; evt < Lc_nevents;evt++) {
    Lc_tree->GetEntry(evt);

    TLorentzVector proton;
    proton.SetPtEtaPhiM(p_PT,p_ETA,p_PHI,protonmass);
    TLorentzVector kaon;
    kaon.SetPtEtaPhiM(K_PT,K_ETA,K_PHI,kaonmass);
    TLorentzVector pion;
    pion.SetPtEtaPhiM(pi_PT,pi_ETA,pi_PHI,pionmass);

    p_beta  = (-proton.P()+kaon.P()+pion.P())/(proton.P()+kaon.P()+pion.P());
    K_beta  = ( proton.P()-kaon.P()+pion.P())/(proton.P()+kaon.P()+pion.P());
    pi_beta = ( proton.P()+kaon.P()-pion.P())/(proton.P()+kaon.P()+pion.P());

    TLorentzVector p_as_pi;
    p_as_pi.SetVectM(proton.Vect(),pionmass);
    TLorentzVector p_as_K;
    p_as_K.SetVectM(proton.Vect(),kaonmass);

    TLorentzVector K_as_pi;
    K_as_pi.SetVectM(kaon.Vect(),pionmass);
    TLorentzVector K_as_p;
    K_as_p.SetVectM(kaon.Vect(),protonmass);

    TLorentzVector pi_as_K;
    pi_as_K.SetVectM(pion.Vect(),kaonmass);
    TLorentzVector pi_as_p;
    pi_as_p.SetVectM(pion.Vect(),protonmass);

    p_as_piKpi_M = (p_as_pi + kaon + pion).M();
    p_as_KKpi_M = (p_as_K + kaon + pion).M();

    pK_as_pipi_M = (proton + K_as_pi + pion).M();
    pK_as_ppi_M = (proton + K_as_p + pion).M();

    pKpi_as_K_M = (proton + kaon + pi_as_K).M();
    pKpi_as_p_M = (proton + kaon + pi_as_p).M();

    TLorentzVector Xc;
    Xc.SetPtEtaPhiM(Xc_PT,Xc_ETA,Xc_PHI,Xc_M);
    Xicc_cand = 0;

    //some cuts before we begin
    bool vetos = !(1860 < p_as_piKpi_M && p_as_piKpi_M < 1880) && !(1860 < p_as_KKpi_M && p_as_KKpi_M < 1880) && !(1955 < p_as_KKpi_M && p_as_KKpi_M < 1985);//&& !(2005 < MisID_D2piKpi_M && MisID_D2piKpi_M < 2025)
    bool PID_cuts = p_ProbNNp > 0.05 && K_ProbNNk > 0.05 && pi_ProbNNpi > 0.05;
    //bool CorrM_cut = 4000 < Xb_M && Xb_M < 7200;
    bool Lc_Mass_cut = fabs(Xc_M - lcmass) < 20;

    if(!(Lc_Mass_cut && vetos && PID_cuts))continue;

    for(int i = 0; i < Added_n_Particles; i++){//let i be a K-, j a pi^+ and k a pi^+
      if( !(((Xc_ID > 0 && Added_H_PROBNNPID[i] == -321) || (Xc_ID < 0 && Added_H_PROBNNPID[i] == 321)) && Added_H_PT[i] > 150 && Added_H_ProbNNk[i] > 0.1 && PV_IPCHI2[i] > 0.6))continue;// && b_IPCHI2[i] < 0.95
      for(int j = 0; j < Added_n_Particles; j++){
        if( j== i || !(((Xc_ID > 0 && Added_H_PROBNNPID[j] == 221) || (Xc_ID < 0 && Added_H_PROBNNPID[j] == -221)) && Added_H_PT[j] > 150 && Added_H_ProbNNpi[j] > 0.1  && PV_IPCHI2[j] > 0.6))continue;// && b_IPCHI2[j] < 0.95

        K2_ProbNNk[Xicc_cand] = Added_H_ProbNNk[i];
        pi2_ProbNNpi[Xicc_cand] = Added_H_ProbNNpi[j];

        K2_MINIPCHI2[Xicc_cand] = PV_IPCHI2[i];
        K2_BIPCHI2[Xicc_cand] = b_IPCHI2[i];
        pi2_MINIPCHI2[Xicc_cand] = PV_IPCHI2[j];
        pi2_BIPCHI2[Xicc_cand] = b_IPCHI2[j];

        TLorentzVector HK;
        HK.SetPtEtaPhiM(Added_H_PT[i],Added_H_ETA[i],Added_H_PHI[i],kaonmass);
        TLorentzVector Hpi1;
        Hpi1.SetPtEtaPhiM(Added_H_PT[j],Added_H_ETA[j],Added_H_PHI[j],pionmass);

        Xic01_M[Xicc_cand] = (float)(proton + kaon + pion + HK).M();
        Xic02_M[Xicc_cand] = (float)(proton + kaon + Hpi1 + HK).M();

        Lc1_M[Xicc_cand] = (float)(proton + kaon + Hpi1).M();
        Lc3_M[Xicc_cand] = (float)(proton + HK + pion).M();
        Lc4_M[Xicc_cand] = (float)(proton + Hpi1 + HK).M();

        Ds1_M[Xicc_cand] = (float)(p_as_K + kaon + Hpi1).M();
        Ds3_M[Xicc_cand] = (float)(p_as_K + HK + pion).M();
        Ds4_M[Xicc_cand] = (float)(p_as_K + Hpi1 + HK).M();

        MisID_D1_M[Xicc_cand] = (float)(p_as_pi + kaon + Hpi1).M();
        MisID_D3_M[Xicc_cand] = (float)(p_as_pi + HK + pion).M();
        MisID_D4_M[Xicc_cand] = (float)(p_as_pi + Hpi1 + HK).M();

        D1_M[Xicc_cand] = (float)(kaon + Hpi1 + pion).M();
        D4_M[Xicc_cand] = (float)(Hpi1 + HK + pion).M();

        D01_M[Xicc_cand] = (float)(kaon + pion).M();
        D02_M[Xicc_cand] = (float)(Hpi1 + kaon).M();
        D04_M[Xicc_cand] = (float)(HK + pion).M();
        D05_M[Xicc_cand] = (float)(Hpi1 + HK).M();

        Lcpi1_M[Xicc_cand] = (float)((Xc + Hpi1).M() - Xc_M + lcmass);
        Lcpi1_PT[Xicc_cand] = (float)(Xc + Hpi1).Pt();
        Lcpi1_ETA[Xicc_cand] = (float)(Xc + Hpi1).Eta();
        Lcpi1_PHI[Xicc_cand] = (float)(Xc + Hpi1).Phi();
        LcKpi_M[Xicc_cand] = (float)(Xc + HK + Hpi1).M();//(Added_CharmH_M[j] - Xc_M + lcmass) + 2454.0);
        LcKpi_PT[Xicc_cand] = (float)(Xc + HK + Hpi1).Pt();
        LcKpi_ETA[Xicc_cand] = (float)(Xc + HK + Hpi1).Eta();
        LcKpi_PHI[Xicc_cand] = (float)(Xc + HK + Hpi1).Phi();

        Xicc_cand++;

        pi3_ProbNNpi[Xicc_cand] = 0;
        pi2_MINIPCHI2[Xicc_cand] = 0;
        pi2_BIPCHI2[Xicc_cand] = 0;
        Xic03_M[Xicc_cand] = 0;
        Lc2_M[Xicc_cand] = 0;
        Lc5_M[Xicc_cand] = 0;
        Ds2_M[Xicc_cand] = 0;
        Ds5_M[Xicc_cand] = 0;
        MisID_D2_M[Xicc_cand] = 0;
        MisID_D5_M[Xicc_cand] = 0;
        D2_M[Xicc_cand] = 0;
        D3_M[Xicc_cand] = 0;
        D5_M[Xicc_cand] = 0;
        D6_M[Xicc_cand] = 0;
        D03_M[Xicc_cand] = 0;
        D06_M[Xicc_cand] = 0;
        Lcpi2_M[Xicc_cand] = 0;
        Lcpi2_PT[Xicc_cand] = 0;
        Lcpi2_ETA[Xicc_cand] = 0;
        Lcpi2_PHI[Xicc_cand] = 0;
        LcKpipi_M[Xicc_cand] = 0;
        LcKpipi_PT[Xicc_cand] = 0;
        LcKpipi_ETA[Xicc_cand] = 0;
        LcKpipi_PHI[Xicc_cand] = 0;
        if (Xicc_cand == 99) goto get_the_hell_out;
        for(int k = j+1; k < Added_n_Particles; k++){
          if( k == i || !(((Xc_ID > 0 && Added_H_PROBNNPID[k] == 221) || (Xc_ID < 0 && Added_H_PROBNNPID[k] == -221)) && Added_H_PT[k] > 150 && Added_H_ProbNNpi[k] > 0.1 && PV_IPCHI2[k] > 0.6))continue;//&& b_IPCHI2[k] < 0.95

          pi3_ProbNNpi[Xicc_cand] = Added_H_ProbNNpi[k];
          pi3_MINIPCHI2[Xicc_cand] = PV_IPCHI2[k];
          pi3_BIPCHI2[Xicc_cand] = b_IPCHI2[k];

          TLorentzVector Hpi2;
          Hpi2.SetPtEtaPhiM(Added_H_PT[k],Added_H_ETA[k],Added_H_PHI[k],pionmass);

          Xic03_M[Xicc_cand] = (float)(proton + kaon + Hpi2 + HK).M();
          Lc2_M[Xicc_cand] = (float)(proton + kaon + Hpi2).M();
          Lc5_M[Xicc_cand] = (float)(proton + HK + Hpi2).M();
          Ds2_M[Xicc_cand] = (float)(p_as_K + kaon + Hpi2).M();
          Ds5_M[Xicc_cand] = (float)(p_as_K + HK + Hpi2).M();
          MisID_D2_M[Xicc_cand] = (float)(p_as_pi + kaon + Hpi2).M();
          MisID_D5_M[Xicc_cand] = (float)(p_as_pi + HK + Hpi2).M();
          D2_M[Xicc_cand] = (float)(kaon + pion + Hpi2).M();
          D3_M[Xicc_cand] = (float)(kaon + Hpi1 + Hpi2).M();
          D5_M[Xicc_cand] = (float)(HK + pion + Hpi2).M();
          D6_M[Xicc_cand] = (float)(HK + Hpi1 + Hpi2).M();
          D03_M[Xicc_cand] = (float)(kaon + Hpi2).M();
          D06_M[Xicc_cand] = (float)(HK + Hpi2).M();

          Lcpi2_M[Xicc_cand] = (float)((Xc + Hpi2).M() - Xc_M + lcmass);
          Lcpi2_PT[Xicc_cand] = (float)(Xc + Hpi2).Pt();
          Lcpi2_ETA[Xicc_cand] = (float)(Xc + Hpi2).Eta();
          Lcpi2_PHI[Xicc_cand] = (float)(Xc + Hpi2).Phi();
          LcKpipi_M[Xicc_cand] = (float)(Xc + HK + Hpi1 + Hpi2).M();
          LcKpipi_PT[Xicc_cand] = (float)(Xc + HK + Hpi1 + Hpi2).Pt();
          LcKpipi_ETA[Xicc_cand] = (float)(Xc + HK + Hpi1 + Hpi2).Eta();
          LcKpipi_PHI[Xicc_cand] = (float)(Xc + HK + Hpi1 + Hpi2).Phi();

        }
      }
    }
get_the_hell_out:
    added_Lc_tree.Fill();
  }

  added_Lc_tree.Write();

  friendfile->cd();
  temp = decay+"WS";
  TTree added_LcWS_tree(temp,temp);

  added_LcWS_tree.Branch("Added_n_Particles", &Added_n_Particles, "Added_n_Particles/I");
  added_LcWS_tree.Branch("Xicc_cand", &Xicc_cand, "Xicc_cand/I");
  added_LcWS_tree.Branch("LcKpipi_M", &LcKpipi_M, "LcKpipi_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("LcKpipi_PT", &LcKpipi_PT, "LcKpipi_PT[Xicc_cand]/F");
  added_LcWS_tree.Branch("LcKpipi_ETA", &LcKpipi_ETA, "LcKpipi_ETA[Xicc_cand]/F");
  added_LcWS_tree.Branch("LcKpipi_PHI", &LcKpipi_PHI, "LcKpipi_PHI[Xicc_cand]/F");
  added_LcWS_tree.Branch("LcKpi_M", &LcKpi_M, "LcKpi_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("LcKpi_PT", &LcKpi_PT, "LcKpi_PT[Xicc_cand]/F");
  added_LcWS_tree.Branch("LcKpi_ETA", &LcKpi_ETA, "LcKpi_ETA[Xicc_cand]/F");
  added_LcWS_tree.Branch("LcKpi_PHI", &LcKpi_PHI, "LcKpi_PHI[Xicc_cand]/F");
  added_LcWS_tree.Branch("Lcpi1_M", &Lcpi1_M, "Lcpi1_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("Lcpi1_PT", &Lcpi1_PT, "Lcpi1_PT[Xicc_cand]/F");
  added_LcWS_tree.Branch("Lcpi1_ETA", &Lcpi1_ETA, "Lcpi1_ETA[Xicc_cand]/F");
  added_LcWS_tree.Branch("Lcpi1_PHI", &Lcpi1_PHI, "Lcpi1_PHI[Xicc_cand]/F");
  added_LcWS_tree.Branch("Lcpi2_M", &Lcpi2_M, "Lcpi2_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("Lcpi2_PT", &Lcpi2_PT, "Lcpi2_PT[Xicc_cand]/F");
  added_LcWS_tree.Branch("Lcpi2_ETA", &Lcpi2_ETA, "Lcpi2_ETA[Xicc_cand]/F");
  added_LcWS_tree.Branch("Lcpi2_PHI", &Lcpi2_PHI, "Lcpi2_PHI[Xicc_cand]/F");
  added_LcWS_tree.Branch("Lc_M", &Xc_M, "Lc_M/D");
  added_LcWS_tree.Branch("Lb_M", &Xb_M, "Lb_M/D");
  added_LcWS_tree.Branch("p_ProbNNp", &p_ProbNNp, "p_ProbNNp/D");
  added_LcWS_tree.Branch("K1_ProbNNk", &K_ProbNNk, "K1_ProbNNk/D");
  added_LcWS_tree.Branch("pi1_ProbNNpi", &pi_ProbNNpi, "pi1_ProbNNpi/D");
  added_LcWS_tree.Branch("K2_ProbNNk", &K2_ProbNNk, "K2_ProbNNk[Xicc_cand]/F");
  added_LcWS_tree.Branch("pi2_ProbNNpi", &pi2_ProbNNpi, "pi2_ProbNNpi[Xicc_cand]/F");
  added_LcWS_tree.Branch("pi3_ProbNNpi", &pi3_ProbNNpi, "pi3_ProbNNpi[Xicc_cand]/F");
  added_LcWS_tree.Branch("K2_MINIPCHI2", &K2_MINIPCHI2, "K2_MINIPCHI2[Xicc_cand]/F");
  added_LcWS_tree.Branch("pi2_MINIPCHI2", &pi2_MINIPCHI2, "pi2_MINIPCHI2[Xicc_cand]/F");
  added_LcWS_tree.Branch("pi3_MINIPCHI2", &pi3_MINIPCHI2, "pi3_MINIPCHI2[Xicc_cand]/F");
  added_LcWS_tree.Branch("K2_BIPCHI2", &K2_BIPCHI2, "K2_BIPCHI2[Xicc_cand]/F");
  added_LcWS_tree.Branch("pi2_BIPCHI2", &pi2_BIPCHI2, "pi2_BIPCHI2[Xicc_cand]/F");
  added_LcWS_tree.Branch("pi3_BIPCHI2", &pi3_BIPCHI2, "pi3_BIPCHI2[Xicc_cand]/F");
  added_LcWS_tree.Branch("p_beta", &p_beta, "p_beta/D");
  added_LcWS_tree.Branch("K_beta", &K_beta, "K_beta/D");
  added_LcWS_tree.Branch("pi_beta", &pi_beta, "pi_beta/D");
  added_LcWS_tree.Branch("p_as_piKpi_M", &p_as_piKpi_M, "p_as_piKpi_M/D");
  added_LcWS_tree.Branch("p_as_KKpi_M", &p_as_KKpi_M, "p_as_KKpi_M/D");
  added_LcWS_tree.Branch("pK_as_pipi_M", &pK_as_pipi_M, "pK_as_pipi_M/D");
  added_LcWS_tree.Branch("pK_as_ppi_M", &pK_as_ppi_M, "pK_as_ppi_M/D");
  added_LcWS_tree.Branch("pKpi_as_K_M", &pKpi_as_K_M, "pKpi_as_K_M/D");
  added_LcWS_tree.Branch("pKpi_as_p_M", &pKpi_as_p_M, "pKpi_as_p_M/D");
  added_LcWS_tree.Branch("Xic01_M", &Xic01_M, "Xic01_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("Xic02_M", &Xic02_M, "Xic02_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("Xic03_M", &Xic03_M, "Xic03_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("Lc1_M", &Lc1_M, "Lc1_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("Lc2_M", &Lc2_M, "Lc2_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("Lc3_M", &Lc3_M, "Lc3_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("Lc4_M", &Lc4_M, "Lc4_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("Lc5_M", &Lc5_M, "Lc5_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("Ds1_M", &Ds1_M, "Ds1_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("Ds2_M", &Ds2_M, "Ds2_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("Ds3_M", &Ds3_M, "Ds3_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("Ds4_M", &Ds4_M, "Ds4_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("Ds5_M", &Ds5_M, "Ds5_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("MisID_D1_M", &MisID_D1_M, "MisID_D1_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("MisID_D2_M", &MisID_D2_M, "MisID_D2_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("MisID_D3_M", &MisID_D3_M, "MisID_D3_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("MisID_D4_M", &MisID_D4_M, "MisID_D4_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("MisID_D5_M", &MisID_D5_M, "MisID_D5_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("D1_M", &D1_M, "D1_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("D2_M", &D2_M, "D2_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("D3_M", &D3_M, "D3_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("D4_M", &D4_M, "D4_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("D5_M", &D5_M, "D5_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("D6_M", &D6_M, "D6_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("D01_M", &D01_M, "D01_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("D02_M", &D02_M, "D02_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("D03_M", &D03_M, "D03_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("D04_M", &D04_M, "D04_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("D05_M", &D05_M, "D05_M[Xicc_cand]/F");
  added_LcWS_tree.Branch("D06_M", &D06_M, "D06_M[Xicc_cand]/F");

  for (UInt_t evt = 0; evt < Lc_nevents;evt++) {
    Lc_tree->GetEntry(evt);

    TLorentzVector proton;
    proton.SetPtEtaPhiM(p_PT,p_ETA,p_PHI,protonmass);
    TLorentzVector kaon;
    kaon.SetPtEtaPhiM(K_PT,K_ETA,K_PHI,kaonmass);
    TLorentzVector pion;
    pion.SetPtEtaPhiM(pi_PT,pi_ETA,pi_PHI,pionmass);

    p_beta  = (-proton.P()+kaon.P()+pion.P())/(proton.P()+kaon.P()+pion.P());
    K_beta  = ( proton.P()-kaon.P()+pion.P())/(proton.P()+kaon.P()+pion.P());
    pi_beta = ( proton.P()+kaon.P()-pion.P())/(proton.P()+kaon.P()+pion.P());

    TLorentzVector p_as_pi;
    p_as_pi.SetVectM(proton.Vect(),pionmass);
    TLorentzVector p_as_K;
    p_as_K.SetVectM(proton.Vect(),kaonmass);

    TLorentzVector K_as_pi;
    K_as_pi.SetVectM(kaon.Vect(),pionmass);
    TLorentzVector K_as_p;
    K_as_p.SetVectM(kaon.Vect(),protonmass);

    TLorentzVector pi_as_K;
    pi_as_K.SetVectM(pion.Vect(),kaonmass);
    TLorentzVector pi_as_p;
    pi_as_p.SetVectM(pion.Vect(),protonmass);

    p_as_piKpi_M = (p_as_pi + kaon + pion).M();
    p_as_KKpi_M = (p_as_K + kaon + pion).M();

    pK_as_pipi_M = (proton + K_as_pi + pion).M();
    pK_as_ppi_M = (proton + K_as_p + pion).M();

    pKpi_as_K_M = (proton + kaon + pi_as_K).M();
    pKpi_as_p_M = (proton + kaon + pi_as_p).M();

    TLorentzVector Xc;
    Xc.SetPtEtaPhiM(Xc_PT,Xc_ETA,Xc_PHI,Xc_M);
    Xicc_cand = 0;

    //some cuts before we begin
    bool vetos = !(1860 < p_as_piKpi_M && p_as_piKpi_M < 1880) && !(1860 < p_as_KKpi_M && p_as_KKpi_M < 1880) && !(1955 < p_as_KKpi_M && p_as_KKpi_M < 1985);//&& !(2005 < MisID_D2piKpi_M && MisID_D2piKpi_M < 2025)
    bool PID_cuts = p_ProbNNp > 0.05 && K_ProbNNk > 0.05 && pi_ProbNNpi > 0.05;
    //bool CorrM_cut = 4000 < Xb_M && Xb_M < 7200;
    bool Lc_Mass_cut = fabs(Xc_M - lcmass) < 20;

    if(!(Lc_Mass_cut && vetos && PID_cuts))continue;

    for(int i = 0; i < Added_n_Particles; i++){//let i be a K+, j a pi^+ and k a pi^+
      if( !(((Xc_ID > 0 && Added_H_PROBNNPID[i] == 321) || (Xc_ID < 0 && Added_H_PROBNNPID[i] == -321)) && Added_H_PT[i] > 150 && Added_H_ProbNNk[i] > 0.1 && PV_IPCHI2[i] > 0.6))continue;// && b_IPCHI2[i] < 0.95
      for(int j = 0; j < Added_n_Particles; j++){
        if( j == i || !(((Xc_ID > 0 && Added_H_PROBNNPID[j] == 221) || (Xc_ID < 0 && Added_H_PROBNNPID[j] == -221)) && Added_H_PT[j] > 150 && Added_H_ProbNNpi[j] > 0.1  && PV_IPCHI2[j] > 0.6))continue;// && b_IPCHI2[j] < 0.95

        K2_ProbNNk[Xicc_cand] = Added_H_ProbNNk[i];
        pi2_ProbNNpi[Xicc_cand] = Added_H_ProbNNpi[j];

        K2_MINIPCHI2[Xicc_cand] = PV_IPCHI2[i];
        K2_BIPCHI2[Xicc_cand] = b_IPCHI2[i];
        pi2_MINIPCHI2[Xicc_cand] = PV_IPCHI2[j];
        pi2_BIPCHI2[Xicc_cand] = b_IPCHI2[j];

        TLorentzVector HK;
        HK.SetPtEtaPhiM(Added_H_PT[i],Added_H_ETA[i],Added_H_PHI[i],kaonmass);
        TLorentzVector Hpi1;
        Hpi1.SetPtEtaPhiM(Added_H_PT[j],Added_H_ETA[j],Added_H_PHI[j],pionmass);

        Xic01_M[Xicc_cand] = (float)(proton + kaon + pion + HK).M();
        Xic02_M[Xicc_cand] = (float)(proton + kaon + Hpi1 + HK).M();

        Lc1_M[Xicc_cand] = (float)(proton + kaon + Hpi1).M();
        Lc3_M[Xicc_cand] = (float)(proton + HK + pion).M();
        Lc4_M[Xicc_cand] = (float)(proton + Hpi1 + HK).M();

        Ds1_M[Xicc_cand] = (float)(p_as_K + kaon + Hpi1).M();
        Ds3_M[Xicc_cand] = (float)(p_as_K + HK + pion).M();
        Ds4_M[Xicc_cand] = (float)(p_as_K + Hpi1 + HK).M();

        MisID_D1_M[Xicc_cand] = (float)(p_as_pi + kaon + Hpi1).M();
        MisID_D3_M[Xicc_cand] = (float)(p_as_pi + HK + pion).M();
        MisID_D4_M[Xicc_cand] = (float)(p_as_pi + Hpi1 + HK).M();

        D1_M[Xicc_cand] = (float)(kaon + Hpi1 + pion).M();
        D4_M[Xicc_cand] = (float)(Hpi1 + HK + pion).M();

        D01_M[Xicc_cand] = (float)(kaon + pion).M();
        D02_M[Xicc_cand] = (float)(Hpi1 + kaon).M();
        D04_M[Xicc_cand] = (float)(HK + pion).M();
        D05_M[Xicc_cand] = (float)(Hpi1 + HK).M();

        Lcpi1_M[Xicc_cand] = (float)((Xc + Hpi1).M() - Xc_M + lcmass);
        Lcpi1_PT[Xicc_cand] = (float)(Xc + Hpi1).Pt();
        Lcpi1_ETA[Xicc_cand] = (float)(Xc + Hpi1).Eta();
        Lcpi1_PHI[Xicc_cand] = (float)(Xc + Hpi1).Phi();
        LcKpi_M[Xicc_cand] = (float)(Xc + HK + Hpi1).M();//(Added_CharmH_M[j] - Xc_M + lcmass) + 2454.0);
        LcKpi_PT[Xicc_cand] = (float)(Xc + HK + Hpi1).Pt();
        LcKpi_ETA[Xicc_cand] = (float)(Xc + HK + Hpi1).Eta();
        LcKpi_PHI[Xicc_cand] = (float)(Xc + HK + Hpi1).Phi();

        Xicc_cand++;

        pi3_ProbNNpi[Xicc_cand] = 0;
        pi2_MINIPCHI2[Xicc_cand] = 0;
        pi2_BIPCHI2[Xicc_cand] = 0;
        Xic03_M[Xicc_cand] = 0;
        Lc2_M[Xicc_cand] = 0;
        Lc5_M[Xicc_cand] = 0;
        Ds2_M[Xicc_cand] = 0;
        Ds5_M[Xicc_cand] = 0;
        MisID_D2_M[Xicc_cand] = 0;
        MisID_D5_M[Xicc_cand] = 0;
        D2_M[Xicc_cand] = 0;
        D3_M[Xicc_cand] = 0;
        D5_M[Xicc_cand] = 0;
        D6_M[Xicc_cand] = 0;
        D03_M[Xicc_cand] = 0;
        D06_M[Xicc_cand] = 0;
        Lcpi2_M[Xicc_cand] = 0;
        Lcpi2_PT[Xicc_cand] = 0;
        Lcpi2_ETA[Xicc_cand] = 0;
        Lcpi2_PHI[Xicc_cand] = 0;
        LcKpipi_M[Xicc_cand] = 0;
        LcKpipi_PT[Xicc_cand] = 0;
        LcKpipi_ETA[Xicc_cand] = 0;
        LcKpipi_PHI[Xicc_cand] = 0;
        if (Xicc_cand == 99) goto get_the_hell_out_WS;
        for(int k = j+1; k < Added_n_Particles; k++){
          if( k == i || !(((Xc_ID > 0 && Added_H_PROBNNPID[k] == 221) || (Xc_ID < 0 && Added_H_PROBNNPID[k] == -221)) && Added_H_PT[k] > 150 && Added_H_ProbNNpi[k] > 0.1 && PV_IPCHI2[k] > 0.6))continue;//&& b_IPCHI2[k] < 0.95

          pi3_ProbNNpi[Xicc_cand] = Added_H_ProbNNpi[k];
          pi3_MINIPCHI2[Xicc_cand] = PV_IPCHI2[k];
          pi3_BIPCHI2[Xicc_cand] = b_IPCHI2[k];

          TLorentzVector Hpi2;
          Hpi2.SetPtEtaPhiM(Added_H_PT[k],Added_H_ETA[k],Added_H_PHI[k],pionmass);

          Xic03_M[Xicc_cand] = (float)(proton + kaon + Hpi2 + HK).M();
          Lc2_M[Xicc_cand] = (float)(proton + kaon + Hpi2).M();
          Lc5_M[Xicc_cand] = (float)(proton + HK + Hpi2).M();
          Ds2_M[Xicc_cand] = (float)(p_as_K + kaon + Hpi2).M();
          Ds5_M[Xicc_cand] = (float)(p_as_K + HK + Hpi2).M();
          MisID_D2_M[Xicc_cand] = (float)(p_as_pi + kaon + Hpi2).M();
          MisID_D5_M[Xicc_cand] = (float)(p_as_pi + HK + Hpi2).M();
          D2_M[Xicc_cand] = (float)(kaon + pion + Hpi2).M();
          D3_M[Xicc_cand] = (float)(kaon + Hpi1 + Hpi2).M();
          D5_M[Xicc_cand] = (float)(HK + pion + Hpi2).M();
          D6_M[Xicc_cand] = (float)(HK + Hpi1 + Hpi2).M();
          D03_M[Xicc_cand] = (float)(kaon + Hpi2).M();
          D06_M[Xicc_cand] = (float)(HK + Hpi2).M();

          Lcpi2_M[Xicc_cand] = (float)((Xc + Hpi2).M() - Xc_M + lcmass);
          Lcpi2_PT[Xicc_cand] = (float)(Xc + Hpi2).Pt();
          Lcpi2_ETA[Xicc_cand] = (float)(Xc + Hpi2).Eta();
          Lcpi2_PHI[Xicc_cand] = (float)(Xc + Hpi2).Phi();
          LcKpipi_M[Xicc_cand] = (float)(Xc + HK + Hpi1 + Hpi2).M();
          LcKpipi_PT[Xicc_cand] = (float)(Xc + HK + Hpi1 + Hpi2).Pt();
          LcKpipi_ETA[Xicc_cand] = (float)(Xc + HK + Hpi1 + Hpi2).Eta();
          LcKpipi_PHI[Xicc_cand] = (float)(Xc + HK + Hpi1 + Hpi2).Phi();

        }
      }
    }
get_the_hell_out_WS:
    added_LcWS_tree.Fill();
  }
  added_LcWS_tree.Write();
  Lc_tree->SetDirectory(0);
  //delete Lc_tree;
  return;
}
