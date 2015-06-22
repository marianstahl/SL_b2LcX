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
void Xichh(){

  TStopwatch *clock = new TStopwatch();
  clock->Start(1);

  gErrorIgnoreLevel = kError;
  TFile *fSLBS = new TFile("/auto/data/mstahl/SLBaryonSpectroscopy/SLBaryonSpectroscopyStrp21.root","read");
  TFile *f1 = new TFile("/auto/data/mstahl/SLBaryonSpectroscopy/Xichh.root","RECREATE");
  fill_friend("Xib02XicMuNu/Xic2pKpi",fSLBS,f1);
  fill_friend("Xib2Xic0MuNu/Xic02pKKpi",fSLBS,f1);
  fSLBS->Print();
  //fill_friend("Lc_WSmu",fSLBS,f1);

  clock->Stop();clock->Print();delete clock;
  return;
}

void fill_friend(TString decay, TFile *origin, TFile *friendfile){

  bool Xic0 = decay.Contains("Xic0");

  const double protonmass = 938.272013; //MeV
  const double pionmass = 139.57018; //MeV
  const double kaonmass = 493.677; //MeV
  const double xicmass = 2467.8; //MeV
  const double xic0mass = 2470.88; //MeV
  //const double muonmass = 105.6583715; //MeV

  double p_PT, p_ETA, p_PHI;
  double K_PT = 0, K_ETA = 0, K_PHI = 0;
  double SSK1_PT = 0, SSK1_ETA = 0, SSK1_PHI = 0;
  double SSK2_PT = 0, SSK2_ETA = 0, SSK2_PHI = 0;
  double pi_PT, pi_ETA, pi_PHI;
  /*double Xb_OWNPV_X, Xb_OWNPV_Y, Xb_OWNPV_Z;
  double Xb_ENDVERTEX_X, Xb_ENDVERTEX_Y, Xb_ENDVERTEX_Z;
  double Xb_PT, Xb_ETA, Xb_PHI, Xb_M;*/
  double Xc_PT, Xc_ETA, Xc_PHI, Xc_M;
  double p_ProbNNp, K_ProbNNk = 0, pi_ProbNNpi;
  double KStar_M2, Lambda1520_M2;
  double p_CosTheta;
  double SSK1_ProbNNk = 0,SSK2_ProbNNk = 0;
  float Added_H_PT[200], Added_H_ETA[200], Added_H_PHI[200], Added_CharmH_M[200];
  int Added_n_Particles, Xc_ID;
  float Added_H_PROBNNPID[200], Added_H_ProbNNpi[200], Added_H_ProbNNk[200];
  float b_IPCHI2[200], PV_IPCHI2[200];

  origin->cd();
  temp = decay+ "/DecayTree";
  TTree *Xc_tree = (TTree*)gDirectory->Get(temp);

  Xc_tree->SetBranchStatus("*",0); //disable all branches
  //now switch on the ones we need (saves a lot of time)
  Xc_tree->SetBranchStatus("Xib_OWNPV_X",1);
  Xc_tree->SetBranchStatus("Xib_OWNPV_Y",1);
  Xc_tree->SetBranchStatus("Xib_OWNPV_Z",1);
  Xc_tree->SetBranchStatus("Xib_ENDVERTEX_X",1);
  Xc_tree->SetBranchStatus("Xib_ENDVERTEX_Y",1);
  Xc_tree->SetBranchStatus("Xib_ENDVERTEX_Z",1);

  Xc_tree->SetBranchStatus("Xic_M",1);
  Xc_tree->SetBranchStatus("Xic_ID",1);
  Xc_tree->SetBranchStatus("Xic_PT",1);
  Xc_tree->SetBranchStatus("Xic_ETA",1);
  Xc_tree->SetBranchStatus("Xic_PHI",1);

  Xc_tree->SetBranchStatus("Added_n_Particles",1);
  Xc_tree->SetBranchStatus("Added_H_PT",1);
  Xc_tree->SetBranchStatus("Added_H_ETA",1);
  Xc_tree->SetBranchStatus("Added_H_PHI",1);
  Xc_tree->SetBranchStatus("Added_H_PROBNNPID",1);
  Xc_tree->SetBranchStatus("Added_H_ProbNNpi",1);
  Xc_tree->SetBranchStatus("Added_H_ProbNNk",1);

  Xc_tree->SetBranchStatus("Added_CharmH_LOGIPCHI2_NEW",1);
  Xc_tree->SetBranchStatus("Added_CharmH_LOGMINIPCHI2",1);

  Xc_tree->SetBranchStatus("p_ProbNNp",1);
  if(Xic0){
    Xc_tree->SetBranchStatus("SSK1_ProbNNk",1);
    Xc_tree->SetBranchStatus("SSK2_ProbNNk",1);
  }
  else
    Xc_tree->SetBranchStatus("K_ProbNNk",1);
  Xc_tree->SetBranchStatus("pi_ProbNNpi",1);
  Xc_tree->SetBranchStatus("p_PT",1);
  Xc_tree->SetBranchStatus("p_ETA",1);
  Xc_tree->SetBranchStatus("p_PHI",1);
  if(Xic0){
    Xc_tree->SetBranchStatus("SSK1_PT",1);
    Xc_tree->SetBranchStatus("SSK1_ETA",1);
    Xc_tree->SetBranchStatus("SSK1_PHI",1);
    Xc_tree->SetBranchStatus("SSK2_PT",1);
    Xc_tree->SetBranchStatus("SSK2_ETA",1);
    Xc_tree->SetBranchStatus("SSK2_PHI",1);
  }
  else{
    Xc_tree->SetBranchStatus("K_PT",1);
    Xc_tree->SetBranchStatus("K_ETA",1);
    Xc_tree->SetBranchStatus("K_PHI",1);
    Xc_tree->SetBranchStatus("Xic_Dalitz_Kminus_piplus_M2",1);
    Xc_tree->SetBranchStatus("Xic_Dalitz_Kminus_pplus_M2",1);
    Xc_tree->SetBranchStatus("p_CosTheta",1);
  }
  Xc_tree->SetBranchStatus("pi_PT",1);
  Xc_tree->SetBranchStatus("pi_ETA",1);
  Xc_tree->SetBranchStatus("pi_PHI",1);

  //set the branch addresses
  /*Xc_tree->SetBranchAddress("Lb_OWNPV_X",&Xb_OWNPV_X);
  Xc_tree->SetBranchAddress("Lb_OWNPV_Y",&Xb_OWNPV_Y);
  Xc_tree->SetBranchAddress("Lb_OWNPV_Z",&Xb_OWNPV_Z);
  Xc_tree->SetBranchAddress("Lb_ENDVERTEX_X",&Xb_ENDVERTEX_X);
  Xc_tree->SetBranchAddress("Lb_ENDVERTEX_Y",&Xb_ENDVERTEX_Y);
  Xc_tree->SetBranchAddress("Lb_ENDVERTEX_Z",&Xb_ENDVERTEX_Z);*/

  Xc_tree->SetBranchAddress("Xic_M",&Xc_M);
  Xc_tree->SetBranchAddress("Xic_ID",&Xc_ID);
  Xc_tree->SetBranchAddress("Xic_PT",&Xc_PT);
  Xc_tree->SetBranchAddress("Xic_ETA",&Xc_ETA);
  Xc_tree->SetBranchAddress("Xic_PHI",&Xc_PHI);

  Xc_tree->SetBranchAddress("Added_n_Particles",&Added_n_Particles);
  Xc_tree->SetBranchAddress("Added_H_PT",&Added_H_PT);
  Xc_tree->SetBranchAddress("Added_H_ETA",&Added_H_ETA);
  Xc_tree->SetBranchAddress("Added_H_PHI",&Added_H_PHI);
  Xc_tree->SetBranchAddress("Added_H_PROBNNPID",&Added_H_PROBNNPID);
  Xc_tree->SetBranchAddress("Added_H_ProbNNpi",&Added_H_ProbNNpi);
  Xc_tree->SetBranchAddress("Added_H_ProbNNk",&Added_H_ProbNNk);

  Xc_tree->SetBranchAddress("Added_CharmH_M",&Added_CharmH_M);
  Xc_tree->SetBranchAddress("Added_CharmH_LOGIPCHI2_NEW",&b_IPCHI2);
  Xc_tree->SetBranchAddress("Added_CharmH_LOGMINIPCHI2",&PV_IPCHI2);

  Xc_tree->SetBranchAddress("p_ProbNNp",&p_ProbNNp);
  if(Xic0){
    Xc_tree->SetBranchAddress("SSK1_ProbNNk",&SSK1_ProbNNk);
    Xc_tree->SetBranchAddress("SSK2_ProbNNk",&SSK2_ProbNNk);
  }
  else
    Xc_tree->SetBranchAddress("K_ProbNNk",&K_ProbNNk);
  Xc_tree->SetBranchAddress("pi_ProbNNpi",&pi_ProbNNpi);

  Xc_tree->SetBranchAddress("p_PT",&p_PT);
  Xc_tree->SetBranchAddress("p_ETA",&p_ETA);
  Xc_tree->SetBranchAddress("p_PHI",&p_PHI);
  if(decay.Contains("Xic0")){
    Xc_tree->SetBranchAddress("SSK1_PT",&SSK1_PT);
    Xc_tree->SetBranchAddress("SSK1_ETA",&SSK1_ETA);
    Xc_tree->SetBranchAddress("SSK1_PHI",&SSK1_PHI);
    Xc_tree->SetBranchAddress("SSK2_PT",&SSK2_PT);
    Xc_tree->SetBranchAddress("SSK2_ETA",&SSK2_ETA);
    Xc_tree->SetBranchAddress("SSK2_PHI",&SSK2_PHI);
  }
  else{
    Xc_tree->SetBranchAddress("K_PT",&K_PT);
    Xc_tree->SetBranchAddress("K_ETA",&K_ETA);
    Xc_tree->SetBranchAddress("K_PHI",&K_PHI);
    Xc_tree->SetBranchAddress("Xic_Dalitz_Kminus_piplus_M2",&KStar_M2);
    Xc_tree->SetBranchAddress("Xic_Dalitz_Kminus_pplus_M2",&Lambda1520_M2);
    Xc_tree->SetBranchAddress("p_CosTheta",&p_CosTheta);
  }
  Xc_tree->SetBranchAddress("pi_PT",&pi_PT);
  Xc_tree->SetBranchAddress("pi_ETA",&pi_ETA);
  Xc_tree->SetBranchAddress("pi_PHI",&pi_PHI);
  //SLBS_tree->AddBranchToCache("*");
  //SLBS_tree->LoadBaskets(1000000000);//Load baskets up to 1 GB to memory

  //double p_beta, K_beta = 0, pi_beta;
  float Xicpipi_M[100],XicKpi_M[100],XicpiK_M[100],Xicpi1_M[100],Xicpi2_M[100],XicK1_M[100],XicK2_M[100];
  float Xichh_PT[100],Xich1_PT[100],Xich2_PT[100];
  float Xichh_ETA[100],Xich1_ETA[100],Xich2_ETA[100];
  float Xichh_PHI[100],Xich1_PHI[100],Xich2_PHI[100];
  float h1_MINIPCHI2[100],h2_MINIPCHI2[100],h1_BIPCHI2[100],h2_BIPCHI2[100];
  float h1_PROBNNPID[100],h2_PROBNNPID[100];
  double p_as_piKpi_M, p_as_KKpi_M, p_as_KKKpi_M;//, pK_as_pipi_M, pK_as_ppi_M, pKpi_as_K_M, pKpi_as_p_M;
  float h1_ProbNNpi[100],h2_ProbNNpi[100],h1_ProbNNk[100],h2_ProbNNk[100];
  int Xicc_cand;

  //backgrounds
  /*float Xic01_M[100],Xic02_M[100],Xic03_M[100];
  float Lc1_M[100],Lc2_M[100],Lc3_M[100],Lc4_M[100],Lc5_M[100];
  float Ds1_M[100],Ds2_M[100],Ds3_M[100],Ds4_M[100],Ds5_M[100];
  float MisID_D1_M[100],MisID_D2_M[100],MisID_D3_M[100],MisID_D4_M[100],MisID_D5_M[100];
  float D1_M[100],D2_M[100],D3_M[100],D4_M[100],D5_M[100],D6_M[100];
  float D01_M[100],D02_M[100],D03_M[100],D04_M[100],D05_M[100],D06_M[100];*/

  friendfile->cd();
  TTree added_Xc_tree(decay,decay);

  added_Xc_tree.Branch("Added_n_Particles", &Added_n_Particles, "Added_n_Particles/I");
  added_Xc_tree.Branch("Xicc_cand", &Xicc_cand, "Xicc_cand/I");
  added_Xc_tree.Branch("Xicpipi_M", &Xicpipi_M, "Xicpipi_M[Xicc_cand]/F");
  added_Xc_tree.Branch("XicKpi_M", &XicKpi_M, "XicKpi_M[Xicc_cand]/F");
  added_Xc_tree.Branch("XicpiK_M", &XicpiK_M, "XicpiK_M[Xicc_cand]/F");
  added_Xc_tree.Branch("Xichh_PT", &Xichh_PT, "Xichh_PT[Xicc_cand]/F");
  added_Xc_tree.Branch("Xichh_ETA", &Xichh_ETA, "Xichh_ETA[Xicc_cand]/F");
  added_Xc_tree.Branch("Xichh_PHI", &Xichh_PHI, "Xichh_PHI[Xicc_cand]/F");  
  added_Xc_tree.Branch("Xicpi1_M", &Xicpi1_M, "Xicpi1_M[Xicc_cand]/F");
  added_Xc_tree.Branch("Xicpi2_M", &Xicpi2_M, "Xicpi2_M[Xicc_cand]/F");
  added_Xc_tree.Branch("XicK1_M", &XicK1_M, "XicK1_M[Xicc_cand]/F");
  added_Xc_tree.Branch("XicK2_M", &XicK2_M, "XicK2_M[Xicc_cand]/F");
  added_Xc_tree.Branch("Xich1_PT", &Xich1_PT, "Xich1_PT[Xicc_cand]/F");
  added_Xc_tree.Branch("Xich1_ETA", &Xich1_ETA, "Xich1_ETA[Xicc_cand]/F");
  added_Xc_tree.Branch("Xich1_PHI", &Xich1_PHI, "Xich1_PHI[Xicc_cand]/F");
  added_Xc_tree.Branch("Xich2_PT", &Xich2_PT, "Xich2_PT[Xicc_cand]/F");
  added_Xc_tree.Branch("Xich2_ETA", &Xich2_ETA, "Xich2_ETA[Xicc_cand]/F");
  added_Xc_tree.Branch("Xich2_PHI", &Xich2_PHI, "Xich2_PHI[Xicc_cand]/F");
  added_Xc_tree.Branch("Xic_ID", &Xc_ID, "Xic_ID/I");
  added_Xc_tree.Branch("Xic_M", &Xc_M, "Xic_M/D");
  added_Xc_tree.Branch("p_ProbNNp", &p_ProbNNp, "p_ProbNNp/D");
  if(Xic0){
    added_Xc_tree.Branch("SSK1_ProbNNk", &SSK1_ProbNNk, "SSK1_ProbNNk/D");
    added_Xc_tree.Branch("SSK2_ProbNNk", &SSK2_ProbNNk, "SSK2_ProbNNk/D");
  }
  else
    added_Xc_tree.Branch("K_ProbNNk", &K_ProbNNk, "K_ProbNNk/D");
  added_Xc_tree.Branch("pi_ProbNNpi", &pi_ProbNNpi, "pi_ProbNNpi/D");
  added_Xc_tree.Branch("h1_PROBNNPID", &h1_PROBNNPID, "h1_PROBNNPID[Xicc_cand]/F");
  added_Xc_tree.Branch("h2_PROBNNPID", &h2_PROBNNPID, "h2_PROBNNPID[Xicc_cand]/F");
  added_Xc_tree.Branch("h1_ProbNNpi", &h1_ProbNNpi, "h1_ProbNNpi[Xicc_cand]/F");
  added_Xc_tree.Branch("h2_ProbNNpi", &h2_ProbNNpi, "h2_ProbNNpi[Xicc_cand]/F");
  added_Xc_tree.Branch("h1_ProbNNk", &h1_ProbNNk, "h1_ProbNNk[Xicc_cand]/F");
  added_Xc_tree.Branch("h2_ProbNNk", &h2_ProbNNk, "h2_ProbNNk[Xicc_cand]/F");
  added_Xc_tree.Branch("h1_MINIPCHI2", &h1_MINIPCHI2, "h1_MINIPCHI2[Xicc_cand]/F");
  added_Xc_tree.Branch("h2_MINIPCHI2", &h2_MINIPCHI2, "h2_MINIPCHI2[Xicc_cand]/F");
  added_Xc_tree.Branch("h1_BIPCHI2", &h1_BIPCHI2, "h1_BIPCHI2[Xicc_cand]/F");
  added_Xc_tree.Branch("h2_BIPCHI2", &h2_BIPCHI2, "h2_BIPCHI2[Xicc_cand]/F");
  /*added_Xc_tree.Branch("p_beta", &p_beta, "p_beta/D");
  added_Xc_tree.Branch("K_beta", &K_beta, "K_beta/D");
  added_Xc_tree.Branch("pi_beta", &pi_beta, "pi_beta/D");*/
  if(Xic0)
    added_Xc_tree.Branch("p_as_KKKpi_M", &p_as_KKKpi_M, "p_as_KKKpi_M/D");
  else{
    added_Xc_tree.Branch("p_as_piKpi_M", &p_as_piKpi_M, "p_as_piKpi_M/D");
    added_Xc_tree.Branch("p_as_KKpi_M", &p_as_KKpi_M, "p_as_KKpi_M/D");
  }
  /*added_Xc_tree.Branch("pK_as_pipi_M", &pK_as_pipi_M, "pK_as_pipi_M/D");
  added_Xc_tree.Branch("pK_as_ppi_M", &pK_as_ppi_M, "pK_as_ppi_M/D");
  added_Xc_tree.Branch("pKpi_as_K_M", &pKpi_as_K_M, "pKpi_as_K_M/D");
  added_Xc_tree.Branch("pKpi_as_p_M", &pKpi_as_p_M, "pKpi_as_p_M/D");
  added_Xc_tree.Branch("Xic01_M", &Xic01_M, "Xic01_M[Xicc_cand]/F");
  added_Xc_tree.Branch("Xic02_M", &Xic02_M, "Xic02_M[Xicc_cand]/F");
  added_Xc_tree.Branch("Xic03_M", &Xic03_M, "Xic03_M[Xicc_cand]/F");
  added_Xc_tree.Branch("Lc1_M", &Lc1_M, "Lc1_M[Xicc_cand]/F");
  added_Xc_tree.Branch("Lc2_M", &Lc2_M, "Lc2_M[Xicc_cand]/F");
  added_Xc_tree.Branch("Lc3_M", &Lc3_M, "Lc3_M[Xicc_cand]/F");
  added_Xc_tree.Branch("Lc4_M", &Lc4_M, "Lc4_M[Xicc_cand]/F");
  added_Xc_tree.Branch("Lc5_M", &Lc5_M, "Lc5_M[Xicc_cand]/F");
  added_Xc_tree.Branch("Ds1_M", &Ds1_M, "Ds1_M[Xicc_cand]/F");
  added_Xc_tree.Branch("Ds2_M", &Ds2_M, "Ds2_M[Xicc_cand]/F");
  added_Xc_tree.Branch("Ds3_M", &Ds3_M, "Ds3_M[Xicc_cand]/F");
  added_Xc_tree.Branch("Ds4_M", &Ds4_M, "Ds4_M[Xicc_cand]/F");
  added_Xc_tree.Branch("Ds5_M", &Ds5_M, "Ds5_M[Xicc_cand]/F");
  added_Xc_tree.Branch("MisID_D1_M", &MisID_D1_M, "MisID_D1_M[Xicc_cand]/F");
  added_Xc_tree.Branch("MisID_D2_M", &MisID_D2_M, "MisID_D2_M[Xicc_cand]/F");
  added_Xc_tree.Branch("MisID_D3_M", &MisID_D3_M, "MisID_D3_M[Xicc_cand]/F");
  added_Xc_tree.Branch("MisID_D4_M", &MisID_D4_M, "MisID_D4_M[Xicc_cand]/F");
  added_Xc_tree.Branch("MisID_D5_M", &MisID_D5_M, "MisID_D5_M[Xicc_cand]/F");
  added_Xc_tree.Branch("D1_M", &D1_M, "D1_M[Xicc_cand]/F");
  added_Xc_tree.Branch("D2_M", &D2_M, "D2_M[Xicc_cand]/F");
  added_Xc_tree.Branch("D3_M", &D3_M, "D3_M[Xicc_cand]/F");
  added_Xc_tree.Branch("D4_M", &D4_M, "D4_M[Xicc_cand]/F");
  added_Xc_tree.Branch("D5_M", &D5_M, "D5_M[Xicc_cand]/F");
  added_Xc_tree.Branch("D6_M", &D6_M, "D6_M[Xicc_cand]/F");
  added_Xc_tree.Branch("D01_M", &D01_M, "D01_M[Xicc_cand]/F");
  added_Xc_tree.Branch("D02_M", &D02_M, "D02_M[Xicc_cand]/F");
  added_Xc_tree.Branch("D03_M", &D03_M, "D03_M[Xicc_cand]/F");
  added_Xc_tree.Branch("D04_M", &D04_M, "D04_M[Xicc_cand]/F");
  added_Xc_tree.Branch("D05_M", &D05_M, "D05_M[Xicc_cand]/F");
  added_Xc_tree.Branch("D06_M", &D06_M, "D06_M[Xicc_cand]/F");*/

  UInt_t Xc_nevents = Xc_tree->GetEntries();
  cout << "Entries in " + decay + " tree: " << Xc_nevents << endl;

  for (UInt_t evt = 0; evt < Xc_nevents;evt++) {
    Xc_tree->GetEntry(evt);

    TLorentzVector proton;
    proton.SetPtEtaPhiM(p_PT,p_ETA,p_PHI,protonmass);
    TLorentzVector kaon;
    kaon.SetPtEtaPhiM(K_PT,K_ETA,K_PHI,kaonmass);
    TLorentzVector SSkaon1;
    SSkaon1.SetPtEtaPhiM(SSK1_PT,SSK1_ETA,SSK1_PHI,kaonmass);
    TLorentzVector SSkaon2;
    SSkaon2.SetPtEtaPhiM(SSK2_PT,SSK2_ETA,SSK2_PHI,kaonmass);
    TLorentzVector pion;
    pion.SetPtEtaPhiM(pi_PT,pi_ETA,pi_PHI,pionmass);

    /*p_beta  = (-proton.P()+kaon.P()+pion.P())/(proton.P()+kaon.P()+pion.P());
    K_beta  = ( proton.P()-kaon.P()+pion.P())/(proton.P()+kaon.P()+pion.P());
    pi_beta = ( proton.P()+kaon.P()-pion.P())/(proton.P()+kaon.P()+pion.P());*/

    TLorentzVector p_as_pi;
    p_as_pi.SetVectM(proton.Vect(),pionmass);
    TLorentzVector p_as_K;
    p_as_K.SetVectM(proton.Vect(),kaonmass);

    /*TLorentzVector K_as_pi;
    K_as_pi.SetVectM(kaon.Vect(),pionmass);
    TLorentzVector K_as_p;
    K_as_p.SetVectM(kaon.Vect(),protonmass);

    TLorentzVector pi_as_K;
    pi_as_K.SetVectM(pion.Vect(),kaonmass);
    TLorentzVector pi_as_p;
    pi_as_p.SetVectM(pion.Vect(),protonmass);*/

    p_as_piKpi_M = (p_as_pi + kaon + pion).M();
    p_as_KKpi_M = (p_as_K + kaon + pion).M();        
    p_as_KKKpi_M = (p_as_K + SSkaon1 + SSkaon2 + pion).M();

    /*pK_as_pipi_M = (proton + K_as_pi + pion).M();
    pK_as_ppi_M = (proton + K_as_p + pion).M();

    pKpi_as_K_M = (proton + kaon + pi_as_K).M();
    pKpi_as_p_M = (proton + kaon + pi_as_p).M();*/

    TLorentzVector Xc;
    Xc.SetPtEtaPhiM(Xc_PT,Xc_ETA,Xc_PHI,Xc_M);
    Xicc_cand = 0;

    //some cuts before we begin
    bool vetos = true;
    bool PID_cuts = true;
    bool Dalitz_region = true;
    bool additional_cuts = true;
    bool Xc_Mass_cut = true;
    if(Xic0){
      vetos = !(1855 < p_as_KKKpi_M && p_as_KKKpi_M < 1875);
      PID_cuts = p_ProbNNp > 0.1 && SSK1_ProbNNk > 0.25 && SSK2_ProbNNk > 0.25 && pi_ProbNNpi > 0.05;
      Xc_Mass_cut = fabs(Xc_M - xic0mass) < 12;
    }
    else{
     vetos = !(1860 < p_as_piKpi_M && p_as_piKpi_M < 1880) && !(2005 < p_as_piKpi_M && p_as_piKpi_M < 2025) && !(1860 < p_as_KKpi_M && p_as_KKpi_M < 1880) && !(1955 < p_as_KKpi_M && p_as_KKpi_M < 1985);
     PID_cuts = p_ProbNNp > 0.28 && K_ProbNNk > 0.15 && pi_ProbNNpi > 0.2;
     Dalitz_region = (842 < sqrt(KStar_M2) && sqrt(KStar_M2) < 942) || (1505 < sqrt(Lambda1520_M2) && sqrt(Lambda1520_M2) < 1535);
     additional_cuts = p_CosTheta > -0.9;
     Xc_Mass_cut = fabs(Xc_M - xicmass) < 15;
    }

    if(!(Xc_Mass_cut && vetos && PID_cuts && Dalitz_region && additional_cuts))continue;

    for(int i = 0; i < Added_n_Particles; i++){//let i be h1, j h2
      if( !(Added_H_PT[i] > 150 && b_IPCHI2[i] < 1.2 && PV_IPCHI2[i] > 0.4 && (fabs(Added_H_PROBNNPID[i]) == 221 || fabs(Added_H_PROBNNPID[i]) == 321) ) )continue;// && b_IPCHI2[i] < 0.95
      for(int j = i+1 ; j < Added_n_Particles; j++){
        if( !(Added_H_PT[j] > 150 && b_IPCHI2[j] < 1.2 && PV_IPCHI2[j] > 0.4 && (fabs(Added_H_PROBNNPID[j]) == 221 || fabs(Added_H_PROBNNPID[j]) == 321) ) )continue;// && b_IPCHI2[j] < 0.95

        h1_PROBNNPID[Xicc_cand] = Added_H_PROBNNPID[i];
        h2_PROBNNPID[Xicc_cand] = Added_H_PROBNNPID[j];

        h1_ProbNNk[Xicc_cand] = Added_H_ProbNNk[i];
        h2_ProbNNk[Xicc_cand] = Added_H_ProbNNk[j];
        h1_ProbNNpi[Xicc_cand] = Added_H_ProbNNpi[i];
        h2_ProbNNpi[Xicc_cand] = Added_H_ProbNNpi[j];

        h1_MINIPCHI2[Xicc_cand] = PV_IPCHI2[i];
        h1_BIPCHI2[Xicc_cand] = b_IPCHI2[i];
        h2_MINIPCHI2[Xicc_cand] = PV_IPCHI2[j];
        h2_BIPCHI2[Xicc_cand] = b_IPCHI2[j];

        TLorentzVector H1pi;
        H1pi.SetPtEtaPhiM(Added_H_PT[i],Added_H_ETA[i],Added_H_PHI[i],pionmass);
        TLorentzVector H1K;
        H1K.SetPtEtaPhiM(Added_H_PT[i],Added_H_ETA[i],Added_H_PHI[i],kaonmass);
        TLorentzVector H2pi;
        H2pi.SetPtEtaPhiM(Added_H_PT[j],Added_H_ETA[j],Added_H_PHI[j],pionmass);
        TLorentzVector H2K;
        H2K.SetPtEtaPhiM(Added_H_PT[j],Added_H_ETA[j],Added_H_PHI[j],kaonmass);

        /*Xic01_M[Xicc_cand] = (float)(proton + kaon + pion + HK).M();
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
        D05_M[Xicc_cand] = (float)(Hpi1 + HK).M();*/

        if(Xic0){
          Xicpi1_M[Xicc_cand] = (float)((Xc + H1pi).M() - Xc_M + xic0mass);
          XicK1_M[Xicc_cand] = (float)((Xc + H1K).M() - Xc_M + xic0mass);
          Xicpi2_M[Xicc_cand] = (float)((Xc + H2pi).M() - Xc_M + xic0mass);
          XicK2_M[Xicc_cand] = (float)((Xc + H2K).M() - Xc_M + xic0mass);
        }
        else{
          Xicpi1_M[Xicc_cand] = (float)((Xc + H1pi).M() - Xc_M + xicmass);
          XicK1_M[Xicc_cand] = (float)((Xc + H1K).M() - Xc_M + xicmass);
          Xicpi2_M[Xicc_cand] = (float)((Xc + H2pi).M() - Xc_M + xicmass);
          XicK2_M[Xicc_cand] = (float)((Xc + H2K).M() - Xc_M + xicmass);
        }
        Xich1_PT[Xicc_cand] = (float)(Xc + H1pi).Pt();
        Xich1_ETA[Xicc_cand] = (float)(Xc + H1pi).Eta();
        Xich1_PHI[Xicc_cand] = (float)(Xc + H1pi).Phi();
        Xich2_PT[Xicc_cand] = (float)(Xc + H2pi).Pt();
        Xich2_ETA[Xicc_cand] = (float)(Xc + H2pi).Eta();
        Xich2_PHI[Xicc_cand] = (float)(Xc + H2pi).Phi();

        Xicpipi_M[Xicc_cand] = (float)(Xc + H1pi + H2pi).M();
        XicKpi_M[Xicc_cand] = (float)(Xc + H1K + H2pi).M();
        XicpiK_M[Xicc_cand] = (float)(Xc + H2K + H1pi).M();

        Xicc_cand++;
        if (Xicc_cand == 99) goto get_the_hell_out;        
      }
    }
get_the_hell_out:
    added_Xc_tree.Fill();
  }
  Xc_tree->SetDirectory(0);
  added_Xc_tree.Write();
  //delete Xc_tree;
  return;
}


