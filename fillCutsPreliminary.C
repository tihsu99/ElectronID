#include "VarCut.hh"
#include "TFile.h"
#include "TString.h"

void fillCutsPreliminary(){

  // Preliminary safe cuts from Giovanni Zevi Della Porta

  // 
  // All barrel working points
  //

  TString fname_bveto = "cut_repository/cuts_barrel_preliminary_gzdp_50ns_WP_Veto.root";
  TFile *fout_bveto = new TFile(fname_bveto, "recreate");
  VarCut *varCut_bveto = new VarCut();
  varCut_bveto->setCutValue("full5x5_sigmaIetaIeta"   ,0.012); 
  varCut_bveto->setCutValue("dEtaIn"                  ,0.012); 
  varCut_bveto->setCutValue("dPhiIn"                  ,0.8  ); 
  varCut_bveto->setCutValue("hOverE"                  ,0.15 ); 
  varCut_bveto->setCutValue("relIsoWithDBeta"         ,0.20 ); 
  varCut_bveto->setCutValue("ooEmooP"                 ,1e30 ); 
  varCut_bveto->setCutValue("d0"                      ,0.04 ); 
  varCut_bveto->setCutValue("dz"                      ,0.2  ); 
  varCut_bveto->setCutValue("expectedMissingInnerHits",100  );
  varCut_bveto->Write("cuts");
  fout_bveto->Close();

  TString fname_bloose = "cut_repository/cuts_barrel_preliminary_gzdp_50ns_WP_Loose.root";
  TFile *fout_bloose = new TFile(fname_bloose, "recreate");
  VarCut *varCut_bloose = new VarCut();
  varCut_bloose->setCutValue("full5x5_sigmaIetaIeta"   ,0.012); 
  varCut_bloose->setCutValue("dEtaIn"                  ,0.012); 
  varCut_bloose->setCutValue("dPhiIn"                  ,0.15 ); 
  varCut_bloose->setCutValue("hOverE"                  ,0.12 ); 
  varCut_bloose->setCutValue("relIsoWithDBeta"         ,0.20 ); 
  varCut_bloose->setCutValue("ooEmooP"                 ,0.05 ); 
  varCut_bloose->setCutValue("d0"                      ,0.02 ); 
  varCut_bloose->setCutValue("dz"                      ,0.2  ); 
  varCut_bloose->setCutValue("expectedMissingInnerHits",1.02 );
  varCut_bloose->Write("cuts");
  fout_bloose->Close();

  TString fname_bmedium = "cut_repository/cuts_barrel_preliminary_gzdp_50ns_WP_Medium.root";
  TFile *fout_bmedium = new TFile(fname_bmedium, "recreate");
  VarCut *varCut_bmedium = new VarCut();
  varCut_bmedium->setCutValue("full5x5_sigmaIetaIeta"   ,0.01 ); 
  varCut_bmedium->setCutValue("dEtaIn"                  ,0.009); 
  varCut_bmedium->setCutValue("dPhiIn"                  ,0.06 ); 
  varCut_bmedium->setCutValue("hOverE"                  ,0.12 ); 
  varCut_bmedium->setCutValue("relIsoWithDBeta"         ,0.20 ); 
  varCut_bmedium->setCutValue("ooEmooP"                 ,0.05 ); 
  varCut_bmedium->setCutValue("d0"                      ,0.02 ); 
  varCut_bmedium->setCutValue("dz"                      ,0.1  ); 
  varCut_bmedium->setCutValue("expectedMissingInnerHits",1.02 );
  varCut_bmedium->Write("cuts");
  fout_bmedium->Close();

  TString fname_btight = "cut_repository/cuts_barrel_preliminary_gzdp_50ns_WP_Tight.root";
  TFile *fout_btight = new TFile(fname_btight, "recreate");
  VarCut *varCut_btight = new VarCut();
  varCut_btight->setCutValue("full5x5_sigmaIetaIeta"   ,0.01 ); 
  varCut_btight->setCutValue("dEtaIn"                  ,0.009); 
  varCut_btight->setCutValue("dPhiIn"                  ,0.03 ); 
  varCut_btight->setCutValue("hOverE"                  ,0.12 ); 
  varCut_btight->setCutValue("relIsoWithDBeta"         ,0.15 ); 
  varCut_btight->setCutValue("ooEmooP"                 ,0.05 ); 
  varCut_btight->setCutValue("d0"                      ,0.02 ); 
  varCut_btight->setCutValue("dz"                      ,0.1  ); 
  varCut_btight->setCutValue("expectedMissingInnerHits",0.02 );
  varCut_btight->Write("cuts");
  fout_btight->Close();

  // 
  // All endcap working points
  //

  TString fname_eveto = "cut_repository/cuts_endcap_preliminary_gzdp_50ns_WP_Veto.root";
  TFile *fout_eveto = new TFile(fname_eveto, "recreate");
  VarCut *varCut_eveto = new VarCut();
  varCut_eveto->setCutValue("full5x5_sigmaIetaIeta"   ,0.033); 
  varCut_eveto->setCutValue("dEtaIn"                  ,0.022); 
  varCut_eveto->setCutValue("dPhiIn"                  ,0.7  ); 
  varCut_eveto->setCutValue("hOverE"                  ,1e30 ); 
  varCut_eveto->setCutValue("relIsoWithDBeta"         ,0.20 ); 
  varCut_eveto->setCutValue("ooEmooP"                 ,1e30 ); 
  varCut_eveto->setCutValue("d0"                      ,0.04 ); 
  varCut_eveto->setCutValue("dz"                      ,0.2  ); 
  varCut_eveto->setCutValue("expectedMissingInnerHits",100  );
  varCut_eveto->Write("cuts");
  fout_eveto->Close();

  TString fname_eloose = "cut_repository/cuts_endcap_preliminary_gzdp_50ns_WP_Loose.root";
  TFile *fout_eloose = new TFile(fname_eloose, "recreate");
  VarCut *varCut_eloose = new VarCut();
  varCut_eloose->setCutValue("full5x5_sigmaIetaIeta"   ,0.033); 
  varCut_eloose->setCutValue("dEtaIn"                  ,0.021); 
  varCut_eloose->setCutValue("dPhiIn"                  ,0.10 ); 
  varCut_eloose->setCutValue("hOverE"                  ,0.12 ); 
  varCut_eloose->setCutValue("relIsoWithDBeta"         ,0.20 ); 
  varCut_eloose->setCutValue("ooEmooP"                 ,0.05 ); 
  varCut_eloose->setCutValue("d0"                      ,0.02 ); 
  varCut_eloose->setCutValue("dz"                      ,0.2  ); 
  varCut_eloose->setCutValue("expectedMissingInnerHits",1.02 );
  varCut_eloose->Write("cuts");
  fout_eloose->Close();

  TString fname_emedium = "cut_repository/cuts_endcap_preliminary_gzdp_50ns_WP_Medium.root";
  TFile *fout_emedium = new TFile(fname_emedium, "recreate");
  VarCut *varCut_emedium = new VarCut();
  varCut_emedium->setCutValue("full5x5_sigmaIetaIeta"   ,0.031); 
  varCut_emedium->setCutValue("dEtaIn"                  ,0.019); 
  varCut_emedium->setCutValue("dPhiIn"                  ,0.03 ); 
  varCut_emedium->setCutValue("hOverE"                  ,0.12 ); 
  varCut_emedium->setCutValue("relIsoWithDBeta"         ,0.18 ); 
  varCut_emedium->setCutValue("ooEmooP"                 ,0.05 ); 
  varCut_emedium->setCutValue("d0"                      ,0.02 ); 
  varCut_emedium->setCutValue("dz"                      ,0.1  ); 
  varCut_emedium->setCutValue("expectedMissingInnerHits",1.02 );
  varCut_emedium->Write("cuts");
  fout_emedium->Close();

  TString fname_etight = "cut_repository/cuts_endcap_preliminary_gzdp_50ns_WP_Tight.root";
  TFile *fout_etight = new TFile(fname_etight, "recreate");
  VarCut *varCut_etight = new VarCut();
  varCut_etight->setCutValue("full5x5_sigmaIetaIeta"   ,0.031); 
  varCut_etight->setCutValue("dEtaIn"                  ,0.017); 
  varCut_etight->setCutValue("dPhiIn"                  ,0.02 ); 
  varCut_etight->setCutValue("hOverE"                  ,0.12 ); 
  varCut_etight->setCutValue("relIsoWithDBeta"         ,0.13 ); 
  varCut_etight->setCutValue("ooEmooP"                 ,0.05 ); 
  varCut_etight->setCutValue("d0"                      ,0.02 ); 
  varCut_etight->setCutValue("dz"                      ,0.1  ); 
  varCut_etight->setCutValue("expectedMissingInnerHits",0.02 );
  varCut_etight->Write("cuts");
  fout_etight->Close();


}
