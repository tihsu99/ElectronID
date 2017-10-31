#include "VarCut.hh"
#include "TFile.h"
#include "TString.h"

void exampleFillCuts(){

  TString fname = "cut_repository/cuts_barrel_test_eff_0999.root";
  TFile *fout = new TFile(fname, "recreate");

  // barrel 99.9% efficiency per cut
  VarCut *varCut = new VarCut();
  varCut->setCutValue("full5x5_sigmaIetaIeta"   ,0.016798); 
  varCut->setCutValue("dEtaIn"                  ,0.021603); 
  varCut->setCutValue("dPhiIn"                  ,0.277912); 
  varCut->setCutValue("hOverE"                  ,0.360530); 
  varCut->setCutValue("relIsoWithDBeta"         ,1.685940); 
  varCut->setCutValue("ooEmooP"                 ,0.329440); 
  varCut->setCutValue("d0"                      ,0.112254); 
  varCut->setCutValue("dz"                      ,0.863300); 
  varCut->setCutValue("expectedMissingInnerHits",2.000050 );
 
  varCut->Write("cuts");

  fout->Close();

  TString fnameEE = "cut_repository/cuts_endcap_test_eff_0999.root";
  TFile *foutEE = new TFile(fnameEE, "recreate");

  // endcap 99.9% efficiency per cut
  VarCut *varCutEE = new VarCut();
  varCutEE->setCutValue("full5x5_sigmaIetaIeta"   ,0.048); 
  varCutEE->setCutValue("dEtaIn"                  ,0.035); 
  varCutEE->setCutValue("dPhiIn"                  ,0.257); 
  varCutEE->setCutValue("hOverE"                  ,0.738); 
  varCutEE->setCutValue("relIsoWithDBeta"         ,1.321); 
  varCutEE->setCutValue("ooEmooP"                 ,0.170); 
  varCutEE->setCutValue("d0"                      ,0.385); 
  varCutEE->setCutValue("dz"                      ,0.957); 
  varCutEE->setCutValue("expectedMissingInnerHits",3.020);
 
  varCutEE->Write("cuts");

  foutEE->Close();

}
