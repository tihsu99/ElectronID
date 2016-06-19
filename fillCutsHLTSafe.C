#include "VarCut.hh"
#include "TFile.h"
#include "TString.h"

void fillCutsHLTSafe(){

  // Preliminary safe cuts from Giovanni Zevi Della Porta, summer 2016

  // 
  // All barrel working points
  //

  TString fname_bsafe = "cut_repository/cuts_barrel_HLTSafe_20160616_200000.root";
  TFile *fout_bsafe = new TFile(fname_bsafe, "recreate");
  VarCut *varCut_bsafe = new VarCut();
  varCut_bsafe->setCutValue("full5x5_sigmaIetaIeta"   ,0.011); 
  varCut_bsafe->setCutValue("dEtaSeed"                ,0.004); 
  varCut_bsafe->setCutValue("dPhiIn"                  ,0.020); 
  varCut_bsafe->setCutValue("hOverE"                  ,0.06 ); 
  varCut_bsafe->setCutValue("relIsoWithEA"            ,0.10 ); 
  varCut_bsafe->setCutValue("ooEmooP"                 ,0.013); 
  varCut_bsafe->setCutValue("d0"                      ,1e30 ); 
  varCut_bsafe->setCutValue("dz"                      ,1e30 ); 
  // varCut_bsafe->setCutValue("expectedMissingInnerHits",100  );
  varCut_bsafe->Write("cuts");
  fout_bsafe->Close();


  // 
  // All endcap working points
  //

  TString fname_esafe = "cut_repository/cuts_endcap_HLTSafe_20160616_200000.root";
  TFile *fout_esafe = new TFile(fname_esafe, "recreate");
  VarCut *varCut_esafe = new VarCut();
  varCut_esafe->setCutValue("full5x5_sigmaIetaIeta"   ,0.031); 
  varCut_esafe->setCutValue("dEtaSeed"                ,0.007); 
  varCut_esafe->setCutValue("dPhiIn"                  ,0.020); 
  varCut_esafe->setCutValue("hOverE"                  ,0.065); 
  varCut_esafe->setCutValue("relIsoWithEA"            ,0.10 ); 
  varCut_esafe->setCutValue("ooEmooP"                 ,0.013); 
  varCut_esafe->setCutValue("d0"                      ,1e30 ); 
  varCut_esafe->setCutValue("dz"                      ,1e30 ); 
  // varCut_esafe->setCutValue("expectedMissingInnerHits",2  );
  varCut_esafe->Write("cuts");
  fout_esafe->Close();

}
