#ifndef OPTIMIZATIONCONSTANTS_HH
#define OPTIMIZATIONCONSTANTS_HH

#include "TFile.h"
#include "TString.h"
#include "TCut.h"

namespace Opt {

  // Global parameters
  const bool useBarrel = false;
  // Events to test and train. 
  // To use ALL available events, split 50/50, set
  // all of the nTrain and nTest to 0
  const int nTrain_Signal     = 300000;
  const int nTrain_Background = 300000;
  const int nTest_Signal      =1000000;
  const int nTest_Background  =1000000;
  // Cut repository directory
  const TString cutRepositoryDir = "./cut_repository";
  
  // TMVA options for MethodCuts
  const TString methodCutsBaseOptions = "!H:!V:FitMethod=GA:EffMethod=EffSel";
  
  //
  // Constants related to working points of interest
  // The sequence needs to be ordered from most loose to most tight
  const int nWP = 4;
  const float effBarrel[nWP]       = {0.95      , 0.90       , 0.80      , 0.70     };
  // Make it possible to have a lower target efficieny for endcap:
  const float effEndcap[nWP]       = {0.95      , 0.90       , 0.80      , 0.70     };
  const float *eff = useBarrel ? effBarrel : effEndcap;
  const TString wpNames[nWP] = {"WP_Veto", "WP_Loose", "WP_Medium", "WP_Tight"};
  enum WorkingPointIndex       {WP_VETO  , WP_LOOSE  , WP_MEDIUM  , WP_TIGHT  };
  
  // 
  // Input ntuples. Note that the signal and background
  // trees can come from the same file and/or can have the same name.
  // Signal and background are separated by cuts.
  //
  TFile *fileSignal = 0;
  TFile *fileBackground = 0;
  const TString fnameSignalBarrel = "./DYJetsToLL_may29_flat_ntuple_full_barrel.root";
  const TString fnameSignalEndcap = "./DYJetsToLL_may29_flat_ntuple_full_endcap.root";
  const TString fnameSignal = useBarrel ? fnameSignalBarrel : fnameSignalEndcap;
  const TString signalTreeName = "electronTree";
  const TString fnameBackgroundBarrel = "./TTJets_may29_flat_ntuple_full_barrel.root";
  const TString fnameBackgroundEndcap = "./TTJets_may29_flat_ntuple_full_endcap.root";
  const TString fnameBackground = useBarrel ? fnameBackgroundBarrel : fnameBackgroundEndcap;
  const TString backgroundTreeName = "electronTree";
  
  //
  // Various preselection and truth matching cuts
  //
  // ( 0: fake, 1: true prompt, 2: true from tau, 3: true non-prompt)
  const TCut trueEleCut = "isTrueEle == 1";
  const TCut fakeEleCut = "isTrueEle == 0 || isTrueEle == 3";
  
  // Kinematics
  const TCut ptCut = "pt>=20";
  const TCut etaCutBarrel = " abs(etaSC) < 1.4442 " ;
  const TCut etaCutEndcap = " abs(etaSC) > 1.566 && abs(etaSC)<2.5" ;
  
  // Anything else
  const TCut otherPreselectionCuts = "passConversionVeto && abs(dz)<1";
}

#endif
