#include "TCut.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"

const TCut trueEleCut = "isTrue == 1";
const TCut fakeEleCut = "isTrue == 0 || isTrue == 3";

const  TCut preselectionCuts = "passConversionVeto && abs(dz)<1";

const float etaMin = -2.5;
const float etaMax = +2.5;
const int nEtaBins = 50;

// We are interested to start with pt from 20 GeV, but we actually start
// from 14 GeV so that there is data for interpolation
const float ptMin = 14; 
float dptMin = 2; // ptMin + dptMin should be 20 (REVISED: not necessarily?)
const float ptMax = 200;
const int nPtBins = 35; // Variable pt bins are actually used

const bool smallEventCount = false;
const int smallMaxEvents = 100000;

//  Files IN 
const TString fileNameS = "~/DYJetsToLL_cutID_tuning_92X_v1.root";
const TString fileNameB = "~/TTJets_cutID_92X_v1.root";
// Tree Name (file IN):
const TString treeName = "ntupler/ElectronTree";
// File with weights OUT:
const TString fileNameWeights = "kinematicWeights.root";

// Forward declarations
TTree *getTree(TString fname, TString tname);

// Main function
void computeKinematicWeights(){

  gStyle->SetPalette(1);

  //
  // Load the signal and background trees
  //
  TTree *treeS = getTree(fileNameS, treeName);
  TTree *treeB = getTree(fileNameB, treeName);

  // 
  // Book weight histograms
  //
  TString hPtEtaSignalName = "hPtEtaSignal";
  TString hPtEtaBackgroundName = "hPtEtaBackground";
  // Set up variable bin limits for pt
  Double_t ptBinLimits[nPtBins+1];
  ptBinLimits[0] = ptMin;
  float dpt =dptMin;
  for(int ipt = 1; ipt<=nPtBins; ipt++){
    float ptVal = ptBinLimits[ipt-1];
    if( ptVal >= ptMin && ptVal<50)
      dpt = dptMin;
    else if( ptVal >= 50 && ptVal<100)
      dpt = 5;
    else if ( ptVal >=100 && ptVal < 150 )
      dpt = 10;
    else
      dpt = 25;
    
    ptBinLimits[ipt] = ptVal + dpt;
  }    
  
  TH2D *hPtEtaSignal = new TH2D(hPtEtaSignalName,"",nPtBins, ptBinLimits, nEtaBins, etaMin, etaMax);
  TH2D *hPtEtaBackground = new TH2D(hPtEtaBackgroundName,"",nPtBins, ptBinLimits, nEtaBins, etaMin, etaMax);

  // Define cuts
  const TCut preselectionS = trueEleCut && preselectionCuts;
  const TCut preselectionB = fakeEleCut && preselectionCuts;
  
  // Fill kinematic histograms
  int maxEventsS = treeS->GetEntries();
  if( smallEventCount )
    maxEventsS = smallMaxEvents;
  TString drawCommandS = TString::Format("etaSC:pt>>%s", hPtEtaSignalName.Data());
  printf("Process the signal tree\n");
  treeS->Draw(drawCommandS, preselectionS, "", maxEventsS);
  treeS->Print();
  printf("%s\n", drawCommandS.Data());
  preselectionS.Print();
  hPtEtaSignal->Print();

  int maxEventsB = treeB->GetEntries();
  if( smallEventCount )
    maxEventsB = smallMaxEvents;
  TString drawCommandB = TString::Format("etaSC:pt>>%s", hPtEtaBackgroundName.Data());
  printf("Process the background tree\n");
  treeB->Draw(drawCommandB, preselectionB, "", maxEventsB);

  //
  // Compute the weights histogram
  //
  // Normalize signal and background distributions to unit area
  TH2D *hPtEtaSignalNorm = (TH2D*)hPtEtaSignal->Clone("hPtEtaSignalNorm");
  hPtEtaSignalNorm->Scale(1.0/hPtEtaSignalNorm->GetSumOfWeights());
  TH2D *hPtEtaBackgroundNorm = (TH2D*)hPtEtaBackground->Clone("hPtEtaBackgroundNorm");
  hPtEtaBackgroundNorm->Scale(1.0/hPtEtaBackgroundNorm->GetSumOfWeights());
  // Find the weights
  TH2D *hKinematicWeights = (TH2D*)hPtEtaBackgroundNorm->Clone("hKinematicWeights");
  hKinematicWeights->Divide(hPtEtaSignalNorm);
  hKinematicWeights->Draw("colz");

  // Write out the result
  TFile *fout = new TFile(fileNameWeights, "recreate");
  fout->cd();
  hPtEtaSignal->Write();
  hPtEtaBackground->Write();
  hKinematicWeights->Write();
  fout->Close();
}

TTree *getTree(TString fname, TString tname){

  TFile *fin = new TFile(fname);
  if( !fin ){
    printf("Failed to open file %s\n", fname.Data());
    assert(0);
  }

  TTree *tree = (TTree*)fin->Get(tname);
  if( !tree ){
    printf("Failed to find tree %s in file %s\n", tname.Data(), fname.Data());
    assert(0);
  }

  return tree;
}

