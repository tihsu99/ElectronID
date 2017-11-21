#include "TCut.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include <cassert>
#include <TROOT.h>
#include <stdio.h>
#include <stdlib.h>

const TString tagDir = "2017-11-16";

const TCut trueEleCut       = "isTrue == 1";
const TCut fakeEleCut       = "isTrue == 0 || isTrue == 3";
const TCut preselectionCuts = "passConversionVeto && abs(dz)<1";

const float etaMin = -2.5;
const float etaMax = +2.5;
const int nEtaBins = 50;

// We are interested to start with pt from 20 GeV, but we actually start
// from 14 GeV so that there is data for interpolation
const float ptMin = 14; 
float dptMin = 2;
const float ptMax = 200;
const int nPtBins = 35; // Variable pt bins are actually used

const bool smallEventCount = false;
const int smallMaxEvents = 100000;

const TString getFileName(TString type){
  return "/user/tomc/eleIdTuning/tuples/" + type + "_cutID_tuning_92X_v1.root";
}
// Tree Name (file IN):
const TString treeName = "ntupler/ElectronTree";
// File with weights OUT:
const TString fileNameWeights = "kinematicWeights.root";

// Forward declarations
TTree *getTree(TString fname, TString tname);


TH2D *getPtEtaHist(TString name, TTree* tree, double* ptBinLimits, TCut preselection){
  TH2D *hPtEta = new TH2D(name,"",nPtBins, ptBinLimits, nEtaBins, etaMin, etaMax);
  int maxEvents = tree->GetEntries();
  if(smallEventCount) maxEvents = smallMaxEvents;
  TString drawCommand = TString::Format("etaSC:pt>>%s", name.Data());
  printf(TString::Format("Processing %s\n", name.Data()));
  tree->Draw(drawCommand, preselection, "", maxEvents);
  tree->Print();
  printf("%s\n", drawCommand.Data());
  preselection.Print();
  hPtEta->Print();
  return hPtEta;
}

// Main function
void computeKinematicWeights(){

  gStyle->SetPalette(1);

  //
  // Load the signal and background trees
  //
  TTree *treeS = getTree(getFileName("DYJetsToLL"), treeName);
  TTree *treeB = getTree(getFileName("TTJets"),     treeName);

  // 
  // Book weight histograms
  //
  TString hPtEtaSignalName     = "hPtEtaSignal";
  TString hPtEtaBackgroundName = "hPtEtaBackground";
  // Set up variable bin limits for pt
  Double_t ptBinLimits[nPtBins+1];
  ptBinLimits[0] = ptMin;
  float dpt =dptMin;
  for(int ipt = 1; ipt<=nPtBins; ipt++){
    float ptVal = ptBinLimits[ipt-1];
    if( ptVal >= ptMin && ptVal<50)        dpt = dptMin;
    else if( ptVal >= 50 && ptVal<100)     dpt = 5;
    else if ( ptVal >=100 && ptVal < 150 ) dpt = 10;
    else                                   dpt = 25;
    
    ptBinLimits[ipt] = ptVal + dpt;
  }    
 
  TH2D* hPtEtaSignal     = getPtEtaHist("hPtEtaSignal",     treeS, ptBinLimits, trueEleCut && preselectionCuts); 
  TH2D* hPtEtaBackground = getPtEtaHist("hPtEtaBackground", treeB, ptBinLimits, fakeEleCut && preselectionCuts); 

  //
  // Compute the weights histogram
  //
  // Normalize signal and background distributions to unit area
  TH2D *hPtEtaSignalNorm     = (TH2D*) hPtEtaSignal->Clone("hPtEtaSignalNorm");
  TH2D *hPtEtaBackgroundNorm = (TH2D*) hPtEtaBackground->Clone("hPtEtaBackgroundNorm");
  hPtEtaSignalNorm->Scale(1.0/hPtEtaSignalNorm->GetSumOfWeights());
  hPtEtaBackgroundNorm->Scale(1.0/hPtEtaBackgroundNorm->GetSumOfWeights());

  // Find the weights
  TH2D *hKinematicWeights = (TH2D*) hPtEtaBackgroundNorm->Clone("hKinematicWeights");
  hKinematicWeights->Divide(hPtEtaSignalNorm);
  hKinematicWeights->Draw("colz");

  // Write out the result
  system("mkdir -p " + tagDir);
  TFile *fout = new TFile(tagDir + "/" + fileNameWeights, "recreate");
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

// Compiled
int main(int argc, char *argv[]){
  gROOT->SetBatch();
  computeKinematicWeights();
}
