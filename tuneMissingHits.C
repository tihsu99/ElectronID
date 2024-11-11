#include "TArrow.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TColor.h"
#include "TGaxis.h"
#include "TMarker.h"

#include "OptimizationConstants.hh"
#include "VarCut.hh"


const bool smallEventCount = false;

const int nHitsValues = 4;
const int hitsMin = 0;
const int hitsMax = 3;

const TString cutFileNamesBarrel[Opt::nWP] = {
  "cut_repository/cuts_barrel_2019-08-23_WP_Veto.root",
  "cut_repository/cuts_barrel_2019-08-23_WP_Loose.root",
  "cut_repository/cuts_barrel_2019-08-23_WP_Medium.root",
  "cut_repository/cuts_barrel_2019-08-23_WP_Tight.root"
};
  
const TString cutFileNamesEndcap[Opt::nWP] = {
  "cut_repository/cuts_endcap_2019-08-23_WP_Veto.root",
  "cut_repository/cuts_endcap_2019-08-23_WP_Loose.root",
  "cut_repository/cuts_endcap_2019-08-23_WP_Medium.root",
  "cut_repository/cuts_endcap_2019-08-23_WP_Tight.root"
};

const TString cutFileNamesExtend[Opt::nWP] = {
  "cut_repository/cuts_extend_2019-08-23_WP_Veto.root",
  "cut_repository/cuts_extend_2019-08-23_WP_Loose.root",
  "cut_repository/cuts_extend_2019-08-23_WP_Medium.root",
  "cut_repository/cuts_extend_2019-08-23_WP_Tight.root"
};

const int markerColors[nHitsValues] = {kRed, kBlue, kMagenta, kBlack};
const int markerStyles[Opt::nWP]    = {20, 21, 22, 23};

// Forward declarations
TTree *getTreeFromFile(TString fname, TString tname);
void findEfficiencies(int region, TTree *signalTree, TTree *backgroundTree,
		      float &effSignal, float &effBackground, VarCut *cutObject, TCut hitsCut);

// Main function
void tuneMissingHits(){

  // Load trees
  TTree *signalTreeBarrel = getTreeFromFile( Opt::fnameSignalBarrel, Opt::signalTreeName);
  TTree *backgroundTreeBarrel = getTreeFromFile( Opt::fnameBackgroundBarrel, Opt::backgroundTreeName);

  TTree *signalTreeEndcap = getTreeFromFile( Opt::fnameSignalEndcap, Opt::signalTreeName);
  TTree *backgroundTreeEndcap = getTreeFromFile( Opt::fnameBackgroundEndcap, Opt::backgroundTreeName);

  TTree *signalTreeExtend = getTreeFromFile( Opt::fnameSignalExtend, Opt::signalTreeName);
  TTree *backgroundTreeExtend = getTreeFromFile( Opt::fnameBackgroundExtend, Opt::backgroundTreeName);

  // Set up the main canvas
  TCanvas *c1 = new TCanvas("c1","",10,10,600,600);
  TH2F *dummy = new TH2F("dummy","",100, 0.9, 1.0, 100, 0.0, 0.2);
  dummy->SetStats(0);
  dummy->GetXaxis()->SetTitle("signal efficiency");
  dummy->GetYaxis()->SetTitle("background rejection");
  dummy->GetYaxis()->SetTitleOffset(1.4);
  dummy->Draw();
    
  // Loop over working points
  for(int iWP=0; iWP<Opt::nWP; iWP++){

    // Load the working point from a ROOT file
    TFile *cutFile = new TFile(cutFileNamesExtend[iWP]);
    if( !cutFile )
      assert(0);
    VarCut *cutObject = (VarCut*)cutFile->Get("cuts");
    if( !cutObject )
      assert(0);
    
    // Compute the efficiencies
    float effSignal, effBackground;
    for(int ihits = hitsMin; ihits <= hitsMax; ihits++){
    
      TCut hitsCut = TString::Format("expectedMissingInnerHits<=%d",ihits).Data();
      
      int region = 2;
      findEfficiencies(region, signalTreeExtend, backgroundTreeExtend, effSignal, effBackground,
		       cutObject, hitsCut);
      printf("Eff for cut %s with base ID %s      effS= %.5f effB= %.5f\n",
	     hitsCut.GetTitle(), cutFileNamesExtend[iWP].Data(), effSignal, effBackground);  

      // Make a marker and draw it.
      TMarker *marker = new TMarker(effSignal, 1.0-effBackground, 20);
      marker->SetMarkerSize(2);
      marker->SetMarkerColor(markerColors[ihits]);
      marker->SetMarkerStyle(markerStyles[iWP]);
      marker->Draw("same");

    }// end loop over hits cut value

  } // end loop over working points

}

// Compute signal and background efficiencies for given cuts
void findEfficiencies(int region, TTree *signalTree, TTree *backgroundTree,
		      float &effSignal, float &effBackground, VarCut *cutObject, TCut hitsCut){

  TCut etaCut = "";
  if(region == 0){
    etaCut = Opt::etaCutBarrel;
  }
  else if(region == 1){
    etaCut = Opt::etaCutEndcap;
  }
  else etaCut = Opt::etaCutExtend;

  TCut kinematicCuts = Opt::ptCut && etaCut;

  TCut preselectionCuts = kinematicCuts && Opt::otherPreselectionCuts;
  
  TCut signalCuts = preselectionCuts && Opt::trueEleCut;
  TCut backgroundCuts = preselectionCuts && Opt::fakeEleCut;  
 
  TCut selectionCuts = *(cutObject->getCut());

  TH1F *hS_num = new TH1F("hS_num","",100,0.,10000.);
  TH1F *hS_den = new TH1F("hS_den","",100,0.,10000.);
  TH1F *hBG_num = new TH1F("hBG_num","",100,0.,10000.);
  TH1F *hBG_den = new TH1F("hBG_den","",100,0.,10000.);
 

  // draw the histogram
  int maxEntriesS = signalTree->GetEntries();
  if(smallEventCount)
    maxEntriesS = 1000000;
  signalTree->Draw("pt>>hS_num", "genWeight"*(selectionCuts && signalCuts && hitsCut), "goff", maxEntriesS);
  signalTree->Draw("pt>>hS_den", "genWeight"*(selectionCuts && signalCuts), "goff", maxEntriesS);
  
  effSignal = hS_num->GetSumOfWeights()/ hS_den->GetSumOfWeights();
  
  // draw the histogram
  int maxEntriesB = backgroundTree->GetEntries();
  if(smallEventCount)
    maxEntriesB = 1000000;
  backgroundTree->Draw("pt>>hBG_num", "genWeight"*(selectionCuts && backgroundCuts && hitsCut), "goff", maxEntriesB);
  backgroundTree->Draw("pt>>hBG_den", "genWeight"*(selectionCuts && backgroundCuts), "goff", maxEntriesB);
  
  effBackground = hBG_num->GetSumOfWeights()/ hBG_den->GetSumOfWeights();

  delete hS_num ; hS_num = nullptr;
  delete hS_den ; hS_den = nullptr;

  delete hBG_num ; hBG_num = nullptr;
  delete hBG_den ; hBG_den = nullptr;
  
  return;
}

// Get a given tree from a given file name.
TTree *getTreeFromFile(TString fname, TString tname){

  TFile *file = new TFile( fname );
  TTree *tree     = (TTree*) file->Get(tname);
  
  return tree;
}
