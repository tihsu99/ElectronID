#include <iostream>
#include <fstream>
#include "TSystem.h"
#include "TStyle.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TFile.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TTree.h"
#include "TList.h"
#include "TString.h"
#include "TLatex.h"
#include "TLorentzVector.h"
#include "TLegend.h"
#include "TInterpreter.h"
#include "TCut.h"
#include "TBenchmark.h"
#include <signal.h>
#include "TMath.h"

const TString fname = "DYJetsToLL_HLT_study_500K.root";

const float boundaryEBEE = 1.479;
const float maxEta       = 2.5;

const float ptMin = 20;
const float ptMax = 200;

const bool smallEventCount = false;

const int nIsoBins = 1000;
const float isoMin = 0;
const float isoMax = 1;

// Offline effective areas
const int nEtaBins = 7;
const float etaBinLimits[nEtaBins+1] = {
  0.0, 1.0, 1.479, 2.0, 2.2, 2.3, 2.4, 2.5
};
const float effectiveAreaValues[nEtaBins] = { 
  0.1703, 0.1715, 0.1213, 0.1230, 0.1635, 0.1937, 0.2393
};
// HLT effective areas
const int nHltEtaBins = 2;
const float hltEtaBinLimits[nHltEtaBins+1] = {0.0, 1.479, 2.5};
const float hltEAEcal[nHltEtaBins] = {0.165, 0.132};
const float hltEAHcal[nHltEtaBins] = {0.060, 0.131};
// HLT cuts: WP Loose
const float hltEcalIsoCut[nHltEtaBins] = {0.145, 0.135}; 
const float hltHcalIsoCut[nHltEtaBins] = {0.150, 0.130}; 
const float hltTrkIsoCut [nHltEtaBins] = {0.080, 0.080}; 

// Forward declarations
int findOfflineEtaBin(float eta);
int findHltEtaBin(float eta);
TGraph *runUCCOM(TH1F *hnum, TH1F *hden);

// Main function
void computeHLTBounds(bool doBarrel = true){

  //  Generate  a dictionary, so that CINT will have all the information 
  // it needs about type or variable at anytime.
  gROOT->ProcessLine("#include <vector>"); 

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  // HLT cuts
  int hltEtaIndex = 0;
  if( !doBarrel ) 
    hltEtaIndex = 1;

  // Get the tree
  TFile *fin = new TFile(fname);
  TTree *tree = (TTree*)fin->Get("ntupler/ElectronTree");
  
  TH2F *histEcalIso = new TH2F("histEcalIso","",nIsoBins,isoMin,isoMax,
			       nIsoBins,isoMin,isoMax);
  TH1F *histEcalNum = new TH1F("histEcalNum","",nIsoBins,isoMin,isoMax);
  TH1F *histEcalDen = new TH1F("histEcalDen","",nIsoBins,isoMin,isoMax);
  //
  TH1F *histHcalNum = new TH1F("histHcalNum","",nIsoBins,isoMin,isoMax);
  TH1F *histHcalDen = new TH1F("histHcalDen","",nIsoBins,isoMin,isoMax);
  //
  TH1F *histTrkNum = new TH1F("histTrkNum","",nIsoBins,isoMin,isoMax);
  TH1F *histTrkDen = new TH1F("histTrkDen","",nIsoBins,isoMin,isoMax);
  //
  TH1F *histFullNum = new TH1F("histFullNum","",nIsoBins,isoMin,isoMax);
  TH1F *histFullDen = new TH1F("histFullDen","",nIsoBins,isoMin,isoMax);

  // Declare variables to retrieve
  int nEle;
  float rho;
  float rhoCalo;
  std::vector<float> *elePt              = 0;
  std::vector<float> *eleEtaSC           = 0;
  std::vector<int>   *isTrue             = 0;
  //
  std::vector<float> *isoChargedHadrons  = 0;
  std::vector<float> *isoNeutralHadrons  = 0;
  std::vector<float> *isoPhotons         = 0;
  //
  std::vector<float> *isoPFClusterEcal   = 0;
  std::vector<float> *isoPFClusterHcal   = 0;
  std::vector<float> *isoTrk             = 0;

  // Declare branches
  TBranch *b_nEle               = 0;
  TBranch *b_rho                = 0;
  TBranch *b_rhoCalo            = 0;
  TBranch *b_elePt              = 0;
  TBranch *b_eleEtaSC           = 0;
  TBranch *b_isTrue             = 0;
	                       
  TBranch *b_isoChargedHadrons  = 0;
  TBranch *b_isoNeutralHadrons  = 0;
  TBranch *b_isoPhotons         = 0;
	                       
  TBranch *b_isoPFClusterEcal   = 0;
  TBranch *b_isoPFClusterHcal   = 0;
  TBranch *b_isoTrk             = 0;

  // Connect variables and branches to the tree with the data
  tree->SetBranchAddress("nEle", &nEle, &b_nEle);
  tree->SetBranchAddress("rho", &rho, &b_rho);
  tree->SetBranchAddress("rhoCalo", &rhoCalo, &b_rhoCalo);

  tree->SetBranchAddress("pt"    , &elePt   , &b_elePt);
  tree->SetBranchAddress("etaSC" , &eleEtaSC, &b_eleEtaSC);
  tree->SetBranchAddress("isTrue", &isTrue  , &b_isTrue);

  tree->SetBranchAddress("isoChargedHadrons", &isoChargedHadrons, &b_isoChargedHadrons);
  tree->SetBranchAddress("isoNeutralHadrons", &isoNeutralHadrons, &b_isoNeutralHadrons);
  tree->SetBranchAddress("isoPhotons"       , &isoPhotons       , &b_isoPhotons);

  tree->SetBranchAddress("isoPFClusterEcal", &isoPFClusterEcal, &b_isoPFClusterEcal);
  tree->SetBranchAddress("isoPFClusterHcal", &isoPFClusterHcal, &b_isoPFClusterHcal);
  tree->SetBranchAddress("isoTrk"          , &isoTrk          , &b_isoTrk);

  UInt_t maxEvents = tree->GetEntries();
  if( smallEventCount )
    maxEvents = 10000;

  printf("\nStart drawing into hists with total # of events = %u\n", maxEvents );

  //
  // Loop over events
  //

  for(UInt_t ievent = 0; ievent < maxEvents; ievent++){

    Long64_t tentry = tree->LoadTree(ievent);
    // Load the value of the number of the electrons in the event    
    b_nEle->GetEntry(tentry);
    b_rho->GetEntry(tentry);
    b_rhoCalo->GetEntry(tentry);
    
    // Get data for all electrons in this event, only vars of interest
    b_elePt->GetEntry(tentry);
    b_eleEtaSC->GetEntry(tentry);
    b_isTrue->GetEntry(tentry);
    //
    b_isoChargedHadrons->GetEntry(tentry);
    b_isoNeutralHadrons->GetEntry(tentry);
    b_isoPhotons       ->GetEntry(tentry);
    //
    b_isoPFClusterEcal ->GetEntry(tentry);
    b_isoPFClusterHcal ->GetEntry(tentry);
    b_isoTrk           ->GetEntry(tentry);

    for(int iele = 0; iele < nEle; iele++){
     
      // Preselection
      float pt = elePt->at(iele);
      if( !(pt > ptMin && pt < ptMax ) ) continue;
      float abseta = fabs(eleEtaSC->at(iele));
      if( !( (doBarrel && abseta < boundaryEBEE)
	     || ( !doBarrel && abseta >= boundaryEBEE && abseta <= maxEta) )) continue;
      if( isTrue->at(iele) != 1 ) continue;

      // Compute offline corrected isolation
      int offlineEtaBin = findOfflineEtaBin( eleEtaSC->at(iele) );
      float relCombIsoWithEA = 
	(isoChargedHadrons->at(iele) 
	 + max( (float)0.0, isoNeutralHadrons->at(iele) + isoPhotons->at(iele)
		- rho* effectiveAreaValues[offlineEtaBin]) )
	/ elePt->at(iele);

      // Compute HLT corrected isolations
      int hltEtaBin = findHltEtaBin( eleEtaSC->at(iele) );
      float relEcalIso = max((float)0.0, isoPFClusterEcal->at(iele) - rhoCalo*hltEAEcal[hltEtaBin]) 
	/ elePt->at(iele);
      float relHcalIso = max((float)0.0, isoPFClusterHcal->at(iele) - rhoCalo*hltEAHcal[hltEtaBin]) 
	/ elePt->at(iele);
      float relTrkIso = isoTrk->at(iele) / elePt->at(iele);

      // if( relEcalIso + relHcalIso + relTrkIso < 1e-10 ) continue;

      histEcalIso->Fill( relCombIsoWithEA, relEcalIso);

      histEcalDen->Fill( relCombIsoWithEA );
      if( !(relEcalIso < hltEcalIsoCut[hltEtaIndex] ) )
	histEcalNum->Fill( relCombIsoWithEA );
      
      histHcalDen->Fill( relCombIsoWithEA );
      if( !(relHcalIso < hltHcalIsoCut[hltEtaIndex] ) )
	histHcalNum->Fill( relCombIsoWithEA );
      
      histTrkDen->Fill( relCombIsoWithEA );
      if( !(relTrkIso < hltTrkIsoCut[hltEtaIndex] ) )
	histTrkNum->Fill( relCombIsoWithEA );
      
      histFullDen->Fill( relCombIsoWithEA );
      if( !(relEcalIso < hltEcalIsoCut[hltEtaIndex] 
	    && relHcalIso < hltHcalIsoCut[hltEtaIndex]
	    && relTrkIso < hltTrkIsoCut[hltEtaIndex] ) )
	histFullNum->Fill( relCombIsoWithEA );
      
    } // end loop over electrons

  } // end loop over events

  TCanvas *c0 = new TCanvas("c0","c0",10,10,600,600);
  histEcalIso->Draw("colz");

  TCanvas *c1 = new TCanvas("c1","c1",100,10,600,600);
  TGraph *grEcal = runUCCOM(histEcalNum, histEcalDen);
  grEcal->Draw("ALP");

  TCanvas *c2 = new TCanvas("c2","c2",200,10,600,600);
  TGraph *grHcal = runUCCOM(histHcalNum, histHcalDen);
  grHcal->Draw("ALP");

  TCanvas *c3 = new TCanvas("c3","c3",300,10,600,600);
  TGraph *grTrk = runUCCOM(histTrkNum, histTrkDen);
  grTrk->Draw("ALP");

  TCanvas *c4 = new TCanvas("c4","c4",400,10,600,600);
  TGraph *grFull = runUCCOM(histFullNum, histFullDen);
  grFull->Draw("ALP");

}


int findOfflineEtaBin(float eta){

  int result = -1;
  float abseta = fabs(eta);
  for( int i=0; i<nEtaBins; i++){
    if( abseta >= etaBinLimits[i] && abseta < etaBinLimits[i+1] ){
      result = i;
      break;
    }
  }
  return result;
}

int findHltEtaBin(float eta){

  int result = -1;
  float abseta = fabs(eta);
  for( int i=0; i<nHltEtaBins; i++){
    if( abseta >= hltEtaBinLimits[i] && abseta < hltEtaBinLimits[i+1] ){
      result = i;
      break;
    }
  }
  return result;
}

TGraph *runUCCOM(TH1F *hnum, TH1F *hden){

  // Build the UCCOM curve
  TGraph *gr = new TGraph();
  // loop over offline cut values
  for(int ix = 1; ix <= nIsoBins; ix++){ // include overflows
    float offlineCut = hden->GetXaxis()->GetBinUpEdge(ix);
    float passOffline = hden->Integral(1, ix);
    float passOfflineAndNotHLT = hnum->Integral(1, ix);
    float fraction = -1;
    if( passOffline > 0)
      fraction = passOfflineAndNotHLT / passOffline;
    gr->SetPoint(gr->GetN(), offlineCut, fraction);
  }

  return gr;
}
