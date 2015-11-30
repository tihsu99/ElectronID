#include "TLatex.h"
#include "TTree.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TCanvas.h"
#include "TString.h"
#include "TH1F.h"
#include "TFile.h"

#include "OptimizationConstants.hh"
#include "VarCut.hh"

#include <iostream>
using std::cout;

// Draw barrel or endcap
const bool drawBarrel = true;
const int nROCs = 4;

const bool verbose = true;

// Files with signal and background trees (ideally the ntuples
// that were used for TMVA optimization
const TString fnameSignal = "/afs/cern.ch/user/r/rkamalie/workspace/public/DY_Run2Asympt25ns_miniAOD_sept18_flat_ntuple_withWeights.root";
const TString signalTreeName = "electronTree";
const TString fnameBackground = "/afs/cern.ch/user/r/rkamalie/workspace/public/TT_Run2Asympt25ns_miniAOD_sept21_flat_ntuple_withWeights.root";
const TString backgroundTreeName = "electronTree";

// Name TMVA output file that contains the pre-computed ROC, etc
const TString tmvaFileNameBarrel[nROCs] = {
  //only 1st file was optimized and will be used
  "trainingData/training_results_barrel_pass1_20151101_100000/TMVA_training_results_barrel_pass1_20151101_100000.root",
  "trainingData/training_results_barrel_pass1_20151101_100000/TMVA_training_results_barrel_pass1_20151101_100000.root",
  "trainingData/training_results_barrel_pass1_20151101_100000/TMVA_training_results_barrel_pass1_20151101_100000.root",
  "trainingData/training_results_barrel_pass1_20151101_100000/TMVA_training_results_barrel_pass1_20151101_100000.root"

};

const TString tmvaFileNameEndcap[nROCs] ={
   //no optimization was done for endcap, only for aestetic purpose
  "trainingData/training_results_barrel_pass1_20151101_100000/TMVA_training_results_barrel_pass1_20151101_100000.root",
  "trainingData/training_results_barrel_pass1_20151101_100000/TMVA_training_results_barrel_pass1_20151101_100000.root",
  "trainingData/training_results_barrel_pass1_20151101_100000/TMVA_training_results_barrel_pass1_20151101_100000.root",
  "trainingData/training_results_barrel_pass1_20151101_100000/TMVA_training_results_barrel_pass1_20151101_100000.root"

};

//
// Names of ROOT files that contain working points
//
const int nWorkingPointSets = 4;
// Set 1
const int markerColorSet1 = kRed;
const int markerStyleSet1 = 20;
const TString legendSet1 = "WP_Veto";
const int nWP = 1;
const TString cutFileNamesBarrelSet1[nWP] = { 
  "cut_repository/cuts_barrel_20151101_100000_WP_Veto.root"

};
const TString cutFileNamesEndcapSet1[nWP] = {
  "cut_repository/cuts_barrel_20151101_100000_WP_Veto.root"

};

// Set 2
const int markerColorSet2 = kOrange;
const int markerStyleSet2 = 20;
const TString legendSet2 = "WP_Loose";
const TString cutFileNamesBarrelSet2[nWP] = { 
  "cut_repository/cuts_barrel_20151101_100000_WP_Veto.root"

};
const TString cutFileNamesEndcapSet2[nWP] = {
  "cut_repository/cuts_barrel_20151101_100000_WP_Veto.root"
};

// Set 3
const int markerColorSet3 = kBlue;
const int markerStyleSet3 = 20;
const TString legendSet3 = "WP_Medium";
const TString cutFileNamesBarrelSet3[nWP] = { 
  "cut_repository/cuts_barrel_20151101_100000_WP_Veto.root"

};
const TString cutFileNamesEndcapSet3[nWP] = {
  "cut_repository/cuts_barrel_20151101_100000_WP_Veto.root"
};



// Set 4
const int markerColorSet4 = kGreen;
const int markerStyleSet4 = 20;
const TString legendSet4 = "WP_Tight";
const TString cutFileNamesBarrelSet4[nWP] = { 
  "cut_repository/cuts_barrel_20151101_100000_WP_Veto.root"
};
const TString cutFileNamesEndcapSet4[nWP] = {
  "cut_repository/cuts_barrel_20151101_100000_WP_Veto.root"
};

void bazinga (std::string message){
  if (verbose)
    std::cout<<message<<std::endl;
}

// Forward declarations
TTree *getTreeFromFile(TString fname, TString tname);
void overlayWorkingPoints(TCanvas *c1, 
			  TTree *signalTree, TTree *backgroundTree, 
			  const TString *cutFileNames,
			  int markerColor, int markerStyle, 
			  TLegend *leg, const TString legendText);
void   findEfficiencies(TTree *signalTree, TTree *backgroundTree,
			float &effSignal, float &effBackground, 
			VarCut *cutObject);
//
// Main function
//
void drawROCandWPv4(){
  
  bazinga("Step #0");
  TCanvas *c1 = new TCanvas("c1","",10,10,600,600);
  //  c1->cd();

  TString tmvaFileName1 = tmvaFileNameBarrel[0];
  TString tmvaFileName2 = tmvaFileNameBarrel[1];
  TString tmvaFileName3 = tmvaFileNameBarrel[2];
  TString tmvaFileName4 = tmvaFileNameBarrel[3];
  bazinga("Step #1");
  const TString *cutFileNamesSet1 = cutFileNamesBarrelSet1;
  const TString *cutFileNamesSet2 = cutFileNamesBarrelSet2;
  const TString *cutFileNamesSet3 = cutFileNamesBarrelSet3;
  const TString *cutFileNamesSet4 = cutFileNamesBarrelSet4;
  if( !drawBarrel ){
    tmvaFileName1 = tmvaFileNameEndcap[0];
    tmvaFileName2 = tmvaFileNameEndcap[1];
    tmvaFileName3 = tmvaFileNameEndcap[2];
    tmvaFileName4 = tmvaFileNameEndcap[3];

    cutFileNamesSet1 = cutFileNamesEndcapSet1;
    cutFileNamesSet2 = cutFileNamesEndcapSet2;
    cutFileNamesSet3 = cutFileNamesEndcapSet3;
    cutFileNamesSet4 = cutFileNamesEndcapSet4;
  }
  bazinga("Step #2");
  // Open the file with TMVA output
  TFile *tmvaFile1 = new TFile(tmvaFileName1);
  TFile *tmvaFile2 = new TFile(tmvaFileName2);
  TFile *tmvaFile3 = new TFile(tmvaFileName3);
  TFile *tmvaFile4 = new TFile(tmvaFileName4);
  if( !tmvaFile1 || !tmvaFile2 || !tmvaFile3 || !tmvaFile4)
    assert(0);
  bazinga("Step #3");
  // 
  //  Draw the ROC curve
  //
  TH1F *hROC1 = (TH1F*)tmvaFile1->Get("Method_Cuts/Cuts/MVA_Cuts_rejBvsS");
  TH1F *hROC2 = (TH1F*)tmvaFile2->Get("Method_Cuts/Cuts/MVA_Cuts_rejBvsS");
  TH1F *hROC3 = (TH1F*)tmvaFile3->Get("Method_Cuts/Cuts/MVA_Cuts_rejBvsS");
  TH1F *hROC4 = (TH1F*)tmvaFile4->Get("Method_Cuts/Cuts/MVA_Cuts_rejBvsS");
  if( !hROC1 || !hROC2 || !hROC3|| !hROC3 )
    assert(0);

  //c1->cd();
  //hROC1->Draw();
  bazinga("Step #4");
  
  TH1F * hROC = new TH1F("hROC", "", 100, 0, 1);
  hROC->SetDirectory(0);

  for(int iBin = 1; iBin <=100; ++iBin){

    if (drawBarrel){
    // barrrel
      if(1<=iBin && iBin<=75)
	hROC->SetBinContent(iBin, hROC3->GetBinContent(iBin)) ;
      else if(76 <= iBin && iBin<=89)
	hROC->SetBinContent(iBin, hROC2->GetBinContent(iBin)) ;
      else if(90 <=iBin && iBin<=100)
	hROC->SetBinContent(iBin, hROC1->GetBinContent(iBin)) ;
    }
    else {
      // for endcap
      if(1<=iBin && iBin<=70)
	hROC->SetBinContent(iBin, hROC3->GetBinContent(iBin)) ;
      else if(71 <= iBin && iBin<=89)
	hROC->SetBinContent(iBin, hROC2->GetBinContent(iBin)) ;
      else if(90 <=iBin && iBin<=100)
	hROC->SetBinContent(iBin, hROC1->GetBinContent(iBin)) ;
    }
    
  }
  
  bazinga("Step #5");
  
  //Set histogram attributes and draw the ROC curve
  hROC->SetStats(0);
  hROC->SetLineWidth(2);
  hROC->SetTitle("");
  hROC->GetXaxis()->SetTitle("signal efficiency");
  hROC->GetYaxis()->SetTitle("background rejection");
  hROC->GetYaxis()->SetTitleOffset(1.4);
  if( drawBarrel ){
    hROC->GetXaxis()->SetRangeUser(0.6, 1.0);
    hROC->GetYaxis()->SetRangeUser(0.951, 1.0);
  }else{
    hROC->GetXaxis()->SetRangeUser(0.6, 1.0);
    hROC->GetYaxis()->SetRangeUser(0.8, 1.0);
  }

  c1->cd();
  // hROC->Draw("L");
  hROC1->SetLineColor(kRed);
  hROC2->SetLineColor(kOrange);
  hROC3->SetLineColor(kBlue);
  hROC4->SetLineColor(kGreen);
  hROC1->GetYaxis()->SetRangeUser(0.9, 1.0);
  hROC1->SetStats(0);
  hROC1->Draw("L");
  //hROC1->Draw("Lsame");
  //hROC2->Draw("Lsame");
  //hROC3->Draw("Lsame");
  //hROC4->Draw("Lsame");
  bazinga("Step #6");
  TString commentText = "barrel electrons";
  if( !drawBarrel )
    commentText = "endcap electrons";
  TLatex *comment = new TLatex(0.2, 0.2, commentText);
  comment->SetNDC(kTRUE);
  comment->Draw();

  c1->Update();
  bazinga("Step #7");
  // 
  // Overlay the cuts
  //

  // First find the TestTree for measuring efficiency and rejection
  printf("\n Take true electrons from %s   tree %s\n", 
	 fnameSignal.Data(), signalTreeName.Data());
  TTree *signalTree = getTreeFromFile( fnameSignal, signalTreeName);
  // Input background tree  
  printf("\n Take background electrons from %s   tree %s\n", 
	 fnameBackground.Data(), backgroundTreeName.Data());
  TTree *backgroundTree = getTreeFromFile( fnameBackground, backgroundTreeName);
  bazinga("Step #8");
  // Next, draw all working point sets
  if( nWorkingPointSets==4 ){

    TLegend *leg = new TLegend(0.15, 0.45, 0.5, 0.7);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    overlayWorkingPoints(c1, signalTree, backgroundTree, cutFileNamesSet1, 
			 markerColorSet1, markerStyleSet1, leg, legendSet1);
    // bazinga("Step #9");
    // overlayWorkingPoints(c1, signalTree, backgroundTree, cutFileNamesSet2, 
    // 			 markerColorSet2, markerStyleSet2, leg, legendSet2);
    // bazinga("Step #10");
    // overlayWorkingPoints(c1, signalTree, backgroundTree, cutFileNamesSet3,
    // 			 markerColorSet3, markerStyleSet3, leg, legendSet3);
    // bazinga("Step #11");
    // overlayWorkingPoints(c1, signalTree, backgroundTree, cutFileNamesSet4,  			   
    // 			 markerColorSet4, markerStyleSet4, leg, legendSet4);
    // bazinga("Step #12");
    leg->Draw("same");
    
  }
  
  
    // Save the figure into a file
    TString outname = "figures/plot_ROCandWP_barrel.png";
    bazinga("Step #13");
    if( !drawBarrel )
      outname = "figures/plot_ROCandWP_endcap.png";
    c1->Print(outname);
    
    return;
};

// Get a given tree from a given file name.
TTree *getTreeFromFile(TString fname, TString tname){

  TFile *file = new TFile( fname );
  TTree *tree     = (TTree*) file->Get(tname);
  
  return tree;
}

// Draw on a given canvas the full set of working points
void   overlayWorkingPoints(TCanvas *c1, 
			    TTree *signalTree, TTree *backgroundTree, 
			    const TString *cutFileNames,
			    int markerColor, int markerStyle, 
			    TLegend *leg, const TString legendText){


  // Now loop over working points
  for(int iwp = 0; iwp<nWP; iwp++){
    
    // Load the working point from a ROOT file
    TFile *cutFile = new TFile(cutFileNames[iwp]);
    if( !cutFile )
      assert(0);
    VarCut *cutObject = (VarCut*)cutFile->Get("cuts");
    if( !cutObject )
      assert(0);
    
    // Compute the efficiencies
    float effSignal, effBackground;
    findEfficiencies(signalTree, backgroundTree, effSignal, effBackground,
		     cutObject);
    printf("Computed eff for cut from %s, effS= %.3f effB= %.3f\n",
	   cutFileNames[iwp].Data(), effSignal, effBackground);
    
    // Make a marker and draw it.
    TMarker *marker = new TMarker(effSignal, 1.0-effBackground, 20);
    marker->SetMarkerSize(2);
    marker->SetMarkerColor(markerColor);
    marker->SetMarkerStyle(markerStyle);
    marker->Draw("same");

    // Add marker to the legend only once. Do not draw the legend here,
    // it is drawn in the main function later
    if( iwp == 0 ){
      if( !leg )
	assert(0);
      leg->AddEntry(marker, legendText, "p");
    }

    c1->Update();
    
    cutFile->Close();
  }

  
}

// Compute signal and background efficiencies for given cuts
void findEfficiencies(TTree *signalTree, TTree *backgroundTree,
		      float &effSignal, float &effBackground, VarCut *cutObject){

  TCut etaCut = "";
  if( drawBarrel ){
    etaCut = Opt::etaCutBarrel;
  }else{
    etaCut = Opt::etaCutEndcap;
  }
  TCut kinematicCuts = Opt::ptCut && etaCut;

  TCut preselectionCuts = kinematicCuts && Opt::otherPreselectionCuts;
  
  TCut signalCuts = preselectionCuts && Opt::trueEleCut;
  TCut backgroundCuts = preselectionCuts && Opt::fakeEleCut;  
 
  TCut selectionCuts = *(cutObject->getCut());

  // printf("\nSelecton cuts:\n");
  // selectionCuts.Print();
  // printf("\nSignal cuts:\n");
  // signalCuts.Print();
  // printf("\nBackground cuts:\n");
  // backgroundCuts.Print();
  
  TH1F *hS_num = new TH1F("hS_num","",100,0.,10000.);
  TH1F *hS_den = new TH1F("hS_den","",100,0.,10000.);
  TH1F *hBG_num = new TH1F("hBG_num","",100,0.,10000.);
  TH1F *hBG_den = new TH1F("hBG_den","",100,0.,10000.);
 

  // draw the histogram
  signalTree->Draw("pt>>hS_num", "genWeight"*(selectionCuts && signalCuts), "goff", 100000);
  signalTree->Draw("pt>>hS_den", "genWeight"*(signalCuts), "goff", 100000);
  effSignal = hS_num->GetSumOfWeights()/ hS_den->GetSumOfWeights();
  
  // draw the histogram
  backgroundTree->Draw("pt>>hBG_num", "genWeight"*(selectionCuts && backgroundCuts), "goff", 100000);
  backgroundTree->Draw("pt>>hBG_den", "genWeight"*(backgroundCuts), "goff", 100000);
  effBackground = hBG_num->GetSumOfWeights()/ hBG_den->GetSumOfWeights();

  delete hS_num ; hS_num = nullptr;
  delete hS_den ; hS_den = nullptr;

  delete hBG_num ; hBG_num = nullptr;
  delete hBG_den ; hBG_den = nullptr;

 
  
  return;
}
