#include "TMath.h"
#include "TStyle.h"
#include "TH2D.h"
#include "TArrow.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TColor.h"
#include "TGaxis.h"
#include "OptimizationConstants.hh"
#include <cassert>
#include <TROOT.h>

// For debug purposes, set the flag below to true, for regular
// computation set it to false
const bool useSmallEventCount = false;

// Note: the DY ntuple has to be prepared for true electrons only to
// have meaningful kinematic weights.
TString tagDir = "2017-11-16";
TString signalFileName     = tagDir + "/DYJetsToLL_flat_ntuple_true_alleta_full.root";
TString backgroundFileName = tagDir + "/TTJets_flat_ntuple_trueAndFake_alleta_full.root";

// Forward declarations
TTree *getTreeFromFile(TString fname, TString tname);

// Main method
void drawKinematics(){

  TTree *signalTree     = getTreeFromFile( signalFileName, Opt::signalTreeName);
  TTree *backgroundTree = getTreeFromFile( backgroundFileName, Opt::backgroundTreeName);

  int maxEventsS = signalTree->GetEntries();
  int maxEventsB = signalTree->GetEntries();
  if( useSmallEventCount ){
    maxEventsS = 100000;
    maxEventsB = 100000;
  }

  //
  // Pt plots
  //

  TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
  gStyle->SetOptStat(0);
  c1->cd();

  TH1F *hptSig  = new TH1F("hptSig","",100, 0, 200);
  TH1F *hptBg   = new TH1F("hptBg","",100, 0, 200);
  TH1F *hptSigW = new TH1F("hptSigW","",100, 0, 200);
  
  signalTree->Draw("pt>>hptSig",
		   "pt>20 && isTrueEle==1 && abs(dz)<1 && passConversionVeto",
		   "goff", maxEventsS);
  backgroundTree->Draw("pt>>hptBg",
		       "pt>20 && (isTrueEle==0 || isTrueEle==3) && abs(dz)<1 && passConversionVeto", 
		       "goff", maxEventsB);
  signalTree->Draw("pt>>hptSigW","kinWeight*(pt>20 && isTrueEle==1 && abs(dz)<1 && passConversionVeto)", 
		   "goff", maxEventsS);

  hptSig->Scale(hptBg->GetSumOfWeights()/hptSig->GetSumOfWeights());
  hptSigW->Scale(hptBg->GetSumOfWeights()/hptSigW->GetSumOfWeights());

  hptBg->SetLineWidth(2);
  hptBg->GetXaxis()->SetTitle("p_{T} [GeV]");
  hptBg->Draw("hist");

  hptSig->SetLineColor(kRed);
  hptSig->SetLineWidth(2);
  hptSig->Draw("same,hist");

  hptSigW->SetLineColor(kBlue);
  hptSigW->SetLineWidth(2);
  hptSigW->SetMarkerStyle(20);
  hptSigW->SetMarkerSize(1);
  hptSigW->SetMarkerColor(kBlue);
  hptSigW->Draw("same,pe");

  TLegend *leg = new TLegend(0.5, 0.5, 0.8, 0.8);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hptSig, "signal", "l");
  leg->AddEntry(hptBg, "background", "l");
  leg->AddEntry(hptSigW, "signal reweighted", "pl");
  leg->Draw("same");

  c1->Print("figures/plot_kinematics_pt.png");

  //
  // Eta plots
  //

  TCanvas *c2 = new TCanvas("c2","c2",10,10,600,600);
  gStyle->SetOptStat(0);
  c2->cd();

  TH1F *hetaSig  = new TH1F("hetaSig","",50, -2.5, 2.5);
  TH1F *hetaBg   = new TH1F("hetaBg","", 50, -2.5, 2.5);
  TH1F *hetaSigW = new TH1F("hetaSigW","",50, -2.5, 2.5);
  
  signalTree->Draw("etaSC>>hetaSig",
		   "pt>20 && isTrueEle==1 && abs(dz)<1 && passConversionVeto",
		   "goff", maxEventsS);
  backgroundTree->Draw("etaSC>>hetaBg",
		       "pt>20 && (isTrueEle==0 || isTrueEle==3) && abs(dz)<1 && passConversionVeto", 
		       "goff", maxEventsB);
  signalTree->Draw("etaSC>>hetaSigW","kinWeight*(pt>20 && isTrueEle==1 && abs(dz)<1 && passConversionVeto)", 
		   "goff", maxEventsS);

  hetaSig->Scale(hetaBg->GetSumOfWeights()/hetaSig->GetSumOfWeights());
  hetaSigW->Scale(hetaBg->GetSumOfWeights()/hetaSigW->GetSumOfWeights());

  hetaBg->SetLineWidth(2);
  hetaBg->GetXaxis()->SetTitle("#eta_{SC}");
  hetaBg->Draw("hist");

  hetaSig->SetLineColor(kRed);
  hetaSig->SetLineWidth(2);
  hetaSig->Draw("same,hist");

  hetaSigW->SetLineColor(kBlue);
  hetaSigW->SetLineWidth(2);
  hetaSigW->SetMarkerStyle(20);
  hetaSigW->SetMarkerSize(1);
  hetaSigW->SetMarkerColor(kBlue);
  hetaSigW->Draw("same,pe");

  TLegend *leg1 = new TLegend(0.3, 0.15, 0.7, 0.5);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->AddEntry(hetaSig, "signal", "l");
  leg1->AddEntry(hetaBg, "background", "l");
  leg1->AddEntry(hetaSigW, "signal reweighted", "pl");
  leg1->Draw("same");

  c2->Print("figures/plot_kinematics_eta.png");
}

// Get a given tree from a given file name.
TTree *getTreeFromFile(TString fname, TString tname){

  TFile *file = new TFile( fname );
  TTree *tree     = (TTree*) file->Get(tname);
  
  return tree;
}

// Compiled
int main(int argc, char *argv[]){
  gROOT->SetBatch();
  drawKinematics();
}
