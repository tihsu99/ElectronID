#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

const TString fileNameSignal = "DYJetsToLL_may29_flat_ntuple_1M_barrel.root";
const TString fileNameBackground = "TTJets_may29_flat_ntuple_1M_barrel.root";
const TString treeName = "electronTree";

// Forward declarations
TTree *getTree(TString fname, TString tname);

// Main function
void testReweighting(){

  TTree *treeS = getTree(fileNameSignal, treeName);
  TTree *treeB = getTree(fileNameBackground, treeName);

  //
  // Test pt reweighting
  //

  TH1F *hptS = new TH1F("hptS","",180,20,180);
  TH1F *hptSweighted = new TH1F("hptSweighted","",180,20,180);
  TH1F *hptBweighted = new TH1F("hptBweighted","",180,20,180);
  
  TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);

  treeS->Draw("pt>>hptS","isTrueEle==1 && passConversionVeto==1 && abs(dz)<1");
  treeS->Draw("pt>>hptSweighted","kinWeight*(isTrueEle==1 && passConversionVeto==1 && abs(dz)<1)");
  treeB->Draw("pt>>hptBweighted",
	      "kinWeight*( (isTrueEle==0 || isTrueEle==3) && passConversionVeto==1 && abs(dz)<1)");

  hptBweighted->Scale(hptSweighted->GetSumOfWeights()/hptBweighted->GetSumOfWeights());

  hptS        ->SetLineColor(kRed);
  hptSweighted->SetLineColor(kBlue);
  hptS        ->SetLineWidth(2);
  hptSweighted->SetLineWidth(2);
  hptBweighted->SetMarkerStyle(20);
  hptBweighted->SetMarkerSize(1);
  hptBweighted->GetXaxis()->SetTitle("p_{T} [GeV]");

  hptBweighted->Draw("pe");
  hptS        ->Draw("hist,same");
  hptSweighted->Draw("hist,same");
  
  //
  // Test eta reweighting
  //

  TH1F *hetaS = new TH1F("hetaS","",50,-2.5,2.5);
  TH1F *hetaSweighted = new TH1F("hetaSweighted","",50,-2.5,2.5);
  TH1F *hetaBweighted = new TH1F("hetaBweighted","",50,-2.5,2.5);
  
  TCanvas *c2 = new TCanvas("c2","c2",500,10,600,600);

  treeS->Draw("etaSC>>hetaS","isTrueEle==1 && passConversionVeto==1 && abs(dz)<1");
  treeS->Draw("etaSC>>hetaSweighted",
	      "kinWeight*(isTrueEle==1 && passConversionVeto==1 && abs(dz)<1)");
  treeB->Draw("etaSC>>hetaBweighted",
	      "kinWeight*( (isTrueEle==0 || isTrueEle==3) && passConversionVeto==1 && abs(dz)<1)");

  hetaBweighted->Scale(hetaSweighted->GetSumOfWeights()/hetaBweighted->GetSumOfWeights());

  hetaS        ->SetLineColor(kRed);
  hetaSweighted->SetLineColor(kBlue);
  hetaS        ->SetLineWidth(2);
  hetaSweighted->SetLineWidth(2);
  hetaBweighted->SetMarkerStyle(20);
  hetaBweighted->SetMarkerSize(1);
  hetaBweighted->GetXaxis()->SetTitle("#eta");

  hetaBweighted->Draw("pe");
  hetaS        ->Draw("hist,same");
  hetaSweighted->Draw("hist,same");
  


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
