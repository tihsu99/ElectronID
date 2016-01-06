#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "stdio.h"  //printf
#include <iostream>  //cout
#include <math.h>
#include "TGraphErrors.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPad.h"
#include "VarCut.hh"
using namespace std;


const int smallCount = 100000;  // 1M is fine, if not - set to 10-100M
const TString variable_for_which_plot_eff = "etaSC";  // use etaSC or nPV, or pt

bool doBarrel = true;  // for etaSC keep it in the "barrel" mode
int nBins = 100; 


const TString etaRegion = (doBarrel)?"barrel":"endcap";
const TString fileOut = TString::Format("DY_Spring15_25ns_eff_vs_%s_%s_withWeight",variable_for_which_plot_eff.Data(),etaRegion.Data() );

enum class WPType { VETO=0, LOOSE, MEDIUM, TIGHT };
enum variableAndOptimizationRange  {ptOrPVBarrel=1 ,ptOrPVEndcap, etaFullRange};
TFile *DYfile= new TFile (fileOut + ".root","recreate");


//======================================================
const TString cutFileNameBarrel[] = {
  "root://cmsxrootd.fnal.gov///store/user/cmsdas/2016/SHORT_EXERCISES/ElectronsAndPhotons/cut_repository/cuts_barrel_20151101_100000_WP_Veto.root"
};

const TString cutFileNameEndcap[] = {
  // for aestetics putpose only, endcap has not been oprimized for the exercises 
  "root://cmsxrootd.fnal.gov///store/user/cmsdas/2016/SHORT_EXERCISES/ElectronsAndPhotons/cut_repository/cuts_barrel_20151101_100000_WP_Veto.root"
};

// for eff vs nPV 
const float binLimits[] = { 0,  7, 8, 9, 10, 11,
			    12, 13, 14, 15, 16,
			    17, 18, 19, 20, 21,
			    22, 23, 24, 25, 26,
			    27, 28, 30, 40 };




void calculateEffAndErrors(TH1F *phist, TH1F *thist, TH1F *&ehist){

  ehist->Divide(thist);
  double effSignal = 0.0;
  double effSignalErr2 = 0.0;
  double effSignalErr = 0.0;
  
  for(int i=1;i<=nBins;++i){
    //
    effSignal = phist->GetBinContent(i) / thist->GetBinContent(i);
    if (!TMath::Finite(ehist->GetBinContent(i) ) ) ehist->SetBinContent(i, 0.0);
    effSignalErr2 =  thist->GetBinError(i)*thist->GetBinError(i) * effSignal*(1-effSignal) /(thist->GetBinContent(i)*thist->GetBinContent(i) );
    if(effSignalErr2 <0.0 ) continue;
    effSignalErr = sqrt(effSignalErr2);
    if(!TMath::Finite(effSignalErr) ) continue;
    ehist->SetBinError(i, effSignalErr);    
  }
  
}

void setHistCosmetics(TH1F *&effh, int color, TString variable_for_which_plot_eff ){
 
  effh->SetTitle("");
  effh->GetYaxis()->SetTitle("Efficiency");
  //   effh->GetYaxis()->SetNdivisions(505);
  effh->GetYaxis()->SetTitleSize(20);
  effh->GetYaxis()->SetTitleFont(43);
  effh->GetYaxis()->SetTitleOffset(1.55);
  effh->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  effh->GetYaxis()->SetLabelSize(15);
  
  if (variable_for_which_plot_eff == "pt")
    effh->GetXaxis()->SetTitle("Pt,GeV");
  else if (variable_for_which_plot_eff == "etaSC")
    effh->GetXaxis()->SetTitle("#eta");
  else if(variable_for_which_plot_eff == "nPV")
    effh->GetXaxis()->SetTitle("nPV");
  effh->GetXaxis()->SetTitleSize(24);
  effh->GetXaxis()->SetTitleFont(43);
  effh->GetXaxis()->SetTitleOffset(1.2);
  effh->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  effh->GetXaxis()->SetLabelSize(15);
  effh->SetLineColor(color);
  effh->SetStats(0);
}


void applyCutsAndFillHist(const TString cutFileName, TH1F * histo, TTree * tree, bool denom, int passNumber, bool isBG, const TString cutFileName_EtaEndcap){

  TFile *cutFile = new TFile(cutFileName);
  if( !cutFile )
    assert(0);
  VarCut *cutObject = (VarCut*)cutFile->Get("cuts");
  if( !cutObject )
    assert(0);
  TCut selectionCuts = *(cutObject->getCut());
  TString drawLine = TString::Format("%s>>%s",variable_for_which_plot_eff.Data() , histo->GetName() );
  std::cout<<"\ndrawLine is as follows: "<<drawLine<<std::endl;
  std::cout<<"\nDenominator? "<<denom<<", passNumber is "<<passNumber<<std::endl;

  TFile *cutFile_EtaEndcap = new TFile(cutFileName_EtaEndcap);
  if( !cutFile_EtaEndcap )
    assert(0);
  VarCut *cutObject_EtaEndcap = (VarCut*)cutFile_EtaEndcap->Get("cuts");
  if( !cutObject_EtaEndcap )
    assert(0);
  TCut selectionCuts_EtaEndcap = *(cutObject_EtaEndcap->getCut());
  TString drawLine_EtaEndcap = TString::Format("%s>>+%s",variable_for_which_plot_eff.Data() , histo->GetName() );
  std::cout<<"\ndrawLine_EtaEndcap is as follows: "<<drawLine_EtaEndcap<<std::endl;

  
  TCut extra = "";
  TCut extra_EtaEndcap = "";

  // process signal first
  if (!isBG){
    if (passNumber == ptOrPVBarrel)
      extra = "( (pt >20) && ( abs(etaSC) < 1.4442  ) && passConversionVeto &&  (isTrueEle == 1) && (abs(dz) < 1) )" ;
    else if (passNumber == ptOrPVEndcap)
      extra = "( (pt >20)  && ( abs(etaSC) > 1.566  ) && ( abs(etaSC) < 2.5  ) && passConversionVeto &&  (isTrueEle == 1) && (abs(dz) < 1) )" ;
    else if (passNumber == etaFullRange){
      extra = "( (pt >20 ) && ( abs(etaSC) < 1.4442  ) &&passConversionVeto &&  (isTrueEle == 1) && (abs(dz) < 1) )" ;
      extra_EtaEndcap = "( (pt >20 ) && (abs(etaSC) > 1.566 && abs(etaSC) < 2.5  ) && passConversionVeto &&  (isTrueEle == 1) && (abs(dz) < 1) )" ;
    }
  } 
  else {
    if (passNumber == ptOrPVBarrel )
      extra = "( (pt >20 ) && ( abs(etaSC) < 1.4442  ) && passConversionVeto &&  (isTrueEle == 0 || isTrueEle == 3) && (abs(dz) < 1) )" ;
    else if (passNumber == ptOrPVEndcap)
      extra = "( (pt >20 ) && ( abs(etaSC) > 1.566  ) && ( abs(etaSC) < 2.5  ) && passConversionVeto &&  (isTrueEle == 0 || isTrueEle == 3) && (abs(dz) < 1) )" ;
    else if (passNumber == etaFullRange){
      extra = "( (pt >20 ) && (abs(etaSC) <  1.4442 ) && passConversionVeto &&  (isTrueEle == 0 || isTrueEle == 3) && (abs(dz) < 1) )" ;
      extra_EtaEndcap = "( (pt >20 ) && (abs(etaSC) > 1.566 && abs(etaSC) < 2.5  ) && passConversionVeto &&  (isTrueEle == 0 || isTrueEle == 3) && (abs(dz) < 1) )" ;

    }
  }
  
  Long64_t num = 0; 
  Long64_t den = 0;
  Double_t integral_num   = 0;
  Double_t effEntries_num = 0;
  Double_t integral_den   = 0;
  Double_t effEntries_den = 0;

  histo->GetDirectory()->cd();
  histo->Sumw2();

  // process numerator first
  if (!denom){ 

    num =  tree->Draw(drawLine,"genWeight"*(selectionCuts && extra), "goff", smallCount);
    if (variable_for_which_plot_eff == "etaSC") tree->Draw(drawLine_EtaEndcap,"genWeight"*(selectionCuts_EtaEndcap && extra_EtaEndcap), "goff", smallCount);
    
  }
  else if (denom){
    std::cout<<"\ncut is as follows: "<<extra<< " and "<<extra<<std::endl;
    den =tree->Draw(drawLine,"genWeight"*(extra), "goff", smallCount);
    if (variable_for_which_plot_eff == "etaSC") tree->Draw(drawLine_EtaEndcap,"genWeight"*(extra_EtaEndcap), "goff", smallCount);
  }
  else std::cout<<"\ncan't really happen!";

  std::cout<<"\nEntries is numerator: "<<num<<" and in numerator "<<den<<endl;

}


void useHistograms( TH1F *  hV,  const TString * pToCutsFile, TTree * tree, bool denom, int passNumber, bool isBG, const TString * pToCutsFile_EtaEndcap ){

  applyCutsAndFillHist(pToCutsFile[static_cast<int>(WPType::VETO)],   hV, tree, denom, passNumber, isBG, pToCutsFile_EtaEndcap[static_cast<int>(WPType::VETO)]);
 
}

void calculateEfficiencyFromNTUPLE_withGenWeights_v4(){

  if (variable_for_which_plot_eff == "etaSC") 
    doBarrel = true; // correct etaMode for eta variable
  
  TFile fi("root://cmsxrootd.fnal.gov///store/user/cmsdas/2016/SHORT_EXERCISES/ElectronsAndPhotons/samples/DY_Run2Asympt25ns_miniAOD_sept18_flat_ntuple_withWeights.root");   
  gDirectory->ls();
  TTree *tr          = (TTree*)gDirectory->Get("electronTree");

  TFile fiBG("root://cmsxrootd.fnal.gov///store/user/cmsdas/2016/SHORT_EXERCISES/ElectronsAndPhotons/samples/TT_Run2Asympt25ns_miniAOD_sept21_flat_ntuple_withWeights.root");   
  gDirectory->ls();
  TTree *trBG          = (TTree*)fiBG.Get("electronTree");
 

  float binLimitLow = 20;
  float binLimitUpper = 170;

  if(variable_for_which_plot_eff == "pt") {
    std::cout<<"\nVarible is  Pt, all is fine!" << std::endl;
  }
  else if(variable_for_which_plot_eff == "etaSC") {
    binLimitLow = -2.5;
    binLimitUpper = 2.5;
  }
  else if (variable_for_which_plot_eff == "nPV") {
    binLimitLow = 0;
    binLimitUpper = 50;
  } 
  else {
    std::cout<<"\n No no, can't happen!"<<"\nVarible is  "<<variable_for_which_plot_eff.Data() << std::endl;
  }
  

  TH1F *histVeto;
  TH1F *histVetoTot;
 
  TH1F *bg_histVeto;
  TH1F *bg_histVetoTot;


  if (variable_for_which_plot_eff != "nPV") {

    histVeto=new TH1F("histVeto","",nBins,binLimitLow, binLimitUpper);
    bg_histVeto=new TH1F("bg_histVeto","",nBins,binLimitLow, binLimitUpper);

    histVetoTot=new TH1F("histVetoTot","",nBins,binLimitLow, binLimitUpper);
    bg_histVetoTot=new TH1F("bg_histVetoTot","",nBins,binLimitLow, binLimitUpper);

    
  }
  else {
    nBins =24;
    histVeto=new TH1F("histVeto","",nBins,binLimits);
    bg_histVeto=new TH1F("bg_histVeto","",nBins,binLimits);
     
    histVetoTot=new TH1F("histVetoTot","",nBins,binLimits);
    bg_histVetoTot=new TH1F("bg_histVetoTot","",nBins,binLimits);
    
  }
  
  TCanvas *c = new TCanvas ("c", "c", 600, 600);
  c->SetFillColor(kWhite);
  gPad->SetLeftMargin(0.12);
  gPad->SetBottomMargin(0.11);

  TLatex Tl05;
  Tl05.SetTextAlign(23);
  Tl05.SetTextSize(0.04);
 

  const TString *pToCutsFile = nullptr;
  const TString *pToCutsFile_EtaEndcap = nullptr;
  if (doBarrel ) {
    pToCutsFile = cutFileNameBarrel;
    pToCutsFile_EtaEndcap = cutFileNameEndcap;
  }
  else {
    pToCutsFile = cutFileNameEndcap;
    pToCutsFile_EtaEndcap = cutFileNameEndcap;  //is not used here, just for aestetics purpose
  }


  if ( (variable_for_which_plot_eff == "pt" || variable_for_which_plot_eff == "nPV" ) && (doBarrel ==true) ){
    //      if (Pt <20 || Pt >200 || abs(EtaSC) > 1.566 ) continue; 
    std::cout<<"\nI'm here 3 "<<std::endl;
    // signal
    // pt/nPV barrel
    useHistograms( histVeto,    pToCutsFile, tr, false, ptOrPVBarrel, false,  pToCutsFile_EtaEndcap);  
    useHistograms( histVetoTot, pToCutsFile, tr, true, ptOrPVBarrel, false,  pToCutsFile_EtaEndcap);

    std::cout<<"\nI'm here 4 "<<std::endl;
    // BG
    // pt/nPV barrel
    useHistograms( bg_histVeto,    pToCutsFile, trBG, false, ptOrPVBarrel, true,  pToCutsFile_EtaEndcap); 
    useHistograms( bg_histVetoTot, pToCutsFile, trBG, true, ptOrPVBarrel, true,  pToCutsFile_EtaEndcap);  


  }
  else if ( (variable_for_which_plot_eff == "pt" || variable_for_which_plot_eff == "nPV" ) && doBarrel ==false){

    // signal
    // pt/nPV endcap
    useHistograms( histVeto,     pToCutsFile, tr, false, ptOrPVEndcap, false,  pToCutsFile_EtaEndcap); 
    useHistograms( histVetoTot,  pToCutsFile, tr, true, ptOrPVEndcap, false,  pToCutsFile_EtaEndcap); 


    // BG
    // pt/nPV endcap
    useHistograms( bg_histVeto,    pToCutsFile, trBG, false, ptOrPVEndcap, true,  pToCutsFile_EtaEndcap);  
    useHistograms( bg_histVetoTot, pToCutsFile, trBG, true, ptOrPVEndcap, true,  pToCutsFile_EtaEndcap);

    
  }
  else if (variable_for_which_plot_eff == "etaSC"){

    // signal
    // etaFullRange
    useHistograms( histVeto,     pToCutsFile, tr, false, etaFullRange, false, pToCutsFile_EtaEndcap); 
    useHistograms( histVetoTot,  pToCutsFile, tr, true, etaFullRange, false, pToCutsFile_EtaEndcap);  

    // BG
    // etaFullRange
    useHistograms( bg_histVeto,    pToCutsFile, trBG, false, etaFullRange, true, pToCutsFile_EtaEndcap);
    useHistograms( bg_histVetoTot, pToCutsFile, trBG, true, etaFullRange, true, pToCutsFile_EtaEndcap); 
    
  }
  else {
    std::cout<<"\nCan't happen! \nVariable is  "<<variable_for_which_plot_eff<<"  and eta is "<< (doBarrel?"barrel":"endcap")<<std::endl;
  }

  std::cout<<"\nI'm here 6 "<<std::endl;
 
  
  TH1F *effV = (TH1F*)histVeto->Clone("effV");
  TH1F *effV_bg = (TH1F*)bg_histVeto->Clone("effV_bg");

  std::cout<<"\nI'm here 7 "<<std::endl;
  
  calculateEffAndErrors(histVeto, histVetoTot, effV);
  calculateEffAndErrors(bg_histVeto, bg_histVetoTot, effV_bg);


  std::cout<<"\nI'm here 8 "<<std::endl;
  //  Usage :
  //  histogram, color_number, variable
  //
  setHistCosmetics(effV,4, variable_for_which_plot_eff);
  setHistCosmetics(effV_bg,38, variable_for_which_plot_eff);
  

  //------------------------------------------------------------------

  //if (variable_for_which_plot_eff == "etaSC" || variable_for_which_plot_eff == "nPV")   

  effV->GetYaxis()->SetRangeUser(0.,1.);
  effV->Draw("L");
  effV_bg->Draw("Lsame");
  effV->SetLineWidth(3);
  effV_bg->SetLineWidth(3);

  TLegend *leg = new TLegend(0.35, 0.52, 0.75, 0.65);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(effV, "WP Veto", "lf"); 
  leg->AddEntry((TObject*)0, "", "");
  leg->AddEntry(effV_bg, "BG WP Veto", "lf"); 
  leg->SetTextSize(0.03);
  leg->Draw("same");


  if (variable_for_which_plot_eff == "pt" ){
    if (doBarrel)
      Tl05.DrawLatex(40,1.1,"barrel electrons");
    else
      Tl05.DrawLatex(40,1.1,"endcap electrons");     
  }
  else if ( variable_for_which_plot_eff == "nPV" ){
    if (doBarrel)
      Tl05.DrawLatex(30,1.05,"barrel electrons");
    else
      Tl05.DrawLatex(30,1.05,"endcap electrons");     
  }


  c->Modified();
  c->Update();
  c->SaveAs( fileOut + ".png");
  
  effV->SetDirectory(0);
  effV_bg->SetDirectory(0);
  
  DYfile->cd();
  effV->Write("WP_Veto");
  effV_bg->Write("BG_WP_Veto");

  c->Write();
  DYfile->Close();
}













