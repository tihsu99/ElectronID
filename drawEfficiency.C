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
#include "VarCut.hh"

// For debug purposes, set the flag below to true, for regular
// computation set it to false
const bool useSmallEventCount = false;

// Draw barrel or endcap
const bool drawBarrel = true;
const bool longPtRange = true;
const bool scaleBackground = false;

enum Mode {EFF_PT=0, EFF_PT_2TEV, EFF_ETA, EFF_NVTX};
TString varName[4] = {"pt", "pt", "etaSC", "nPV"};

Mode mode = EFF_PT_2TEV; // If mode is EFF_ETA, set drawBarrel to true!

TString dateTag = "2017-11-07";

// File name with working point cuts
const TString cutFileNamesBarrel[4] = { 
  "./cut_repository/cuts_barrel_" + dateTag + "_WP_Veto.root",
  "./cut_repository/cuts_barrel_" + dateTag + "_WP_Loose.root",
  "./cut_repository/cuts_barrel_" + dateTag + "_WP_Medium.root",
  "./cut_repository/cuts_barrel_" + dateTag + "_WP_Tight.root"
};
const TString cutFileNamesEndcap[4] = {
  "./cut_repository/cuts_endcap_" + dateTag + "_WP_Veto.root",
  "./cut_repository/cuts_endcap_" + dateTag + "_WP_Loose.root",
  "./cut_repository/cuts_endcap_" + dateTag + "_WP_Medium.root",
  "./cut_repository/cuts_endcap_" + dateTag + "_WP_Tight.root"
};


// Cuts on expected missing hits are separate from VarCut cuts, tuned by hand.
// These are summer 2016 values:
const int missingHitsCutBarrel[Opt::nWP] = { 2, 1, 1, 1};
const int missingHitsCutEndcap[Opt::nWP] = { 3, 1, 1, 1};

// Signal and background files
// TString signalFileName = drawBarrel ? Opt::fnameSignalBarrel : Opt::fnameSignalEndcap;
// TString backgroundFileName = drawBarrel ? Opt::fnameBackgroundBarrel : Opt::fnameBackgroundEndcap;
TString signalFileName     = "DYJetsToLL_flat_ntuple_true_alleta_full.root";
TString backgroundFileName = "TTJets_flat_ntuple_trueAndFake_alleta_full.root";


// Histogram settings
const int nPtBinsMax = 1000; // a large number
double ptBinLimits[nPtBinsMax+1];

const int nEtaBinsMax = 1000; // a large number
double etaBinLimits[nEtaBinsMax+1];

const int nNvtxBins = 28;
double nvtxBinLimits[nNvtxBins+1] = 
  {0, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
   21, 22, 23, 24, 25, 27, 30, 35, 40, 45, 50};

int sigColors[Opt::nWP] = {2, 4, 6, 8};
int  bgColors[Opt::nWP] = {9, 46, 30, 42};

const TString sigLegString[Opt::nWP] = {"WP Veto", "WP Loose", "WP Medium", "WP Tight"};
const TString  bgLegStringScaled[Opt::nWP] = {"WP Veto #epsilon x5", 
					      "WP Loose #epsilon x5", 
					      "WP Medium #epsilon x5", 
					      "WP Tight #epsilon x5"};
const TString  bgLegStringUnscaled[Opt::nWP] = {"WP Veto", 
						"WP Loose", 
						"WP Medium", 
						"WP Tight"};
const TString *bgLegString = scaleBackground ? bgLegStringScaled : bgLegStringUnscaled;


// Forward declarations
TTree * getTreeFromFile(TString fname, TString tname);
TH1F*   calculateEffAndErrors(TH1F *numHist, TH1F *denHist, TString eName, double* binLimits);
double *getPtBinLimits(int &nBins, bool longPtRange);
double *getEtaBinLimits(int &nBins);
double *getNvtxBinLimits(int &nBins);
void    setHistogram(TH1F *hist, int wp, bool isSignal);

// Get histogram for selectionsCuts and varName into histName
TH1F* drawFromTree(TTree* tree, TCut selectionCuts, TString varName, TString histName, int nBins, double* binLimits, int maxEvents){
  TH1F* hist = new TH1F(histName , histName, nBins, binLimits);
  hist->Sumw2();
  TString cutString = TString::Format("genWeight*kinWeight*(%s)", selectionCuts.GetTitle());
  TString command    = TString::Format("%s>>%s", varName.Data(), histName.Data());
  tree->Draw(command, cutString, "goff", maxEvents);
  return hist;
}


// Main function
void drawEfficiency(bool drawBarrel){

  const TString *cutFileNames = drawBarrel ? cutFileNamesBarrel : cutFileNamesEndcap;
  const int *missingHitsCut = drawBarrel ? missingHitsCutBarrel : missingHitsCutEndcap;

  if(mode == EFF_PT_2TEV) signalFileName = "DoubleEleFlat_flat_ntuple_trueAndFake_alleta_full.root";
  TTree *signalTree     = getTreeFromFile( dateTag + "/" + signalFileName, Opt::signalTreeName);
  TTree *backgroundTree = getTreeFromFile( dateTag + "/" + backgroundFileName, Opt::backgroundTreeName);

  int maxEventsS = signalTree->GetEntries();
  int maxEventsB = signalTree->GetEntries();
  if( useSmallEventCount ){
    maxEventsS = 100000;
    maxEventsB = 100000;
  }

  //
  // Draw efficiency
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
  gStyle->SetOptStat(0);
  c1->cd();

  TH1F *sigNum[Opt::nWP];
  TH1F *sigDen[Opt::nWP];
  TH1F *sigEff[Opt::nWP];

  TH1F *bgNum[Opt::nWP];
  TH1F *bgDen[Opt::nWP];
  TH1F *bgEff[Opt::nWP];

  int nBins;
  double *binLimits = nullptr;
  TString axisLabel = "";
  if( mode == EFF_PT or mode == EFF_PT_2TEV){
    binLimits = getPtBinLimits(nBins, mode==EFF_PT_2TEV);
    axisLabel = "p_{T} [GeV]";
  }else if( mode == EFF_ETA ){
    binLimits = getEtaBinLimits(nBins);
    axisLabel = "#eta_{SC}";
  }else if(mode == EFF_NVTX){
    binLimits = getNvtxBinLimits(nBins);
    axisLabel = "Nvtx";
  }else{
    printf("ERROR: unknown mode\n");
    assert(0);
  }

  // for(int i=0; i<=nBins; i++)
  //   printf("%2d    %5.0f\n", i, binLimits[i]);
  TH2D *dummy = new TH2D("dummy","", 100, binLimits[0], binLimits[nBins],
			 100, 0, 1);


  dummy->GetXaxis()->SetTitle(axisLabel);
  dummy->GetXaxis()->SetTitleOffset(1.2);
  dummy->GetYaxis()->SetTitle("efficiency");
  dummy->GetYaxis()->SetTitleOffset(1.2);
  dummy->Draw();

  // Define denominator cuts
  TString comment = "barrel electrons";
  TCut preselectionCuts = Opt::ptCut;
  if( mode == EFF_ETA ){
    comment = "";
    TCut anyEtaCut = "abs(etaSC)<2.5";
    preselectionCuts += anyEtaCut;
  }else{
    if( drawBarrel ){
      comment = "barrel electrons";
      preselectionCuts += Opt::etaCutBarrel;
    }else{
      comment = "endcap electrons";
      preselectionCuts += Opt::etaCutEndcap;
    }
  }
  preselectionCuts += Opt::otherPreselectionCuts;

  TCut signalCuts     = (mode == EFF_PT_2TEV ? preselectionCuts : preselectionCuts && Opt::trueEleCut);
  TCut backgroundCuts = preselectionCuts && Opt::fakeEleCut;

  // Load cut files
  TFile *fileCut[Opt::nWP];
  TFile *fileCutBarrel[Opt::nWP];
  TFile *fileCutEndcap[Opt::nWP];
  VarCut *wpCutsObject[Opt::nWP];
  TCut wpCuts[Opt::nWP];
  TCut wpCutsBarrel[Opt::nWP];
  TCut wpCutsEndcap[Opt::nWP];
  for(int i=0; i<Opt::nWP; i++){
    fileCut[i] = new TFile( cutFileNames[i] );
    wpCutsObject[i] = (VarCut*)fileCut[i]->Get("cuts");
    wpCuts[i] = *(wpCutsObject[i]->getCut());
    fileCut[i]->Close();
    // 
    fileCutBarrel[i] = new TFile( cutFileNamesBarrel[i] );
    wpCutsObject[i] = (VarCut*)fileCutBarrel[i]->Get("cuts");
    wpCutsBarrel[i] = *(wpCutsObject[i]->getCut());
    fileCutBarrel[i]->Close();
    //
    fileCutEndcap[i] = new TFile( cutFileNamesEndcap[i] );
    wpCutsObject[i] = (VarCut*)fileCutEndcap[i]->Get("cuts");
    wpCutsEndcap[i] = *(wpCutsObject[i]->getCut());
    fileCutEndcap[i]->Close();
  }

  for(int i=0; i<Opt::nWP; i++){
  //for(int i=0; i<2; i++){
    printf("Process working point %d\n", i);

    // Book all histograms
    TString sigNumHistName = TString::Format("sigNum_wp%d", i);
    TString sigDenHistName = TString::Format("sigDen_wp%d", i);
    TString sigEffHistName = TString::Format("sigEff_wp%d", i);
    /*
    sigNum[i] = new TH1F(sigNumHistName , sigNumHistName, nBins, binLimits);
    sigDen[i] = new TH1F(sigDenHistName , sigDenHistName, nBins, binLimits);
    sigEff[i] = new TH1F(sigEffHistName , sigEffHistName, nBins, binLimits);
    sigNum[i]->Sumw2();
    sigDen[i]->Sumw2();
    sigEff[i]->Sumw2();
*/
    TString bgNumHistName = TString::Format("bgNum_wp%d", i);
    TString bgDenHistName = TString::Format("bgDen_wp%d", i);
    TString bgEffHistName = TString::Format("bgEff_wp%d", i);
    /*
    bgNum[i] = new TH1F(bgNumHistName, bgNumHistName, nBins, binLimits);
    bgDen[i] = new TH1F(bgDenHistName, bgDenHistName, nBins, binLimits);
    bgEff[i] = new TH1F(bgEffHistName, bgEffHistName, nBins, binLimits);
    bgNum[i]->Sumw2();
    bgDen[i]->Sumw2();
    bgEff[i]->Sumw2();
*/
    // Set up cuts
    TCut selectionCuts;
    TCut missingHits;
    printf("\nNOTE: the missing hits cuts are not taken from optimization, but are added by hand!\n\n");
    if( mode != EFF_ETA ){
      selectionCuts += wpCuts[i];
      missingHits = TString::Format("expectedMissingInnerHits<=%d", missingHitsCut[i]).Data();
      selectionCuts = selectionCuts && missingHits;
    } else {
      TCut selectionCutsBarrel;
      selectionCutsBarrel += wpCutsBarrel[i];
      TCut missingHitsBarrel = TString::Format("expectedMissingInnerHits<=%d", missingHitsCutBarrel[i]).Data();
      selectionCutsBarrel = selectionCutsBarrel && missingHitsBarrel && Opt::etaCutBarrel;

      TCut selectionCutsEndcap = wpCutsEndcap[i];
      TCut missingHitsEndcap = TString::Format("expectedMissingInnerHits<=%d", missingHitsCutEndcap[i]).Data();;
      selectionCutsEndcap = selectionCutsEndcap && missingHitsEndcap && Opt::etaCutEndcap;
      //
      selectionCuts = (selectionCutsBarrel) || (selectionCutsEndcap);
    }    

    TString command;
    // Signal numerator and denominator
    sigNum[i] = drawFromTree(signalTree,     (selectionCuts && signalCuts),     varName[mode], TString::Format("sigNum_wp%d", i), nBins, binLimits, maxEventsS);
    sigDen[i] = drawFromTree(signalTree,     (signalCuts),                      varName[mode], TString::Format("sigDen_wp%d", i), nBins, binLimits, maxEventsS);
    bgNum[i]  = drawFromTree(backgroundTree, (selectionCuts && backgroundCuts), varName[mode], TString::Format("bgNum_wp%d", i),  nBins, binLimits, maxEventsB);
    bgDen[i]  = drawFromTree(backgroundTree, (backgroundCuts),                  varName[mode], TString::Format("bgDen_wp%d", i),  nBins, binLimits, maxEventsB);

    // Compute efficiencies
    if( scaleBackground ){
      printf("\n\nSCALE BACKGROUND EFF BY x5\n\n");
      bgNum[i]->Scale(5);
    }
    sigEff[i] = calculateEffAndErrors(sigNum[i], sigDen[i], TString::Format("sigEff_wp%d", i), binLimits);
    bgEff[i]  = calculateEffAndErrors(bgNum[i],  bgDen[i],  TString::Format("bgEff_wp%d", i),  binLimits);

    setHistogram(sigEff[i], i, true);
    setHistogram( bgEff[i], i, false);
    
    // Draw on canvas
    sigEff[i]->Draw("same,pe");
    if(mode != EFF_PT_2TEV) bgEff[i] ->Draw("same,pe");
    c1->Update();

  } // end loop over working points
  
  TLegend *leg = nullptr;
  if( mode == EFF_PT or mode == EFF_PT_2TEV ) leg = new TLegend(0.2, 0.3, 0.6, 0.7);
  else if( mode == EFF_ETA )                  leg = new TLegend(0.38, 0.2, 0.75, 0.6);
  else if( mode == EFF_NVTX )                 leg = new TLegend(0.2, 0.2, 0.6, 0.6);
  else{
    printf("ERROR: unknown mode requested\n");
    assert(0);
  }

  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry((TObject*)0, "Signal:", "");
  for(int i=0; i<Opt::nWP; i++){
    leg->AddEntry(sigEff[i], sigLegString[i], "pl");
  }
  if(mode != EFF_PT_2TEV){
    leg->AddEntry((TObject*)0, "", "");
    leg->AddEntry((TObject*)0, "Background:", "");
    for(int i=0; i<Opt::nWP; i++){
      leg->AddEntry(bgEff[i], bgLegString[i], "pl");
    }
  }
  leg->Draw("same");

  TLatex *lat = new TLatex(0.5, 0.95, comment); // 0.85
  lat->SetNDC(kTRUE);
  lat->Draw("same");

  TString filename = "figures/plot_eff_";
  if( mode != EFF_ETA ){
    if( drawBarrel ) filename += "barrel_";
    else             filename += "endcap_";
  }

  filename += varName[mode];
  if ( mode == EFF_PT_2TEV) filename += "_2TeV";
  filename += ".png";
  c1->Print(filename);
}

// Get a given tree from a given file name.
TTree *getTreeFromFile(TString fname, TString tname){
  TFile *file = new TFile( fname );
  TTree *tree = (TTree*) file->Get(tname);
  return tree;
}

TH1F* calculateEffAndErrors(TH1F *numHist, TH1F *denHist, TString eName, double* binLimits){
  int nBins = numHist->GetNbinsX();

  TH1F* ehist = new TH1F(eName, eName, nBins, binLimits);
  ehist->Sumw2();

  //ehist->Divide(denomHist);
  for(int i=1;i<=nBins;++i){
    //
    double eff = 0.0;
    double effErr2 = 0.0;
    double effErr = 0.0;
    //
    double num = numHist->GetBinContent(i);
    double numErr = numHist->GetBinError(i);
    double den = denHist->GetBinContent(i);
    double denErr = denHist->GetBinError(i);
    double tot = num+den;

    if( den != 0 ){
      eff = num/den;
      if( num != 0 )    effErr2 = (num*num*den*den)/(tot*tot*tot*tot) * ( numErr*numErr/(num*num) + denErr*denErr/(den*den) );
      if( effErr2 < 0 ) effErr2 = 0;
      effErr = sqrt(effErr2);
    }

    ehist->SetBinContent(i, eff);
    ehist->SetBinError(i, effErr);
  }
  return ehist;
}

void addBins(int &nBins, double step, int times){
  for(int i=0; i < times; ++i){
    ptBinLimits[nBins+1] = ptBinLimits[nBins] + step;
    ++nBins;
  }
}

double *getPtBinLimits(int &nBins, bool longPtRange){
  nBins = 0;
  ptBinLimits[0] = 20;
  addBins(nBins, 1,  20);
  addBins(nBins, 2,  30);
  addBins(nBins, 5,  10);
  addBins(nBins, 10,  5);
  if(longPtRange){
    addBins(nBins, 20,  20);
    addBins(nBins, 40,  10);
    addBins(nBins, 100, 10);
  }
  return ptBinLimits;
}

double *getEtaBinLimits(int &nBins){
  nBins = 50;
  double xmin = -2.5;
  double xmax = 2.5;
  double delta = (xmax - xmin)/nBins;
  
  for(int i=0; i<=nBins; i++) etaBinLimits[i] = xmin + i*delta;
  return etaBinLimits;
}

double *getNvtxBinLimits(int &nBins){
  nBins = nNvtxBins;
  return nvtxBinLimits;
}

void setHistogram(TH1F *hist, int wp, bool isSignal){
  int color = 0;
  if( isSignal ) color = sigColors[wp];
  else           color = bgColors[wp];

  hist->SetLineColor(color);
  hist->SetLineWidth(2);

  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(0.8);
  hist->SetMarkerColor(color);

}
