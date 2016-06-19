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

const bool doOverlayCuts = true;

// File names of the files that contain the working points to display
// (of doOverlayCuts above is false, no cuts will be shown, and the 
// content of this array is ignored).
const TString cutFileNamesBarrel[4] = { 
  // only WP Veto was optimized with the exercise purpose
  "./cut_repository/cuts_barrel_pass1_20160611_200000_WP_Veto.root",
  "./cut_repository/cuts_barrel_pass2_20160611_200000_WP_Loose.root",
  "./cut_repository/cuts_barrel_HLTSafe_20160613_200000.root",
  "./cut_repository/cuts_barrel_pass2_20160611_200000_WP_Loose.root"
};
const TString cutFileNamesEndcap[4] = {
  // for aestetics purpose, switch to real endcap files, once they are made
  "./cut_repository/cuts_barrel_pass1_20160611_200000_WP_Veto.root",
  "./cut_repository/cuts_barrel_pass2_20160611_200000_WP_Loose.root",
  "./cut_repository/cuts_endcap_HLTSafe_20160613_200000.root",
  "./cut_repository/cuts_barrel_pass2_20160611_200000_WP_Loose.root"
};

static Int_t c_Canvas         = TColor::GetColor( "#f0f0f0" );
static Int_t c_FrameFill      = TColor::GetColor( "#fffffd" );
static Int_t c_TitleBox       = TColor::GetColor( "#5D6B7D" );
static Int_t c_TitleBorder    = TColor::GetColor( "#7D8B9D" );
static Int_t c_TitleText      = TColor::GetColor( "#FFFFFF" );
static Int_t c_SignalLine     = TColor::GetColor( "#0000ee" );
static Int_t c_SignalFill     = TColor::GetColor( "#7d99d1" );
static Int_t c_BackgroundLine = TColor::GetColor( "#ff0000" );
static Int_t c_BackgroundFill = TColor::GetColor( "#ff0000" );
static Int_t c_NovelBlue      = TColor::GetColor( "#2244a5" );

// Forward declaraitons
TCanvas * drawOneVariable(TTree *signalTree, 
			  TTree *backgroundTree1, TTree *backgroundTree2,TTree *backgroundTreeAdditional,
			  TCut signalCuts, TCut backgroundCuts,
			  TString var,
			  int nbins, double xlow, double xhigh,
			  TString sigLegend, TString bg1Legend, TString bg2Legend,TString bg3Legend,
			  TString comment);

void overlayCuts(TCanvas *canvas, TString variable);

void setHistogramAttributes(TH1F *hsig, TH1F *hbg1, TH1F *hbg2, TH1F *hbg3);

// Define a helper struct
struct VarPlotSettings {
  TString varName;
  int nbins;
  float xmin;
  float xmax;
  VarPlotSettings(TString nameIn, int nbinsIn, float xminIn, float xmaxIn) :
    varName(nameIn), nbins(nbinsIn), xmin(xminIn), xmax(xmaxIn) {};
};

//
// Main function
//
void drawVariablesAndCuts_HLTsafe(){
  //
  // Open files
  //
  TString fname1 = "DYJetsToLL_jun14_flat_ntuple_trueAndFake_barrel_full.root";
  if( !drawBarrel )
    fname1 = "DYJetsToLL_jun14_flat_ntuple_trueAndFake_endcap_full.root";
  TFile *input1 = new TFile( fname1 );
  //  input1->cd("wp3");
  TTree *signalTree     = (TTree*) input1->Get("electronTree");
  
  TString fname2 = "TTJets_may29_flat_ntuple_full_barrel.root";
  if( !drawBarrel )
    fname2 = "TTJets_may29_flat_ntuple_full_endcap.root";
  TFile *input2 = new TFile( fname2 );
  //  input2->cd("wp3");
  TTree *backgroundTree     = (TTree*) input2->Get("electronTree");
  
  TString fname3 = "GJet_jun14_flat_ntuple_trueAndFake_alleta_full.root";
  TFile *input3 = new TFile( fname3 );
  //  input3->cd("wp3");
  TTree *backgroundTreeAdditional     = 0;
  if( input3 )
    backgroundTreeAdditional = (TTree*) input3->Get("electronTree");
 
 //
  // Define cuts
  //
  TString comment = "barrel electrons";
  TCut preselectionCuts = Opt::ptCut;
  if( drawBarrel ){
    comment = "barrel electrons";
    preselectionCuts += Opt::etaCutBarrel;
  }else{
    comment = "endcap electrons";
    preselectionCuts += Opt::etaCutEndcap;
  }
  preselectionCuts += Opt::otherPreselectionCuts;

  TCut signalCuts     = preselectionCuts && Opt::trueEleCut;
  TCut backgroundCuts = preselectionCuts && Opt::fakeEleCut;

  // 
  // Define details of the plots: names, limits.
  //
  vector<struct VarPlotSettings*> vPlotSettings;
  // for the purpose of the exercise plot only sigma_IeIe
  //
  // The constructor below takes: (<variable name>, <nbins>, <var_min>, <var max for plotting>)

 
  // vPlotSettings.push_back( new VarPlotSettings("full5x5_sigmaIetaIeta", 100, 0, drawBarrel ? 0.018 : 0.06));

  vPlotSettings.push_back( new VarPlotSettings("full5x5_sigmaIetaIeta", 100, 0.004, drawBarrel ? 0.018 : 0.06));
  vPlotSettings.push_back( new VarPlotSettings("dEtaSeed", 100, drawBarrel ? -0.0125 : -0.03 , drawBarrel ? 0.0125 : 0.03 ));
  vPlotSettings.push_back( new VarPlotSettings("dPhiIn", 100, drawBarrel ? -0.1 : -0.4, drawBarrel ? 0.1 : 0.4 ));
  vPlotSettings.push_back( new VarPlotSettings("hOverE", 100, 0.0, 0.1 ));
  vPlotSettings.push_back( new VarPlotSettings("relIsoWithEA", 100, 0, drawBarrel ? 0.25 : 0.3 ));
  vPlotSettings.push_back( new VarPlotSettings("ooEmooP", 100, 0.0, 0.15));
  /*
    vPlotSettings.push_back( new VarPlotSettings("d0", 100, -0.1, 0.1));
    vPlotSettings.push_back( new VarPlotSettings("dz", 100, -0.9, 0.9));
    
    vPlotSettings.push_back( new VarPlotSettings("expectedMissingInnerHits", 5, -0.5, 4.5));
  */
  TCanvas *c1;
  TString variable = "";
  float xmin, xmax;
  int nbins;

  for( uint i=0; i<vPlotSettings.size(); i++){
    variable = vPlotSettings.at(i)->varName;
    nbins = vPlotSettings.at(i)->nbins;
    xmin = vPlotSettings.at(i)->xmin;
    xmax = vPlotSettings.at(i)->xmax;
    c1 = drawOneVariable(signalTree, signalTree, backgroundTree,backgroundTreeAdditional,
			 signalCuts, backgroundCuts,
			 variable, nbins, xmin, xmax,
			 "signal DYJetsToLL", "fakes from DYJetsToLL", "fakes from TTJets","fakes from GJets" ,comment);
    if( doOverlayCuts )
      overlayCuts(c1, variable);

    TString outname = "figures/plot_barrel_HLTSafe_";
    if( !drawBarrel )
      outname = "figures/plot_endcap_HLTSafe_";
    outname += variable;
    outname += ".png";
    c1->Print(outname);
  }


  
}

TCanvas *drawOneVariable(TTree *signalTree, TTree *backgroundTree1, TTree *backgroundTree2,TTree *backgroundTree3,
			 TCut signalCuts, TCut backgroundCuts,
			 TString var,
			 int nbins, double xlow, double xhigh,
			 TString sigLegend, TString bg1Legend, TString bg2Legend,TString bg3Legend,
			 TString comment)
{

  TString cname = "c_";
  cname += var;
  TCanvas *c1 = new TCanvas(cname,cname,10,10,600,600);
  c1->cd();

  TH1F *hsig = new TH1F(TString("hsig_")+var,"",nbins, xlow, xhigh);
  TH1F *hbg1 = new TH1F(TString("hbg1_")+var,"",nbins, xlow, xhigh);
  TH1F *hbg2 = new TH1F(TString("hbg2_")+var,"",nbins, xlow, xhigh);
  TH1F *hbg3 = new TH1F(TString("hbg3_")+var,"",nbins, xlow, xhigh);

  TString projectCommandSig = var+TString(">>hsig_")+var;
  TString projectCommandBg1 = var+TString(">>hbg1_")+var;
  TString projectCommandBg2 = var+TString(">>hbg2_")+var;
  TString projectCommandBg3 = var+TString(">>hbg3_")+var;
  
  if( !useSmallEventCount ){
    signalTree->Draw(projectCommandSig, "genWeight*kinWeight"*signalCuts);
  }else{
    printf("DEBUG MODE: using small event count\n");
    signalTree->Draw(projectCommandSig, "genWeight*kinWeight"*signalCuts, "", 100000);
  }
  TGaxis::SetMaxDigits(3);
  hsig->GetXaxis()->SetTitle(var);
  hsig->SetDirectory(0);

  if(backgroundTree1 != 0){
    if( !useSmallEventCount ){
      backgroundTree1->Draw(projectCommandBg1, "genWeight*kinWeight"*backgroundCuts);
    }else{
      printf("DEBUG MODE: using small event count\n");
      backgroundTree1->Draw(projectCommandBg1, "genWeight*kinWeight"*backgroundCuts, "", 100000);
    }
    hbg1->Scale(hsig->GetSumOfWeights() / hbg1->GetSumOfWeights());
    hbg1->SetDirectory(0);
  } else {
    delete hbg1;
    hbg1 = 0;
  }

  if(backgroundTree2 != 0){
    if( !useSmallEventCount ){
      backgroundTree2->Draw(projectCommandBg2, "genWeight*kinWeight"*backgroundCuts);
    }else{
      printf("DEBUG MODE: using small event count\n");
      backgroundTree2->Draw(projectCommandBg2, "genWeight*kinWeight"*backgroundCuts, "", 100000);
    }
    hbg2->Scale(hsig->GetSumOfWeights() / hbg2->GetSumOfWeights());
    hbg2->SetDirectory(0);
  } else {
    delete hbg2;
    hbg2 = 0;
  }


  if(backgroundTree3 != 0){
    if( !useSmallEventCount ){
      backgroundTree3->Draw(projectCommandBg3, "genWeight*kinWeight"*backgroundCuts);
    }else{
      printf("DEBUG MODE: using small event count\n");
      backgroundTree3->Draw(projectCommandBg3, "genWeight*kinWeight"*backgroundCuts, "", 100000);
    }
    hbg3->Scale(hsig->GetSumOfWeights() / hbg3->GetSumOfWeights());
    hbg3->SetDirectory(0);
  } else {
    delete hbg3;
    hbg3 = 0;
  }


  setHistogramAttributes(hsig, hbg1, hbg2, hbg3);

  c1->Clear();

  hsig->Draw("hist");
  if( hbg1 ){

    hbg1->Draw("same");
  }
  if( hbg2 ){

    hbg2->Draw("same");
  }
  if( hbg3 ){

    hbg3->Draw("same");
  }
  TLegend *leg = new TLegend(0.55, 0.65, 0.95, 0.80); // 0.6 0.9
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hsig, sigLegend, "lf");
  leg->AddEntry(hbg1, bg1Legend, "lf");
  leg->AddEntry(hbg2, bg2Legend, "lf");
  leg->AddEntry(hbg3, bg3Legend, "lf");
  leg->Draw("same");

  TLatex *lat = new TLatex(0.5, 0.95, comment); // 0.85
  lat->SetNDC(kTRUE);
  lat->Draw("same");

  c1->Update();

  return c1;
}

void setHistogramAttributes(TH1F *hsig, TH1F *hbg1, TH1F *hbg2, TH1F * hbg3){

  //signal
  // const Int_t FillColor__S = 38 + 150; // change of Color Scheme in ROOT-5.16.
  // convince yourself with gROOT->GetListOfColors()->Print()
  Int_t FillColor__S = c_SignalFill;
  Int_t FillStyle__S = 1001;
  Int_t LineColor__S = c_SignalLine;
  Int_t LineWidth__S = 2;
  
  // background
  //Int_t icolor = UsePaperStyle ? 2 + 100 : 2;
  Int_t FillColor__B = c_BackgroundFill;
  Int_t FillStyle__B = 3554;
  Int_t LineColor__B = c_BackgroundLine;
  Int_t LineWidth__B = 2;
  
  if (hsig != NULL) {
    hsig->SetStats(0);
    hsig->SetLineColor( LineColor__S );
    hsig->SetLineWidth( LineWidth__S );
    hsig->SetFillStyle( FillStyle__S );
    hsig->SetFillColor( FillColor__S );
  }
  
  if (hbg1 != NULL) {
    hbg1->SetStats(0);
    hbg1->SetLineColor( LineColor__B );
    hbg1->SetLineWidth( LineWidth__B );
    hbg1->SetFillStyle( FillStyle__B );
    hbg1->SetFillColor( FillColor__B );
  }
  
  if (hbg2 != NULL) {
    hbg2->SetStats(0);
    hbg2->SetLineColor( kGreen+3 );
    hbg2->SetLineWidth( LineWidth__B );
    hbg2->SetFillStyle( 3003 );
    hbg2->SetFillColor( kGreen+3 );
  }
 
  if (hbg3 != NULL) {
    hbg3->SetStats(0);
    hbg3->SetLineColor( kMagenta+2 );
    hbg3->SetLineWidth( LineWidth__B );
    hbg3->SetFillStyle( 3018 );
    hbg3->SetFillColor( kMagenta+2 );
  }



  
}

void overlayCuts(TCanvas *canvas, TString variable){

  canvas->cd();
  double ymax = canvas->GetPad(0)->GetUymax();
  // The min/max below refer to the full span of the pad.
  // Normally, 10% of each side is taken by margins, so 
  // recompute that into in-plot limits
  double xminPad = canvas->GetPad(0)->GetX1();
  double xmaxPad = canvas->GetPad(0)->GetX2();
  double windowPad = xmaxPad - xminPad;
  double xmin = xminPad + windowPad*0.1;
  double xmax = xmaxPad - windowPad*0.1;

  // 
  // Read cut values for each working point from the 
  // corresponding file.
  //
  float cutVal[Opt::nWP];
  bool symmetric = false;
  // Only four working points are listed in the array of file names 
  // in the beginning of this file.
  if( Opt::nWP != 4 )
    assert(0);
  // Pull up the cuts
  for(int iwp=0; iwp<Opt::nWP-1; iwp++){ // DEBUG
    TString cutFileName = "";
    if( drawBarrel )
      cutFileName = cutFileNamesBarrel[iwp];
    else
      cutFileName = cutFileNamesEndcap[iwp];
    TFile *fcut = new TFile(cutFileName);
    if( !fcut )
      assert(0);
    VarCut *thisCut = (VarCut*)fcut->Get("cuts");
    if( !thisCut )
      assert(0);
    cutVal[iwp] = thisCut->getCutValue(variable);
    if (variable == "expectedMissingInnerHits") 
      cutVal[iwp] = int(cutVal[iwp]);
    // This flag below is overwritten every iteration, but that's ok
    // since all cuts of the same set should be the same wrt this.
    symmetric = thisCut->isSymmetric(variable);
    fcut->Close();
  }

  // 
  // Draw the cut values as arrows
  //
  TArrow *arrowUpper[Opt::nWP];
  TArrow *arrowLower[Opt::nWP];
  float arXmin, arXmax, arYmin, arYmax;
  for(int i=0; i<Opt::nWP-1; i++){ // DEBUG
    arXmin = arXmax = cutVal[i];
    arYmin = ymax*0.35 + ymax*i*0.05;
    arYmax = 0;
    if( cutVal[i] > xmax ){
      // Make the error horizontal instead
      float xwindow = xmax-xmin;
      arXmin = xmin + (0.92-0.04*i)*xwindow;
      arXmax = xmax;
      arYmin = arYmax = ymax*(0.50 + i*0.05);
    }
    arrowUpper[i] = new TArrow(arXmin, arYmin, arXmax, arYmax, 0.05, "|->");
    arrowUpper[i]->SetLineWidth(2);
    arrowUpper[i]->SetLineColor(kOrange+i);
    arrowUpper[i]->Draw();
    if(symmetric){
      arXmin = arXmax = -cutVal[i];
      arYmin = ymax*0.35 + ymax*i*0.05;
      arYmax = 0;
      if( -cutVal[i] < xmin ){
	// Make the error horizontal instead
	float xwindow = xmax-xmin;
	arXmin = xmin + (0.08+0.04*i)*xwindow;
	arXmax = xmin;
	arYmin = arYmax = ymax*(0.50 + i*0.05);
      }
      arrowLower[i] = new TArrow(arXmin, arYmin, arXmax, arYmax, 0.05, "|->");
      arrowLower[i]->SetLineWidth(2);
      arrowLower[i]->SetLineColor(kOrange+i);
      arrowLower[i]->Draw();
    }
  }

  canvas->Update();
}




