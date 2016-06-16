#include <iostream>
#include <fstream>
#include "TSystem.h"
#include "TStyle.h"
#include "TF1.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TFile.h"
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

#include "OptimizationConstants.hh"

enum MatchType  {MATCH_TRUE, MATCH_FAKE, MATCH_ANY};
enum SampleType {SAMPLE_UNDEF, SAMPLE_DY, SAMPLE_TT, SAMPLE_GJ};
enum EtaRegion  {ETA_EB, ETA_EE, ETA_FULL};

const MatchType  matchType      = MATCH_ANY; // MC truth matching to signal or bg electrons
const EtaRegion  etaRegion      = ETA_FULL;   // barrel, endcap, or all etas
const SampleType sample         = SAMPLE_GJ;

// Preselection cuts: must match or be looser than 
// cuts in OptimizatioConstants.hh
const float ptMin = 20;
const float etaMax = 2.5;
const float dzMax = 1.0;
const float boundaryEBEE = 1.479;

//====================================================
const bool talkativeRegime = false;
const bool smallEventCount = false;   
const int maxEventsSmall = 1000000;

// Files input
const TString fileNameDY = "~/DYJetsToLL_madgraph_80X_v3.root";
const TString fileNameTT = "~/TTJets_amcatnlo_80X_v3.root";
const TString fileNameGJ = "~/GJet_pythia8_80X_v1.root";
// Tree name input 
const TString treeName = "ntupler/ElectronTree";
// File and histogram with kinematic weights
const TString fileNameWeights = "kinematicWeights_20160611.root";
const TString histNameWeights = "hKinematicWeights";

// Files output
const TString flatNtupleFileNameBaseDY = "DYJetsToLL_jun14_flat_ntuple";
const TString flatNtupleFileNameBaseTT = "TTJets_jun14_flat_ntuple";
const TString flatNtupleFileNameBaseGJ = "GJet_jun14_flat_ntuple";


// //  Files IN 
// const TString fileNameS 
// = "/afs/cern.ch/user/i/ikrav/workspace/ntuples/Spring16/DYJetsToLL_madgraph_80X_v3.root";
// const TString fileNameBG 
// = "/afs/cern.ch/user/i/ikrav/workspace/ntuples/Spring16/TTJets_amcatnlo_80X_v3.root";
// // Tree Name (file IN):

// //  Files OUT

// Effective areas for electrons derived by Ilya for Summer16
// https://indico.cern.ch/event/482673/contributions/2187022/attachments/1282446/1905912/talk_electron_ID_spring16.pdf
namespace EffectiveAreas {
  const int nEtaBins = 7;
  const float etaBinLimits[nEtaBins+1] = {
    0.0, 1.0, 1.479, 2.0, 
    2.2, 2.3, 2.4, 2.5
  };
  const float effectiveAreaValues[nEtaBins] = { 
    0.1703,	 0.1715,	 0.1213,	 
    0.1230,	 0.1635,	 0.1937,	 
    0.2393	
  };
}


void bazinga (std::string mes){
  if (talkativeRegime)
    std::cout<<"\n"<<mes<<endl;
}

// Forward declarations
float findKinematicWeight(TH2D *hist, float pt, float etaSC);
bool passPreselection(int isTrue, float pt, float eta, 
		      int passConversionVeto, float dz);
TString eventCountString();
void drawProgressBar(float progress);

//
// Main program
//

void convert_EventStrNtuple_To_FlatNtuple(){

  bazinga("Start main function");
  gBenchmark->Start("Timing");

  //  Generate  a dictionary, so that CINT will have all the information 
  // it needs about type or variable at anytime.
  gROOT->ProcessLine("#include <vector>"); 
  
  // General settings in case one adds drawing of histograms etc
  gStyle->SetOptFit();
  gStyle->SetOptStat(0);

  // ======================================================================
  // Set up input/output files and find/create trees
  // ======================================================================
  bazinga("Set up input/output files\n");

  // Input/output file names
  TString inputFileName = "";
  TString flatNtupleFileNameBase = "Undefined";
  if( sample == SAMPLE_DY ){
    inputFileName = fileNameDY;
    flatNtupleFileNameBase = flatNtupleFileNameBaseDY;
  }else if( sample == SAMPLE_TT ){
    inputFileName = fileNameTT;
    flatNtupleFileNameBase = flatNtupleFileNameBaseTT;
  }else if( sample == SAMPLE_GJ ){
    inputFileName = fileNameGJ;
    flatNtupleFileNameBase = flatNtupleFileNameBaseGJ;
  }else{
    printf("Unknown sample requested\n");
    assert(0);
  }

  TString flatNtupleFileNameTruth = "";
  if( matchType == MATCH_TRUE ){
    flatNtupleFileNameTruth = "_true";
  }else if( matchType == MATCH_FAKE ){
    flatNtupleFileNameTruth = "_fake";
  }else if( matchType == MATCH_ANY ){
    flatNtupleFileNameTruth = "_trueAndFake";
  }

  TString flatNtupleFileNameEtas = "";
  if( etaRegion == ETA_EB ){
    flatNtupleFileNameEtas = "_barrel";
  }else if( etaRegion == ETA_EE ){
    flatNtupleFileNameEtas = "_endcap";
  }else if( etaRegion == ETA_FULL ){
    flatNtupleFileNameEtas = "_alleta";
  }

  TString flatNtupleFileNameEvents = eventCountString();
  TString flatNtupleFileNameEnding = ".root";

  TString flatNtupleFileName = flatNtupleFileNameBase + flatNtupleFileNameTruth 
    + flatNtupleFileNameEtas + flatNtupleFileNameEvents + flatNtupleFileNameEnding;

  // Open input file and find the tree
  TFile *inputFile = new TFile(inputFileName);
  if( !inputFile ){
    printf("Failed to open input file %s\n", inputFileName.Data());
    assert(0);
  }
  TTree *treeIn = (TTree*)inputFile->Get(treeName);
  if( !treeIn ){
    printf("Failed to find tree %s in the file %s\n", treeName.Data(), inputFileName.Data());
    assert(0);
  }

  // Open output file
  TFile *fileOut = new TFile(flatNtupleFileName, "recreate");
  TTree *treeOut = new TTree("electronTree","Flat_ntuple");

  // ======================================================================
  // Set up all variables and branches for the input tree
  // ======================================================================
  //
  // Declare variables and branches
  // Event-level variables:
  int eleNEle; // the number of reconstructed electrons in the event
  float eleRho;
  float genWeight;
  int nPV;
  // Per-eletron variables
  // Kinematics
  std::vector <float> *elePt = 0;         // electron PT
  std::vector <float> *eleEtaSC = 0;      // supercluser eta
  std::vector <float> *elePhiSC = 0;      // supercluser phi
  // Variables for analysis
  std::vector <float> *isoChargedHadrons = 0;
  std::vector <float> *isoNeutralHadrons = 0;
  std::vector <float> *eleIsoPhotons = 0;
  std::vector <int> *eleIsTrueElectron = 0;
  // Other vars  
  // Impact parameters
  std::vector <float> *eleD0 = 0;      // r-phi plane impact parameter
  std::vector <float> *eleDZ = 0;      // r-z plane impact parameter
  // Matching track-supercluster
  std::vector <float> *eleDEtaSeed = 0;  // deltaEtaIn
  std::vector <float> *eleDPhiIn = 0;  // deltaPhiIn
  // Misc ID variables
  std::vector <float> *eleHoverE = 0;  // H/E  
  std::vector <float> *eleFull5x5SigmaIEtaIEta = 0;  
  std::vector <float> *eleOOEMOOP = 0; // |1/E - 1/p|
  // Conversion rejection
  std::vector <int> *eleExpectedMissingInnerHits = 0;
  std::vector <int> *electronPassConversionVeto = 0;

  // Declare branches
  TBranch *b_eleNEle = 0;
  TBranch *b_nPV = 0;
  TBranch *b_genWeight = 0;
  TBranch *b_eleRho = 0;
  TBranch *b_elePt = 0;
  TBranch *b_eleEtaSC = 0;
  TBranch *b_elePhiSC = 0;
  TBranch *b_isoChargedHadrons = 0;
  TBranch *b_isoNeutralHadrons = 0;
  TBranch *b_eleIsoPhotons = 0;
  TBranch *b_eleIsTrueElectron;
  
  // Other vars
  TBranch *b_eleD0 = 0;
  TBranch *b_eleDZ = 0;
  TBranch *b_eleDEtaSeed = 0;
  TBranch *b_eleDPhiIn = 0;
  TBranch *b_eleHoverE = 0;
  TBranch *b_eleFull5x5SigmaIEtaIEta = 0;
  TBranch *b_eleOOEMOOP = 0;
  TBranch *b_eleExpectedMissingInnerHits = 0;
  TBranch *b_electronPassConversionVeto = 0;

  // Connect variables and branches to the tree with the data
  treeIn->SetBranchAddress("nEle", &eleNEle, &b_eleNEle);
  treeIn->SetBranchAddress("nPV", &nPV, &b_nPV);
  treeIn->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
  treeIn->SetBranchAddress("rho", &eleRho, &b_eleRho);
  treeIn->SetBranchAddress("pt", &elePt, &b_elePt);
  treeIn->SetBranchAddress("etaSC", &eleEtaSC, &b_eleEtaSC);
  treeIn->SetBranchAddress("phiSC", &elePhiSC, &b_elePhiSC);
  treeIn->SetBranchAddress("isoChargedHadrons", &isoChargedHadrons, &b_isoChargedHadrons);
  treeIn->SetBranchAddress("isoNeutralHadrons", &isoNeutralHadrons, &b_isoNeutralHadrons);
  treeIn->SetBranchAddress("isoPhotons",        &eleIsoPhotons,        &b_eleIsoPhotons);
  treeIn->SetBranchAddress("isTrue",    &eleIsTrueElectron,    &b_eleIsTrueElectron);
  treeIn->SetBranchAddress("d0",                &eleD0,             &b_eleD0);
  treeIn->SetBranchAddress("dz",                &eleDZ,             &b_eleDZ);
  treeIn->SetBranchAddress("dEtaSeed",            &eleDEtaSeed,         &b_eleDEtaSeed);
  treeIn->SetBranchAddress("dPhiIn",            &eleDPhiIn,         &b_eleDPhiIn);
  treeIn->SetBranchAddress("hOverE",            &eleHoverE,         &b_eleHoverE);
  treeIn->SetBranchAddress("full5x5_sigmaIetaIeta", &eleFull5x5SigmaIEtaIEta,
			  &b_eleFull5x5SigmaIEtaIEta);
  treeIn->SetBranchAddress("ooEmooP",           &eleOOEMOOP,        &b_eleOOEMOOP);
  treeIn->SetBranchAddress("expectedMissingInnerHits", &eleExpectedMissingInnerHits, 
			  &b_eleExpectedMissingInnerHits);
  treeIn->SetBranchAddress("passConversionVeto",       &electronPassConversionVeto,
			  &b_electronPassConversionVeto);

  
  //
  // Get the histogram with kinematic weights
  //  
  TFile *fweights = new TFile(fileNameWeights);
  if( !fweights ){
    printf("File with kinematic weights %s is not found\n", fileNameWeights.Data());
    assert(0);
  }
  TH2D *hKinematicWeights = (TH2D*)fweights->Get("hKinematicWeights");
  if( !hKinematicWeights ){
    printf("The histogram %s is not found in file %s\n", 
	   histNameWeights.Data(), fileNameWeights.Data());
    assert(0);
  }
  
  // ======================================================================
  // Set up all variables and branches for the input tree
  // ======================================================================

  //bazinga("Set up output tree\n");
  fileOut->cd();  //?????????
  
  Float_t gweight_ = 0.0;      // gen Weight
  Float_t kweight_ = 0.0;      // kinematic Weight

  // all electron variables
  Float_t pt_ = 0.0;
  Float_t etaSC_= 0.0;
  Float_t dEtaSeed_= 0.0;
  Float_t dPhiIn_= 0.0;
  Float_t hOverE_= 0.0;
  Float_t full5x5_sigmaIetaIeta_= 0.0;
  Float_t isoChargedHadrons_= 0.0;
  Float_t isoNeutralHadrons_= 0.0;
  Float_t isoPhotons_= 0.0;
  Float_t isoChargedFromPU_= 0.0;
  Float_t relIsoWithEA_= 0.0;
  Float_t relIsoWithDBeta_= 0.0;
  Float_t ooEmooP_= 0.0;
  Float_t d0_= 0.0;
  Float_t dz_= 0.0;
  Int_t   expectedMissingInnerHits_= 0;
  Int_t   passConversionVeto_= 0;
  Int_t   isTrueEle_= 0;
    
  treeOut->Branch("genWeight"  ,  &gweight_ , "gweight/F");  
  treeOut->Branch("kinWeight"  ,  &kweight_ , "kweight/F");  

  treeOut->Branch("pt"    ,  &pt_    , "pt/F");			    
  treeOut->Branch("etaSC" ,  &etaSC_ , "etaSC/F");
  treeOut->Branch("dEtaSeed",  &dEtaSeed_, "dEtaSeed/F");
  treeOut->Branch("dPhiIn",  &dPhiIn_, "dPhiIn/F");
  treeOut->Branch("hOverE",  &hOverE_, "hOverE/F");
  treeOut->Branch("full5x5_sigmaIetaIeta", &full5x5_sigmaIetaIeta_, "full5x5_sigmaIetaIeta/F");
  treeOut->Branch("relIsoWithEA"           , &relIsoWithEA_, "relIsoWithEA/F");
  treeOut->Branch("ooEmooP", &ooEmooP_, "ooEmooP/F");
  treeOut->Branch("d0"     , &d0_,      "d0/F");
  treeOut->Branch("dz"     , &dz_,      "dz/F");
  treeOut->Branch("expectedMissingInnerHits", &expectedMissingInnerHits_, "expectedMissingInnerHits/I");
  treeOut->Branch("passConversionVeto", &passConversionVeto_, "passConversionVeto/I");
  treeOut->Branch("isTrueEle"    , &isTrueEle_,     "isTrueEle/I");
  
  //
  // Loop over events
  //
  UInt_t maxEvents = treeIn->GetEntries();
  UInt_t maxEventsOver10000 =   maxEvents/10000.;

  if( smallEventCount )
    maxEvents = maxEventsSmall;

  printf("\nStart processing events, will run on %u events\n", maxEvents );

  for(UInt_t ievent = 0; ievent < maxEvents; ievent++){

    Long64_t tentry = treeIn->LoadTree(ievent);
    // Load the value of the number of the electrons in the event    
    b_eleNEle->GetEntry(tentry);
    
    if( ievent%100000 == 0 || ievent == maxEvents-1){
      //printf("."); fflush(stdout);
      drawProgressBar( (1.0*ievent+1)/maxEvents);
    }
      
    // Get data for all electrons in this event, only vars of interest
    b_eleRho->GetEntry(tentry);
    b_genWeight->GetEntry(tentry);
    b_elePt->GetEntry(tentry);
    b_eleEtaSC->GetEntry(tentry);
    b_elePhiSC->GetEntry(tentry);
    b_isoChargedHadrons->GetEntry(tentry);
    b_isoNeutralHadrons->GetEntry(tentry);
    b_eleIsoPhotons->GetEntry(tentry);
    b_eleIsTrueElectron->GetEntry(tentry);

    // // Other vars
    b_eleD0->GetEntry(tentry);
    b_nPV->GetEntry(tentry);
    b_eleDZ->GetEntry(tentry);
    b_eleDEtaSeed->GetEntry(tentry);
    b_eleDPhiIn->GetEntry(tentry);
    b_eleHoverE->GetEntry(tentry);
    b_eleFull5x5SigmaIEtaIEta->GetEntry(tentry);
    b_eleOOEMOOP->GetEntry(tentry);
    b_eleExpectedMissingInnerHits->GetEntry(tentry);
    b_electronPassConversionVeto->GetEntry(tentry);

    // Loop over the electrons
    for(int iele = 0; iele < eleNEle; iele++){
     
      gweight_ = genWeight;
      // nEle_ = eleNEle;
      pt_ = elePt->at(iele);
      etaSC_ = eleEtaSC->at(iele);
      kweight_ = findKinematicWeight(hKinematicWeights, pt_, etaSC_);
      
      dEtaSeed_ = eleDEtaSeed->at(iele);
      dPhiIn_ = eleDPhiIn->at(iele);
      full5x5_sigmaIetaIeta_ = eleFull5x5SigmaIEtaIEta->at(iele);
      hOverE_ =  eleHoverE->at(iele);
      d0_ = eleD0->at(iele);
      dz_ = eleDZ->at(iele);
      isoChargedHadrons_ =  isoChargedHadrons->at(iele);
      isoNeutralHadrons_ = isoNeutralHadrons->at(iele);
      isoPhotons_ = eleIsoPhotons->at(iele);
      expectedMissingInnerHits_ = eleExpectedMissingInnerHits->at(iele);
      // nPV_ = nPV;
      // // nPU_ = nPU->at(iele);
      // // nPUTrue_ = nPUTrue->at(iele);
      // rho_ = eleRho; //->at(iele);
      ooEmooP_ = eleOOEMOOP->at(iele);
      passConversionVeto_ = electronPassConversionVeto->at(iele);
      isTrueEle_ =  eleIsTrueElectron->at(iele); 

      // Compute isolation with effective area correction for PU
      // Find eta bin first. If eta>2.5, the last eta bin is used.
      int etaBin = 0; 
      while ( etaBin < EffectiveAreas::nEtaBins-1 
	      && abs(etaSC_) > EffectiveAreas::etaBinLimits[etaBin+1] )
	{ ++etaBin; };
      double area = EffectiveAreas::effectiveAreaValues[etaBin];
      relIsoWithEA_ = (  isoChargedHadrons_
			 + max(0.0, isoNeutralHadrons_ + isoPhotons_ 
			       - eleRho * area ) )/pt_;
      if( ! passPreselection( isTrueEle_, pt_, etaSC_, passConversionVeto_, dz_) )
	  continue;
      treeOut->Fill();// IK will kill me next time, if this line is not at the right place!
      
    } // end loop over the electrons
  
  } // end loop over SIGNAL events
  
  bazinga("I'm here to write signal tree");    
  treeOut->Write();
  fileOut->Write();
  fileOut->Close();
  delete fileOut; 
  fileOut = nullptr;
  delete inputFile;
  inputFile = nullptr;
  
  bazinga("I'm finished with calculation, here is the info from dBenchmark:");  
  printf("\n");  gBenchmark->Show("Timing");  // get timing info
} // end of main fcn

float findKinematicWeight(TH2D *hist, float pt, float etaSC){
  
  // For signal electrons from Drell-Yan, use kinematic weights,
  // for background electrons from ttbar, do not use any weights
  float weight = 1.0;
  if( matchType == MATCH_TRUE ){
    // Retrieve kinematic weight
    int ipt = hist->GetXaxis()->FindBin(pt);
    int npt = hist->GetNbinsX();
    int ieta = hist->GetYaxis()->FindBin(etaSC);
    int neta = hist->GetNbinsY();
    if( ipt < 1 || ipt > npt || ieta < 1 || ieta > neta ){      
      // If pt and eta are outside of the limits of the weight histogram,
      // set the weight to the edge value
      if( ipt < 1 ) 
  	ipt = 1;
      if( ipt > npt )
  	ipt = npt;
      if( ieta < 1 ) 
  	ieta = 1;
      if( ieta > neta ) 
  	ieta = neta;
      weight = hist->GetBinContent(ipt, ieta);
    } else {
      // Normal case, pt and eta are within limits of the weight histogram
      weight = hist->Interpolate(pt, etaSC);
    } // end is within hist limits
  } // end is signal 

  return weight;
}

bool passPreselection(int isTrue, float pt, float eta, 
		      int passConversionVeto, float dz){

  bool pass = true;

  if( matchType == MATCH_TRUE ){
    if( !(isTrue==1) )
      pass = false;
  }else if( matchType == MATCH_FAKE ){
    if( !(isTrue==0 || isTrue==3) )
      pass = false;
  }else if( matchType == MATCH_ANY){
    // No matching needed, do nothing
  }else{
    printf("Unknown truth match requested\n");
    assert(0);
  }
    
  if( pt < ptMin )
    pass = false;

  if( etaRegion == ETA_EB ){
    if( !(abs(eta) <= boundaryEBEE ) )
      pass = false;
  }else if( etaRegion == ETA_EE ){
    if( !(abs(eta) >= boundaryEBEE  && abs(eta) <= etaMax ) )
      pass = false;
  }else if( etaRegion == ETA_FULL ){
    if( !( abs(eta) < etaMax ) )
      pass = false;
  }else{
    printf("Unknown eta region requested\n");
    assert(0);
  }

  if( !passConversionVeto )
    pass = false;
  
  if( ! (abs(dz) < dzMax ) )
    pass = false;
  
  return pass;
  
}

const int maxPowersBin = 6;
TString powers[maxPowersBin] = {"","K","M","G","T","P"};

TString eventCountString(){

  TString result = "_full";

  if( smallEventCount ){
    int powerOfTen = log10(maxEventsSmall);
    int powersBin = powerOfTen/3;
    if( powersBin >= maxPowersBin )
      powersBin = maxPowersBin-1;
    int neventsShort = maxEventsSmall / TMath::Power(10, 3*powersBin);
    result = TString::Format("_%d%s", neventsShort, powers[powersBin].Data());
  }

  return result;
}

void drawProgressBar(float progress){

  const int barWidth = 70;
  
  std::cout << "[";
  int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos) std::cout << "=";
    else if (i == pos) std::cout << ">";
    else std::cout << " ";
  }
  std::cout << "] " << int(progress * 100.0) << " %\r";
  std::cout.flush();
  
  if( progress >= 1.0 )
    std::cout << std::endl;

  return;
}
