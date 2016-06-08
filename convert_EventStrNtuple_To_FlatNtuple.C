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


bool isDY = false;


//====================================================
const bool talkativeRegime = false;
const bool smallEventCount = true;   //5000k
const float etaMin = -2.5;
const float etaMax = +2.5;

const int ptMin = 20;
const int ptMax = 170;
const int etaBins = 120;
const int ptBins = 120;

const TCut trueEleCut = "isTrue == 1";
const TCut fakeEleCut = "isTrue == 0 || isTrue == 3";

const  TCut psCut = "pt>=20 && pt<=170 && abs(etaSC) < 2.5";
const  TCut gapPScut = "abs(etaSC)>1.4442 && abs(etaSC)<1.566";

const  TCut passConv_dZ = "passConversionVeto && abs(dZ)<1";
const  TCut preselectionS = trueEleCut && psCut && !gapPScut && passConv_dZ;
const  TCut preselectionBG = fakeEleCut && psCut && !gapPScut && passConv_dZ;
 
//  Files IN 
const TString fileNameS = "~/DYJetsToLL_madgraph_80X.root";
const TString fileNameBG = "~/TTJets_amcatnlo_80X.root";
// Tree Name (file IN):
const TString treeName = "ntupler/ElectronTree";
//  Files OUT
const TString flatNtupleNameS = "DYJetsToLL_may29_flat_ntuple_withWeights_5M.root";
const TString flatNtupleNameBG = "TTJets_may29_flat_ntuple_withWeights_5M.root";

// Effective areas for electrons derived by Rami and ILya for 25ns Spring16 DY sample in May 2016
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

//
// Main program
//

void convert_EventStrNtuple_To_FlatNtuple(){
  gBenchmark->Start("Timing");

  //  Generate  a dictionary, so that CINT will have all the information it needs about type or variable at anytime.
  gROOT->ProcessLine("#include <vector>"); 
  //  Same as above, since <float> goes by default
  //gInterpreter->GenerateDictionary("vector<float>");
  //  Generate dict for int's
  //gInterpreter->GenerateDictionary("vector<int>","vector");
  
  // General settings in case one adds drawing of histograms etc
  gStyle->SetOptFit();
  gStyle->SetOptStat(0);

  //
  // Open file with event structure 
  //

  TFile *myFile = nullptr;
  TFile *fileOut = nullptr;
  TTree *tree  = nullptr;
  TTree *electronTree_ = nullptr;  

  if (isDY){
    myFile  = new TFile(fileNameS);
    fileOut = TFile::Open(flatNtupleNameS,"RECREATE");
    printf("\n Work with true electrons\n");
  }
  else {
    myFile  = new TFile(fileNameBG);
    fileOut = TFile::Open(flatNtupleNameBG,"RECREATE");
    printf("\n Work with fake electrons\n");
  }
  tree     = (TTree*)myFile->Get(treeName);  
  if( !myFile || !tree){
    printf("Failed to open the input file or find the tree, check\n   %s or %s \n", 
	   myFile->GetName(), tree->GetName() );
    assert(0);
  }
  if (!fileOut) { return; }
  electronTree_  = new TTree("electronTree","Flat_ntuple");

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
  //std::vector <int> *isTrueElectronAlternative = 0;
  // Other vars  
  // Impact parameters
  std::vector <float> *eleD0 = 0;      // r-phi plane impact parameter
  std::vector <float> *eleDZ = 0;      // r-z plane impact parameter
  // Matching track-supercluster
  std::vector <float> *eleDEtaIn = 0;  // deltaEtaIn
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
  TBranch *b_eleDEtaIn = 0;
  TBranch *b_eleDPhiIn = 0;
  TBranch *b_eleHoverE = 0;
  TBranch *b_eleFull5x5SigmaIEtaIEta = 0;
  TBranch *b_eleOOEMOOP = 0;
  TBranch *b_eleExpectedMissingInnerHits = 0;
  TBranch *b_electronPassConversionVeto = 0;

  
  // Connect variables and branches to the tree with the data
  tree->SetBranchAddress("nEle", &eleNEle, &b_eleNEle);
  tree->SetBranchAddress("nPV", &nPV, &b_nPV);
  tree->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
  tree->SetBranchAddress("rho", &eleRho, &b_eleRho);
  tree->SetBranchAddress("pt", &elePt, &b_elePt);
  tree->SetBranchAddress("etaSC", &eleEtaSC, &b_eleEtaSC);
  tree->SetBranchAddress("phiSC", &elePhiSC, &b_elePhiSC);
  tree->SetBranchAddress("isoChargedHadrons", &isoChargedHadrons, &b_isoChargedHadrons);
  tree->SetBranchAddress("isoNeutralHadrons", &isoNeutralHadrons, &b_isoNeutralHadrons);
  tree->SetBranchAddress("isoPhotons",        &eleIsoPhotons,        &b_eleIsoPhotons);
  tree->SetBranchAddress("isTrue",    &eleIsTrueElectron,    &b_eleIsTrueElectron);
  tree->SetBranchAddress("d0",                &eleD0,             &b_eleD0);
  tree->SetBranchAddress("dz",                &eleDZ,             &b_eleDZ);
  tree->SetBranchAddress("dEtaIn",            &eleDEtaIn,         &b_eleDEtaIn);
  tree->SetBranchAddress("dPhiIn",            &eleDPhiIn,         &b_eleDPhiIn);
  tree->SetBranchAddress("hOverE",            &eleHoverE,         &b_eleHoverE);
  tree->SetBranchAddress("full5x5_sigmaIetaIeta", &eleFull5x5SigmaIEtaIEta,
			  &b_eleFull5x5SigmaIEtaIEta);
  tree->SetBranchAddress("ooEmooP",           &eleOOEMOOP,        &b_eleOOEMOOP);
  tree->SetBranchAddress("expectedMissingInnerHits", &eleExpectedMissingInnerHits, 
			  &b_eleExpectedMissingInnerHits);
  tree->SetBranchAddress("passConversionVeto",       &electronPassConversionVeto,
			  &b_electronPassConversionVeto);

  
  

  
  //
  //*************************************************************************//
  //
  // Here declare file, tree and variables for the new flat ntuple
  //
  //************************************************************************//
  //


  fileOut->cd();  //?????????
  
  // Vars for pile-up
  // Int_t nPUTrue_;    // true pile-up
  // Int_t nPU_;        // generated pile-up
  Int_t nPV_ =0;        // number of reconsrtucted primary vertices
  
  Float_t rho_ = 0.0;      // the rho variable
  Float_t gweight_ = 0.0;      // gen Weight

  Int_t nEle_ =0;
  // all electron variables
  Float_t pt_ = 0.0;
  Float_t etaSC_= 0.0;
  Float_t dEtaIn_= 0.0;
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
    
  // electronTree_->Branch("nPU"        ,  &nPU_     , "nPU/I");
  // electronTree_->Branch("nPUTrue"    ,  &nPUTrue_ , "nPUTrue/I");
  // electronTree_->Branch("sigmaIetaIeta",         &sigmaIetaIeta_, "sigmaIetaIeta/F");
  // electronTree_->Branch("isoChargedHadrons"      , &isoChargedHadrons_);
  // electronTree_->Branch("isoNeutralHadrons"      , &isoNeutralHadrons_);
  // electronTree_->Branch("isoPhotons"             , &isoPhotons_);
  // electronTree_->Branch("isoChargedFromPU"       , &isoChargedFromPU_);
  // electronTree_->Branch("relIsoWithDBeta"      , &relIsoWithDBeta_, "relIsoWithDBeta/F");

  electronTree_->Branch("genWeight"  ,  &gweight_ , "gweight/F");  
 
  electronTree_->Branch("nPV"        ,  &nPV_     , "nPV/I");
  electronTree_->Branch("rho"        ,  &rho_ , "rho/F");
  electronTree_->Branch("nEle"       ,  &nEle_ , "nEle/I");
  electronTree_->Branch("pt"    ,  &pt_    , "pt/F");			    
  electronTree_->Branch("etaSC" ,  &etaSC_ , "etaSC/F");
  electronTree_->Branch("dEtaIn",  &dEtaIn_, "dEtaIn/F");
  electronTree_->Branch("dPhiIn",  &dPhiIn_, "dPhiIn/F");
  electronTree_->Branch("hOverE",  &hOverE_, "hOverE/F");
  electronTree_->Branch("full5x5_sigmaIetaIeta", &full5x5_sigmaIetaIeta_, "full5x5_sigmaIetaIeta/F");
  electronTree_->Branch("relIsoWithEA"           , &relIsoWithEA_, "relIsoWithEA/F");
  electronTree_->Branch("ooEmooP", &ooEmooP_, "ooEmooP/F");
  electronTree_->Branch("d0"     , &d0_,      "d0/F");
  electronTree_->Branch("dz"     , &dz_,      "dz/F");
  electronTree_->Branch("expectedMissingInnerHits", &expectedMissingInnerHits_, "expectedMissingInnerHits/I");
  electronTree_->Branch("passConversionVeto", &passConversionVeto_, "passConversionVeto/I");
  electronTree_->Branch("isTrueEle"    , &isTrueEle_,     "isTrueEle/I");
  
  UInt_t maxEvents = tree->GetEntries();
  UInt_t maxEventsOver10000 =   maxEvents/10000.;

  if( smallEventCount )
    maxEvents = 5000000;

  printf("\nStart drawing into hists with total # of events = %u\n", maxEvents );

  //
  // Loop over events
  //

  for(UInt_t ievent = 0; ievent < maxEvents; ievent++){

    Long64_t tentry = tree->LoadTree(ievent);
    // Load the value of the number of the electrons in the event    
    b_eleNEle->GetEntry(tentry);
    
    if( ievent%maxEventsOver10000 == 0){
      printf("> "); fflush(stdout);
      //printf("Event %d, number of electrons %u\n", ievent, eleNEle);
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
    b_eleDEtaIn->GetEntry(tentry);
    b_eleDPhiIn->GetEntry(tentry);
    b_eleHoverE->GetEntry(tentry);
    b_eleFull5x5SigmaIEtaIEta->GetEntry(tentry);
    b_eleOOEMOOP->GetEntry(tentry);
    b_eleExpectedMissingInnerHits->GetEntry(tentry);
    b_electronPassConversionVeto->GetEntry(tentry);

    // Loop over the electrons
    for(int iele = 0; iele < eleNEle; iele++){
     
      gweight_ = genWeight;
      nEle_ = eleNEle;
      pt_ = elePt->at(iele);
      etaSC_ = eleEtaSC->at(iele);
      dEtaIn_ = eleDEtaIn->at(iele);
      dPhiIn_ = eleDPhiIn->at(iele);
      full5x5_sigmaIetaIeta_ = eleFull5x5SigmaIEtaIEta->at(iele);
      hOverE_ =  eleHoverE->at(iele);
      d0_ = eleD0->at(iele);
      dz_ = eleDZ->at(iele);
      isoChargedHadrons_ =  isoChargedHadrons->at(iele);
      isoNeutralHadrons_ = isoNeutralHadrons->at(iele);
      isoPhotons_ = eleIsoPhotons->at(iele);
      expectedMissingInnerHits_ = eleExpectedMissingInnerHits->at(iele);
      nPV_ = nPV;
      // // nPU_ = nPU->at(iele);
      // // nPUTrue_ = nPUTrue->at(iele);
      rho_ = eleRho; //->at(iele);
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
			       - rho_ * area ) )/pt_;
      
      electronTree_->Fill();// IK will kill me next time, if this line is not at the right place!
      
    } // end loop over the electrons
  
  } // end loop over SIGNAL events
  
  bazinga("I'm here to write signal tree");    
  electronTree_->Write();
  fileOut->Write();
  fileOut->Close();
  delete fileOut; 
  fileOut = nullptr;
  delete myFile;
  myFile = nullptr;
  
  bazinga("I'm finished with calculation, here is the info from dBenchmark:");  
  printf("\n");  gBenchmark->Show("Timing");  // get timing info
} // end of main fcn
