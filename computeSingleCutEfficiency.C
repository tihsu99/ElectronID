#include "TTree.h"
#include "TString.h"
#include "TFile.h"

#include "OptimizationConstants.hh"
#include "VarCut.hh"

// Files with signal and background trees (ideally the ntuples
// that were used for TMVA optimization
const TString fnameSignal = "../../ntuples/DYJetsToLL_50ns.root";
const TString signalTreeName = "ntupler/ElectronTree";
const TString fnameBackground = "../../ntuples/TTJets_50ns.root";
const TString backgroundTreeName = "ntupler/ElectronTree";

// Forward declarations
TTree *getTreeFromFile(TString fname, TString tname);
void   findEfficiencies(TTree *signalTree, TTree *backgroundTree, bool doBarrel,
			float &effSignal, float &effBackground, 
			VarCut *cutObject);
//
// Main function
//
void computeSingleCutEfficiency(TString var, float cutVal, bool doBarrel){

  // Get trees
  printf("\n Take true electrons from %s   tree %s\n", 
	 fnameSignal.Data(), signalTreeName.Data());
  TTree *signalTree = getTreeFromFile( fnameSignal, signalTreeName);
  // Input background tree  
  printf("\n Take background electrons from %s   tree %s\n", 
	 fnameBackground.Data(), backgroundTreeName.Data());
  TTree *backgroundTree = getTreeFromFile( fnameBackground, backgroundTreeName);

  VarCut *cutObject = new VarCut();
  for(int ivar= 0; ivar<Vars::nVariables; ivar++){
    if( var == Vars::variables[ivar]->name )
      cutObject->setCutValue( Vars::variables[ivar]->name, cutVal);
    else
    cutObject->setCutValue( Vars::variables[ivar]->name, 1e30);
  }
  printf("Compute the efficiency of the following cut set:\n");
  cutObject->print();
  
  // Compute the efficiencies
  float effSignal, effBackground;
  findEfficiencies(signalTree, backgroundTree, doBarrel,
		   effSignal, effBackground,
		   cutObject);
  printf("Computed eff: effS= %.3f effB= %.3f\n",
	 effSignal, effBackground);
  

}

// Get a given tree from a given file name.
TTree *getTreeFromFile(TString fname, TString tname){

  TFile *file = new TFile( fname );
  TTree *tree     = (TTree*) file->Get(tname);
  
  return tree;
}

// Compute signal and background efficiencies for given cuts
void findEfficiencies(TTree *signalTree, TTree *backgroundTree, bool doBarrel,
		      float &effSignal, float &effBackground, VarCut *cutObject){

  TCut etaCut = "";
  if( doBarrel ){
    etaCut = Opt::etaCutBarrel;
  }else{
    etaCut = Opt::etaCutEndcap;
  }
  TCut kinematicCuts = Opt::ptCut && etaCut;

  TCut preselectionCuts = kinematicCuts && Opt::otherPreselectionCuts;
  
  TCut signalCuts = preselectionCuts && Opt::trueEleCut;
  TCut backgroundCuts = preselectionCuts && Opt::fakeEleCut;  
 
  TCut selectionCuts = *(cutObject->getCut());

  printf("\nSelecton cuts: %s\n",
	 selectionCuts.GetTitle());
  printf("\nSignal cuts: %s\n",
	 signalCuts.GetTitle());
  printf("\nBackground cuts: %s\n",
	 backgroundCuts.GetTitle());

  effSignal = (1.0*signalTree->GetEntries(selectionCuts && signalCuts) )
    / signalTree->GetEntries(signalCuts);

  effBackground = (1.0*backgroundTree->GetEntries(selectionCuts 
						  && backgroundCuts) )
    / backgroundTree->GetEntries(backgroundCuts);
  
  return;
}
