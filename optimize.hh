#ifndef OPTIMIZE_HH
#define OPTIMIZE_HH

#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include "Variables.hh"
#include "VariableLimits.hh"
#include "VarCut.hh"
#include "OptimizationConstants.hh"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Config.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"
#endif

// Forward declarations
// Input
TTree * getTreeFromFile(TString fname, TString tname, TFile **fileHandle);
// Configuration for TMVA
void    configureVariables(TMVA::DataLoader *factory);
TString getTrainAndTestOptions(int region);
void    configureCuts(TCut &signalCuts, TCut &backgroundCuts, int region);
TString getMethodOptions(TString cutMaxFileName, VarLims::VariableLimits **userDefinedCutLimits);

// Output
void writeWorkingPoints(const TMVA::Factory *factory, TString cutsOutFileNameBase, int region);

// Main method
void optimize(TString cutMaxFileName = "cuts_barrel_eff_0999_20140727_165000.root",
	      TString cutsOutFileNameBase = "cuts_barrel_xxxx_99999999_999999",
	      TString trainingDataOutputBase = "training_xxxx_results_99999999_999999",
	      VarLims::VariableLimits **userDefinedCutLimits = VarLims::limitsNoRestrictions, int region=0);

#endif
