#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "OptimizationConstants.hh"
#include "VariableLimits.hh"
#include "optimize.hh"

void simpleOptimization(){
  
  // Define source for the initial cut range
  TString startingCutMaxFileName 
    = "cuts_barrel_eff_0999_20140729_111500.root";
  if( !Opt::useBarrel )
    startingCutMaxFileName 
      = "cuts_endcap_eff_0999_20140729_111500.root";

  TString namePrefix = "cuts_barrel_";
  if( !Opt::useBarrel )
    namePrefix = "cuts_endcap_";
  TString nameTime = "20140730_200000";
  
  // This string is the file name that contains the ROOT file
  // with the VarCut object that defines the range of cut variation.
  TString cutMaxFileName = startingCutMaxFileName;
  
  // The string below is used to construct the file names
  // to save the cut objects
  TString cutOutputBase = namePrefix + nameTime;
    
  // This string will be used to construct the dir for the output
  // of TMVA: the dir for weights and the filename for diagnostics
  TString trainingDataOutputBase = "training_results_";
  if( Opt::useBarrel ){
    trainingDataOutputBase += "barrel_";
  }else{
    trainingDataOutputBase += "endcap_";
  }
  trainingDataOutputBase += nameTime;

  // Use the following user-defined limits
  // (if in doubt, use the no-restrictions one defined in VariableLimits.hh)
  VarLims::VariableLimits **userDefinedCutLimits = VarLims::limitsNoRestrictions;

  printf("\n-----------------------------------------------------------------\n");
  printf("\n");
  printf("    Run optimization, single pass  \n");
  printf("\n");
  printf(" Input file that defines cut optimization limits:");
  printf(" %s\n", cutMaxFileName.Data());
  printf(" Base for the file name of output cuts          :");
  printf(" %s\n", cutOutputBase.Data());
  printf("------------------------------------------------------------------\n\n");
    
  optimize(cutMaxFileName, cutOutputBase, trainingDataOutputBase,
	   userDefinedCutLimits);    
    
}

