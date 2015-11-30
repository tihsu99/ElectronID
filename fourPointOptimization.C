#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "OptimizationConstants.hh"
#include "VariableLimits.hh"
#include "optimize.hh"

const int optimizeTightOnly = 1;

void fourPointOptimization(){

  // Define source for the initial cut range
  TString startingCutMaxFileName 
    = "cuts_barrel_eff_0999_20161101_090000.root";
  if( !Opt::useBarrel )
    startingCutMaxFileName 
      = "cuts_endcap_eff_0999_20161101_090000.root";

  TString namePrefix = "cuts_barrel_";
  if( !Opt::useBarrel )
    namePrefix = "cuts_endcap_";
  TString namePass[Opt::nWP] = {"pass1_","pass2_","pass3_","pass4_"};
  TString nameTime = "20151101_100000";

  for( int ipass = 0; ipass < optimizeTightOnly; ipass++){

    // This string is the file name that contains the ROOT file
    // with the VarCut object that defines the range of cut variation.
    // Note: for each subsequence pass, the previous working point
    // is used. For the first pass, the 99.9% efficient cut range is used.
    TString cutMaxFileName = startingCutMaxFileName;
    if( ipass > 0 ){
      cutMaxFileName = namePrefix + namePass[ipass-1] + nameTime + TString("_")
	+ Opt::wpNames[ipass-1] + TString(".root");
    }
    
    // The string below is used to construct the file names
    // to save the cut objects
    TString cutOutputBase = namePrefix;
    cutOutputBase += namePass[ipass];
    cutOutputBase += nameTime;

    // This string will be used to construct the dir for the output
    // of TMVA: the dir for weights and the filename for diagnostics
    TString trainingDataOutputBase = "training_results_";
    if( Opt::useBarrel ){
      trainingDataOutputBase += "barrel_";
    }else{
      trainingDataOutputBase += "endcap_";
    }
    trainingDataOutputBase += namePass[ipass];
    trainingDataOutputBase += nameTime;

    // Use the following user-defined limits
    // (if in doubt, use the no-restrictions one defined in VariableLimits.hh)
    VarLims::VariableLimits **userDefinedCutLimits = VarLims::limitsNoRestrictions;
    if( ipass > 0 ){
      userDefinedCutLimits = VarLims::limitsWPAnyV1;
    }

    if( ipass == 3 ){
      userDefinedCutLimits = VarLims::limitsWP4endcap;	
      if (Opt::useBarrel)      userDefinedCutLimits = VarLims::limitsWP4barrel;	
    }

    printf("\n-----------------------------------------------------------------\n");
    printf("\n");
    printf("    Run new optimization pass  \n");
    printf("\n");
    printf(" Input file that defines cut optimization limits:");
    printf(" %s\n", cutMaxFileName.Data());
    printf(" Base for the file name of output cuts          :");
    printf(" %s\n", cutOutputBase.Data());
    printf("------------------------------------------------------------------\n\n");
    
    optimize(cutMaxFileName, cutOutputBase, trainingDataOutputBase,
    	     userDefinedCutLimits);    

  }
 
  // Finally, create the files containing the working point
  // by copying the appropriate pass files.
  // The first working point is output of pass1, the second of pass2, etc.
  // All other working poitns of all passes are ignored
  printf("\n");
  printf("====================================================\n");
  printf("Final definition of working points\n");
  printf("====================================================\n");
  printf("Copy files with working points info into the final locations\n");
  for(int i=0; i<optimizeTightOnly /*Opt::nWP*/; i++){

    TString wpPassFileName = Opt::cutRepositoryDir + TString("/")
      + namePrefix + namePass[i] + nameTime + TString("_")
      + Opt::wpNames[i] + TString(".root");

    TString wpFinalFileName =  Opt::cutRepositoryDir + TString("/")
      + namePrefix + nameTime + TString("_")
      + Opt::wpNames[i] + TString(".root");

    TString copyCommand = TString::Format("cp %s %s", 
					  wpPassFileName.Data(), 
					  wpFinalFileName.Data());

    int result = gSystem->Exec(copyCommand); 
    int tries =0;
    while( result != 0 && tries <10) {
      printf("Copy command failed!!! The command: %s\n", copyCommand.Data());
      result = gSystem->Exec(copyCommand); 
      tries += 1;
    }

    printf("\nFinal definition for working point %s\n", Opt::wpNames[i].Data());
    printf(" file name:   %s\n", wpFinalFileName.Data());
    TFile *fwp = new TFile(wpFinalFileName);
    VarCut *thisWP = (VarCut*)fwp->Get("cuts");
    if( thisWP != 0 )
      thisWP->print();
    else
      printf("???? not found????\n");
    fwp->Close();
    
  }
 
}

