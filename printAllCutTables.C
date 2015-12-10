#include "TString.h"
#include "TFile.h"
#include "VarCut.hh"
#include "Variables.hh"

const int nWP = 4;
const TString workingPointNames[nWP] = {
  "veto",
  "loose",
  "medium",
  "tight"
};


const TString barrelCutFiles[nWP] = {
  // only WP Veto was optimized, other files are only for aestetics purpose here
  "root://cmsxrootd.fnal.gov///store/user/cmsdas/2016/SHORT_EXERCISES/ElectronsAndPhotons/cut_repository/cuts_barrel_20151101_100000_WP_Veto.root",
  "root://cmsxrootd.fnal.gov///store/user/cmsdas/2016/SHORT_EXERCISES/ElectronsAndPhotons/cut_repository/cuts_barrel_20151101_100000_WP_Veto.root",
  "root://cmsxrootd.fnal.gov///store/user/cmsdas/2016/SHORT_EXERCISES/ElectronsAndPhotons/cut_repository/cuts_barrel_20151101_100000_WP_Veto.root",
  "root://cmsxrootd.fnal.gov///store/user/cmsdas/2016/SHORT_EXERCISES/ElectronsAndPhotons/cut_repository/cuts_barrel_20151101_100000_WP_Veto.root"

};

const TString endcapCutFiles[nWP] = {
  // change to endcap files once the optimzation is made
  "root://cmsxrootd.fnal.gov///store/user/cmsdas/2016/SHORT_EXERCISES/ElectronsAndPhotons/cut_repository/cuts_barrel_20151101_100000_WP_Veto.root",
  "root://cmsxrootd.fnal.gov///store/user/cmsdas/2016/SHORT_EXERCISES/ElectronsAndPhotons/cut_repository/cuts_barrel_20151101_100000_WP_Veto.root",
  "root://cmsxrootd.fnal.gov///store/user/cmsdas/2016/SHORT_EXERCISES/ElectronsAndPhotons/cut_repository/cuts_barrel_20151101_100000_WP_Veto.root",
  "root://cmsxrootd.fnal.gov///store/user/cmsdas/2016/SHORT_EXERCISES/ElectronsAndPhotons/cut_repository/cuts_barrel_20151101_100000_WP_Veto.root"
};


enum WorkingPointNames {
  WP_VETO = 0,
  WP_LOOSE,
  WP_MEDIUM,
  WP_TIGHT
};

void printAllCutTables(){

  TFile *barrelFiles[nWP];
  TFile *endcapFiles[nWP];
  VarCut *barrelCuts[nWP];
  VarCut *endcapCuts[nWP];
  for(int i=0; i<nWP; i++){

    // Load all barrel cuts
    TString barrelFileName =  barrelCutFiles[i];
    barrelFiles[i] = new TFile(barrelFileName);
    if( !barrelFiles[i] ){
      printf("Could not open the file %s\n", barrelFileName.Data());
      assert(0);
    }
    barrelCuts[i] = (VarCut*) barrelFiles[i]->Get("cuts");
    if( !barrelCuts[i] ){
      printf("Could not find the cuts object in file %s\n", barrelFileName.Data());
      assert(0);
    }

    // Load all endcap cuts
    TString endcapFileName =  endcapCutFiles[i];
    endcapFiles[i] = new TFile(endcapFileName);
    if( !endcapFiles[i] ){
      printf("Could not open the file %s\n", endcapFileName.Data());
      assert(0);
    }
    endcapCuts[i] = (VarCut*) endcapFiles[i]->Get("cuts");
    if( !endcapCuts[i] ){
      printf("Could not find the cuts object in file %s\n", endcapFileName.Data());
      assert(0);
    }

  }

  // Print barrel table for the twiki
  printf("\nAll cuts for barrel in a twiki table format:\n\n");
  printf("|                                    |   Veto     |    Loose   |   Medium   |  Tight   |\n");
  for(int i=0; i< Vars::nVariables; i++){
    TString variableName     = Vars::variables[i]->name.Data();
    TString variableTmvaName = Vars::variables[i]->nameTmva.Data();
    TString sign = (i == Vars::nVariables-1)?"<=":"<";
    printf("|  %30s %s  |  %g  |  %g  |  %g  | %g |\n", 
	   variableTmvaName.Data(), 
	   sign.Data(),
	   barrelCuts[WP_VETO  ]->getCutValue(variableName),
	   barrelCuts[WP_LOOSE ]->getCutValue(variableName),
	   barrelCuts[WP_MEDIUM]->getCutValue(variableName),
	   barrelCuts[WP_TIGHT ]->getCutValue(variableName)
	   );
  }
  printf("\n");
  
  // Print endcap table for the twiki
  printf("\nAll cuts for endcap in a twiki table format:\n\n");
  printf("|                                    |   Veto     |    Loose   |   Medium   |  Tight   |\n");
  for(int i=0; i< Vars::nVariables; i++){
    TString variableName     = Vars::variables[i]->name.Data();
    TString variableTmvaName = Vars::variables[i]->nameTmva.Data();
    TString sign = (i == Vars::nVariables-1)?"<=":"<";
    printf("|  %30s %s  |  %g  |  %g  |  %g  | %g |\n", 
	   variableTmvaName.Data(), 
	   sign.Data(),
	   endcapCuts[WP_VETO  ]->getCutValue(variableName),
	   endcapCuts[WP_LOOSE ]->getCutValue(variableName),
	   endcapCuts[WP_MEDIUM]->getCutValue(variableName),
	   endcapCuts[WP_TIGHT ]->getCutValue(variableName)
	   );
  }
  printf("\n");
  
}
