#include "Variables.hh"
#include "OptimizationConstants.hh"
#include "VarCut.hh"
#include "TFile.h"
#include "TString.h"


const TString cutFileNamesBarrel[Opt::nWP] = {
  "cut_repository/cuts_barrel_20160616_200000_WP_Veto",
  "cut_repository/cuts_barrel_20160616_200000_WP_Loose",
  "cut_repository/cuts_barrel_20160616_200000_WP_Medium",
  "cut_repository/cuts_barrel_20160616_200000_WP_Tight"
};
  
const TString cutFileNamesEndcap[Opt::nWP] = {
  "cut_repository/cuts_endcap_20160616_200000_WP_Veto",
  "cut_repository/cuts_endcap_20160616_200000_WP_Loose",
  "cut_repository/cuts_endcap_20160616_200000_WP_Medium",
  "cut_repository/cuts_endcap_20160616_200000_WP_Tight"
};

// Missing hits is unused
const int missingHitsBarrel[Opt::nWP] = { 2, 1, 1, 1};
const int missingHitsEndcap[Opt::nWP] = { 3, 1, 1, 1};

const float d0Barrel[Opt::nWP] = {0.050, 0.0275, 0.0168, 0.00787};
const float d0Endcap[Opt::nWP] = {0.100, 0.100, 0.0653, 0.0273};

const float dzBarrel[Opt::nWP] = {0.100, 0.100, 0.100, 0.100};
const float dzEndcap[Opt::nWP] = {0.200, 0.200, 0.200, 0.200};

void repackageCuts(){

  // Repack barrel
  for(int iWP = 0; iWP<Opt::nWP; iWP++){

    TString cutFileNameBase = cutFileNamesBarrel[iWP];
    TString cutFileIn  = cutFileNameBase + TString(".root");
    TString cutFileOut = cutFileNameBase + TString("_adjusted.root");
    
    TFile fin(cutFileIn);
    VarCut *originalCut = (VarCut*)fin.Get("cuts");

    TFile fout(cutFileOut, "recreate");
    VarCut *repackedCut = new VarCut();
    for(int ivar = 0; ivar < Vars::nVariables; ivar++){
      TString vname = Vars::variables[ivar]->name;
      if( vname == "d0" ){
	repackedCut->setCutValue("d0", d0Barrel[iWP]);
      }else if (vname == "dz"){
	repackedCut->setCutValue("dz", dzBarrel[iWP]);
      }else{
	repackedCut->setCutValue(vname, originalCut->getCutValue(vname));
      }
    } // end loop over variables
    repackedCut->Write("cuts");
    fout.Close();
    fin.Close();
  } // end loop over working points


  // Repack endcap
  for(int iWP = 0; iWP<Opt::nWP; iWP++){

    TString cutFileNameBase = cutFileNamesEndcap[iWP];
    TString cutFileIn  = cutFileNameBase + TString(".root");
    TString cutFileOut = cutFileNameBase + TString("_adjusted.root");
    
    TFile fin(cutFileIn);
    VarCut *originalCut = (VarCut*)fin.Get("cuts");

    TFile fout(cutFileOut, "recreate");
    VarCut *repackedCut = new VarCut();
    for(int ivar = 0; ivar < Vars::nVariables; ivar++){
      TString vname = Vars::variables[ivar]->name;
      if( vname == "d0" ){
	printf("NOTE: filling d0 values with hardwired numbers! Make sure they are correct!\n");
	repackedCut->setCutValue("d0", d0Endcap[iWP]);
      }else if (vname == "dz"){
	printf("NOTE: filling dz values with hardwired numbers! Make sure they are correct!\n");
	repackedCut->setCutValue("dz", dzEndcap[iWP]);
      }else{
	repackedCut->setCutValue(vname, originalCut->getCutValue(vname));
      }
    } // end loop over variables
    repackedCut->Write("cuts");
    fout.Close();
    fin.Close();
  } // end loop over working points



}

