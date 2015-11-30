//
// This code finds individual cut values that correspond to a given
// cut efficiency (such as 99.9%)
//
#include "TString.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TCut.h"

#include "Variables.hh"
#include "OptimizationConstants.hh"
#include "VarCut.hh"

// Define unique part of the file name for saving the cuts
TString dateTag = "20161101_090000";

// For debug purposes, set the flag below to true, for regular
// computation set it to false
const bool useSmallEventCount = false;

// Forward declarations
void findVarLimits(TString var, bool isBarrel, float &xmin, float &xmax);

//
// A helper class VarInfo
// 
class VarInfo {

public:
  VarInfo(TString varname, double xlow, double xhigh);
  ~VarInfo();

  double findUpperCut(TTree *tree, TCut &preselectionCut, double eff);
  
private:
  
  // Data members
  TString _varname;
  double _xmin;
  double _xmax;
  int    _nbins;
  
  TH1F *_hist;

};

VarInfo::VarInfo(TString varname, double xmin, double xmax) :
  _varname (varname), _xmin (xmin), _xmax (xmax) 
{
  // Number of bins should be large enough
  _nbins = 100000;

  TString hname = "hist_";
  // Get rid of () in variable name to use it in hist name
  TString flatVarname = _varname;
  flatVarname.ReplaceAll("(","_");
  flatVarname.ReplaceAll(")","_");
  hname += flatVarname;
  _hist = new TH1F(hname, "", _nbins, _xmin, _xmax);
}

VarInfo::~VarInfo(){
  delete _hist;
}

double VarInfo::findUpperCut(TTree *tree, TCut &preselectionCut, double eff){

  if( tree == 0 )
    assert(0);

  // Make the histogram
  TString drawCommand = _varname;
  drawCommand += ">>";
  drawCommand += _hist->GetName();
  if( !useSmallEventCount ){
    tree->Draw(drawCommand, preselectionCut);
  }else{
    printf("DEBUG MODE: using small event count\n");
    tree->Draw(drawCommand, preselectionCut, "", 100000);
  }
  
  // DEBUG
  // preselectionCut.Print();

  // Start from the highest bin, and walk downward  
  double sumTotal = _hist->GetSumOfWeights() 
    + _hist->GetBinContent(0) // underflow
    + _hist->GetBinContent(_nbins+1); //overflow
  double sumRunning = _hist->GetBinContent(_nbins+1);
  double effObserved;;
  double cut = _xmax;
  // From the last bin to the second bin (avoid first, in case almost
  // all entries are sitting at var=0, like for H/E
  for(int i=_nbins; i>=2; i--){
    sumRunning += _hist->GetBinContent(i);
    effObserved  = sumRunning/sumTotal;
    if( effObserved > (1-eff) )
      break;
    cut = _hist->GetXaxis()->GetBinLowEdge(i);
  };
  
  // Check
  TCut newCut(TString::Format("%s < %f", _varname.Data(), cut ));
  double num = tree->GetEntries( preselectionCut && newCut );
  double denom = tree->GetEntries( preselectionCut );
  double effCheck = num/denom;
  // printf("Preselectioncut is\n");
  // preselectionCut.Print();
  printf("Found the cut for variable %30s: requested eff= %.4f, observed= %.4f (%.0f/%.0f), cut= %f\n",
	 _varname.Data(), eff, effCheck, num, denom, cut);

  return cut;
}

//
// Main function
//
void findCutLimits(){

  // Get the signal tree
  TFile *input = new TFile( Opt::fnameSignal );
  TTree *tree = (TTree*)input->Get(Opt::signalTreeName);
  if( !tree ) 
    assert(0);

  const float eff = 0.999;

  //
  // Barrel first
  //
  printf("\nFind cuts for barrel electrons\n");
  TCut preselectionCuts = Opt::ptCut && Opt::etaCutBarrel 
    && Opt::otherPreselectionCuts && Opt::trueEleCut;
  bool useBarrel = true;

  VarCut *cutBarrel0999 = new VarCut();
  for(int i=0; i<Vars::nVariables; i++){
    float xlow = 0;
    float xhigh = 1000; // just a large number, overwritten below
    findVarLimits( Vars::variables[i]->name, useBarrel, xlow, xhigh);
    // Note: use nameTmva below, so that the var string will contain
    // the abs() as needed.
    VarInfo var(Vars::variables[i]->nameTmva, xlow, xhigh);
    float cutValue = var.findUpperCut(tree, preselectionCuts, eff);
    cutBarrel0999->setCutValue(Vars::variables[i]->name, cutValue);
  }
  cutBarrel0999->print();

  // //
  // // Endcap second
  // //
  // printf("\nFind cuts for endcap electrons\n");
  // preselectionCuts = Opt::ptCut && Opt::etaCutEndcap
  //   && Opt::otherPreselectionCuts && Opt::trueEleCut;
  // useBarrel = false;

  // VarCut *cutEndcap0999 = new VarCut();
  // for(int i=0; i<Vars::nVariables-8; i++){
  //   float xlow = 0;
  //   float xhigh = 1000; // just a large number, overwritten below
  //   findVarLimits( Vars::variables[i]->name, useBarrel, xlow, xhigh);
  //   // Note: use nameTmva below, so that the var string will contain
  //   // the abs() as needed.
  //   VarInfo var(Vars::variables[i]->nameTmva, xlow, xhigh);
  //   float cutValue = var.findUpperCut(tree, preselectionCuts, eff);
  //   cutEndcap0999->setCutValue(Vars::variables[i]->name, cutValue);
  // }
  // cutEndcap0999->print();

  //
  // Save the cuts
  //
  TString fileBarrelName = Opt::cutRepositoryDir + TString("/")
    + TString("cuts_barrel_eff_0999_") + dateTag + TString(".root");
  printf("\nSave barrel cuts into the file %s\n", fileBarrelName.Data());
  TFile *fileBarrel = new TFile(fileBarrelName, "recreate");
  if(fileBarrel == 0)
    assert(0);
  cutBarrel0999->Write("cuts");
  fileBarrel->Close();

  // TString fileEndcapName = Opt::cutRepositoryDir + TString("/")
  //   + TString("cuts_endcap_eff_0999_") + dateTag + TString(".root");
  // printf("\nSave endcap cuts into the file %s\n", fileEndcapName.Data());
  // TFile *fileEndcap = new TFile(fileEndcapName, "recreate");
  // if(fileEndcap == 0)
  //   assert(0);
  // cutEndcap0999->Write("cuts");
  // fileEndcap->Close();

}

void findVarLimits(TString var, bool useBarrel, float &xlow, float &xhigh){

  xlow = 0;
  if ( var == "full5x5_sigmaIetaIeta" ){
    xhigh = useBarrel ? 0.03 : 0.1;
  } else if ( var == "dEtaIn"         ){      
    xhigh = useBarrel ? 0.05 : 0.1;
  } else if ( var == "dPhiIn"         ){      
    xhigh = 0.4;
  } else if ( var == "hOverE"         ){      
    xhigh = useBarrel ? 0.5 : 1.0;
  } else if ( var == "relIsoWithEA"){      
    xhigh = 2;
  } else if ( var == "ooEmooP"        ){      
    xhigh = 0.5;
  } else if ( var == "d0"             ){      
    xhigh = useBarrel ? 0.2 : 0.4;
  } else if ( var == "dz"             ){      
    xhigh = 5;
  } else if ( var == "expectedMissingInnerHits"){
    xhigh = 5;
  } else {
    printf("ERROR: can not set var limits properly, unknown variable.\n");
    xhigh = 1000;
  };

  return;
}

