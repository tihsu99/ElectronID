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
TString dateTag = "2017-11-16";

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
    double  _xmin,_xmax;
    int     _nbins;
    TH1F*   _hist;
};

VarInfo::VarInfo(TString varname, double xmin, double xmax) :
  _varname (varname), _xmin (xmin), _xmax (xmax){

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

  tree->Draw(_varname + ">>" + _hist->GetName(), preselectionCut, "goff");

  // Start from the highest bin, and walk downward
  double sumTotal   = _hist->GetSumOfWeights() + _hist->GetBinContent(0) + _hist->GetBinContent(_nbins+1); // underflow and overflow
  double sumRunning = _hist->GetBinContent(_nbins+1);
  double cut = _xmax;

  // From the last bin to the second bin (avoid first, in case almost
  // all entries are sitting at var=0, like for H/E
  for(int i=_nbins; i>=2; i--){
    sumRunning += _hist->GetBinContent(i);
    if(sumRunning/sumTotal> (1-eff)) break;
    cut = _hist->GetXaxis()->GetBinLowEdge(i);
  };

  // Check
  TCut newCut(TString::Format("%s < %f", _varname.Data(), cut));
  double num      = tree->GetEntries(preselectionCut && newCut);
  double denom    = tree->GetEntries(preselectionCut);
  double effCheck = num/denom;
  printf("Found the cut for variable %30s: requested eff= %.4f, observed= %.4f (%.0f/%.0f), cut= %f\n", _varname.Data(), eff, effCheck, num, denom, cut);

  return cut;
}

void writeCutAtEff(float eff, bool useBarrel, TTree* tree, TCut preselectionCuts, TString name){
  VarCut *cutAtEff = new VarCut();
  for(int i=0; i<Vars::nVariables; i++){
    float xlow = 0;
    float xhigh = 1000; // just a large number, overwritten below
    findVarLimits(Vars::variables[i]->name, useBarrel, xlow, xhigh);
    // Note: use nameTmva below, so that the var string will contain
    // the abs() as needed.
    VarInfo var(Vars::variables[i]->nameTmva, xlow, xhigh);
    float cutValue = var.findUpperCut(tree, preselectionCuts, eff);
    cutAtEff->setCutValue(Vars::variables[i]->name, cutValue);
  }
  cutAtEff->printCuts();

  TString fileName = Opt::cutRepositoryDir + TString("/")  + name;
  TFile *file = new TFile(fileName, "recreate");
  if(file == 0) assert(0);
  cutAtEff->Write("cuts");
  file->Close();
}

//
// Main function
//
void findCutLimits(){
  // Get the signal trees
  TFile *inputBarrel = new TFile( Opt::fnameSignalBarrel );
  TTree *treeBarrel = (TTree*)inputBarrel->Get(Opt::signalTreeName);
  if( !treeBarrel ) assert(0);
  TFile *inputEndcap = new TFile( Opt::fnameSignalEndcap );
  TTree *treeEndcap = (TTree*)inputEndcap->Get(Opt::signalTreeName);
  if( !treeEndcap ) assert(0);

  //
  // Barrel first
  //
  writeCutAtEff(0.999, true,  treeBarrel, Opt::ptCut && Opt::etaCutBarrel && Opt::otherPreselectionCuts && Opt::trueEleCut, TString("cuts_barrel_eff_0999_") + dateTag + ".root");
  writeCutAtEff(0.999, false, treeEndcap, Opt::ptCut && Opt::etaCutEndcap && Opt::otherPreselectionCuts && Opt::trueEleCut, TString("cuts_endcap_eff_0999_") + dateTag + ".root");
}

void findVarLimits(TString var, bool useBarrel, float &xlow, float &xhigh){
  xlow = 0;
  if ( var == "full5x5_sigmaIetaIeta" )        xhigh = useBarrel ? 0.03 : 0.1;
  else if ( var == "dEtaSeed"         )        xhigh = useBarrel ? 0.05 : 0.1;
  else if ( var == "dPhiIn"         )          xhigh = 0.4;
  else if ( var == "hOverE"         )          xhigh = useBarrel ? 0.5 : 1.0;
  else if ( var == "relIsoWithEA")             xhigh = 2;
  else if ( var == "ooEmooP"        )          xhigh = 0.5;
  else if ( var == "d0"             )          xhigh = useBarrel ? 0.2 : 0.4;
  else if ( var == "dz"             )          xhigh = 5;
  else if ( var == "expectedMissingInnerHits") xhigh = 5;
  else {
    printf("ERROR: can not set var limits properly, unknown variable.\n");
    xhigh = 1000;
  };

  return;
}

