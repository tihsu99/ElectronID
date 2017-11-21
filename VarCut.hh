#ifndef VARCUT_HH
#define VARCUT_HH

#include "TCut.h"
#include <TObject.h>

#include "Variables.hh"

class VarCut :public TObject {

public:
  VarCut();

  // Set value of the variable in the internal array based on its name,
  // using the regular name or the name known to TMVA (may include abs())
  void  setCutValue(TString varName, float val);

  // Look up cut value for given variable (regular name)
  float getCutValue(TString var);

  // Get the full TCut object with cuts on all variables
  // Get the TCut for a given variable if selectVar is specified
  TCut* getCut(TString selectVar = "");

  // Get index of the variable in the internal array from its name,
  // using the regular name or the name known to TMVA (may include abs())
  int getVariableIndex(TString var);

  // Does the cut involve abs()? i.e., is the cut symmetric wrt 0?
  bool isSymmetric(TString variable);

  // Input/output
  void print(); // print to stdout

private:
  // The actual list of variables for which cuts are stored here
  // is found in Variables.hh
  float _cuts[Vars::nVariables];

  ClassDef(VarCut,1)
};

#endif
