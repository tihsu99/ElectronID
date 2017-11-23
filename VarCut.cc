#include "VarCut.hh"
#include <iomanip>
#include <iostream>
#include <sstream>

const int sig_figs = 3;
  
auto numberNew =   [](int number_of_sig_figs, float number){
  std::stringstream myStream;
  myStream << std::setprecision(number_of_sig_figs) << number;
  return std::stod(myStream.str());
};

const int UNDEFCUT = -999;

// All cuts are initialized to an unlikely value
VarCut::VarCut(){
  for(int i=0; i<Vars::nVariables; i++) _cuts[i]= UNDEFCUT;
};

// Construct TCut object for all cuts joined with &&
TCut *VarCut::getCut(TString selectVar){

  TCut *cut = 0;

  // Die if something appears uninitialized
  for(int i=0; i<Vars::nVariables; i++){
    if(_cuts[i] == UNDEFCUT){
      printf("VarCut:: not all cuts are set! Die!\n");
      assert(0);
    }
  }

  cut = new TCut("");
  for(int i=0; i<Vars::nVariables; i++){
    if(selectVar != "" and Vars::variables[i]->name != selectVar) continue;
    // The += adds all cuts with &&:
    (*cut) += TString::Format("%s<%f", Vars::variables[i]->nameTmva.Data(), _cuts[i]);
  }
  
  return cut;
}

void VarCut::setCutValue(TString varName, float val){
  _cuts[getVariableIndex(varName)] = numberNew(sig_figs, val);
}

float VarCut::getCutValue(TString variable){
  return numberNew(sig_figs,  _cuts[getVariableIndex(variable)]);
}

int VarCut::getVariableIndex(TString variable){
  for(int i=0; i<Vars::nVariables; i++){
    if(variable == Vars::variables[i]->name or variable == Vars::variables[i]->nameTmva) return i;
  }
  printf("VarCut::getVariableIndex: requested variable is not known!!!\n");
  exit(1);
}


bool VarCut::isSymmetric(TString variable){
  return Vars::variables[getVariableIndex(variable)]->symmetric;
}

// Print all cut values to stdout
void VarCut::printCuts(){
  printf("VarCut::print: Cut values are\n");
  for(int i=0; i<Vars::nVariables; i++){
    printf("  %30s < %g\n", Vars::variables[i]->nameTmva.Data(), numberNew(sig_figs, _cuts[i]));
  }
}
