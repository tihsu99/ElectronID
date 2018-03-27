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
  for(int i=0; i<Vars::nVariables; i++) _cuts[i]      = UNDEFCUT;
  for(int i=0; i<Vars::nConstants; i++) _constants[i] = UNDEFCUT;
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
    if(Vars::variables[i]->nameTmva == "hOverE" and _constants[getConstantIndex("C_rho")] > 0){
      (*cut) += TString::Format("%s<(%f+%f/eSC+%f*rho/eSC)", Vars::variables[i]->nameTmva.Data(), _cuts[i], _constants[getConstantIndex("C_E")], _constants[getConstantIndex("C_rho")]);
    } else if(Vars::variables[i]->nameTmva == "relIsoWithEA" and _constants[getConstantIndex("C_pt")] > 0){
      (*cut) += TString::Format("%s<(%f+%f/pt)", Vars::variables[i]->nameTmva.Data(), _cuts[i], _constants[getConstantIndex("C_pt")]);
    } else {
      (*cut) += TString::Format("%s<%f", Vars::variables[i]->nameTmva.Data(), _cuts[i]);
    }
  }
  
  return cut;
}

void VarCut::setCutValue(TString varName, float val){
  _cuts[getVariableIndex(varName)] = numberNew(sig_figs, val);
}

float VarCut::getCutValue(TString variable){
  return numberNew(sig_figs,  _cuts[getVariableIndex(variable)]);
}

void VarCut::setConstantValue(TString varName, float val){
  _constants[getConstantIndex(varName)] = numberNew(sig_figs, val);
}

float VarCut::getConstantValue(TString variable){
  return numberNew(sig_figs,  _constants[getConstantIndex(variable)]);
}


int VarCut::getVariableIndex(TString variable){
  for(int i=0; i<Vars::nVariables; i++){
    if(variable == Vars::variables[i]->name or variable == Vars::variables[i]->nameTmva) return i;
  }
  printf("VarCut::getVariableIndex: requested variable is not known!!!\n");
  exit(1);
}

int VarCut::getConstantIndex(TString constant){
  for(int i=0; i<Vars::nConstants; i++){
    if(constant == Vars::constants[i]->name) return i;
  }
  printf("VarCut::getConstantIndex: requested variable is not known!!!\n");
  exit(1);
}
bool VarCut::isSymmetric(TString variable){
  return Vars::variables[getVariableIndex(variable)]->symmetric;
}

// Print all cut values to stdout
void VarCut::printCuts(){
  printf("VarCut::print: Cut values are\n");
  for(int i=0; i<Vars::nVariables; i++){
    if(Vars::variables[i]->nameTmva == "hOverE" and _constants[getConstantIndex("C_rho")] > 0){
      printf("  %30s < %g + %g/E + %g*rho/E \n", Vars::variables[i]->nameTmva.Data(), numberNew(sig_figs, _cuts[i]), numberNew(sig_figs, _constants[getConstantIndex("C_E")]), numberNew(sig_figs, _constants[getConstantIndex("C_rho")]));
    } else if(Vars::variables[i]->nameTmva == "relIsoWithEA" and _constants[getConstantIndex("C_pt")] > 0){
      printf("  %30s < %g + %g/pt\n", Vars::variables[i]->nameTmva.Data(), numberNew(sig_figs, _cuts[i]), numberNew(sig_figs, _constants[getConstantIndex("C_pt")]));
    } else {
      printf("  %30s < %g\n", Vars::variables[i]->nameTmva.Data(), numberNew(sig_figs, _cuts[i]));
    }
  }
}
