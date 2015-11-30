#include "VarCut.hh"
#include <iomanip>
#include <iostream>
#include <sstream>

const  int sig_figs = 3;
  
auto numberNew =   [](int number_of_sig_figs, float number) 
{
  std::stringstream myStream;
  myStream << std::setprecision(number_of_sig_figs) << number;
  return std::stod(myStream.str());
};

const int UNDEFCUT = -999;

VarCut::VarCut()
{
  // All cuts are initialized to an unlikely value
  for(int i=0; i<Vars::nVariables; i++)
    _cuts[i]= UNDEFCUT;
};

// Construct TCut object for all cuts joined with &&
TCut *VarCut::getCut(){

  TCut *cut = 0;

  // Die if something appears uninitialized
  for(int i=0; i<Vars::nVariables; i++){
    if( _cuts[i] == UNDEFCUT ){
      printf("VarCut:: not all cuts are set! Die!\n");
      assert(0);
    }
  }

  cut = new TCut("");
  for(int i=0; i<Vars::nVariables; i++){
    // The += adds all cuts with &&:
    (*cut) += TString::Format(" %s < %f ",
			      Vars::variables[i]->nameTmva.Data(),
			      _cuts[i]);
  }
  
  return cut;
}

void VarCut::setCutValue(TString varName, float val){

  int index = getVariableIndex(varName);

  if( index != -1 ){
    _cuts[index] = numberNew(sig_figs, val);
  }else{
    printf("VarCut::setCutValue: requested variable is not known!!!\n");
  }

  return;
}

void VarCut::setCutValueTmvaName(TString varNameTmva, float val){

  int index = getVariableIndexTmvaName(varNameTmva);

  if( index != -1 ){
    _cuts[index] = numberNew(sig_figs, val);
  }else{
    printf("VarCut::setCutValue: requested variable is not known!!!\n");
  }

  return;
}

float VarCut::getCutValue(TString variable){

  float cutVal = UNDEFCUT;
  int index = getVariableIndex(variable);

  if( index != -1 ){
    cutVal = numberNew(sig_figs,  _cuts[index] );
  }else{
    printf("VarCut::getCutValue: requested variable is not known!!!\n");
  }
  if (variable == "expectedMissingInnerHits") cutVal = (int) cutVal;
  return cutVal;
}

int VarCut::getVariableIndex(TString variable){

  int index = -1;
  for(int i=0; i<Vars::nVariables; i++){
    if( variable == Vars::variables[i]->name ){
      index = i;
      break;
    } 
  }

  return index;
}

int VarCut::getVariableIndexTmvaName(TString variableTmva){

  int index = -1;
  for(int i=0; i<Vars::nVariables; i++){
    if( variableTmva == Vars::variables[i]->nameTmva ){
      index = i;
      break;
    } 
  }

  return index;
}

bool VarCut::isSymmetric(TString variable){

  bool result = false;
  int index = getVariableIndex(variable);

  if( index != -1 ){
    result = Vars::variables[index]->symmetric;
  }else{
    printf("VarCut::isSymmetric: requestd variable is not known!!!\n");
  }

  return result;
}

// Print all cut values to stdout
void VarCut::print(){

  printf("VarCut::print: Cut values are\n");
  for(int i=0; i<Vars::nVariables; i++){
    if (i == 8)
      printf("  %30s <= %g\n", Vars::variables[i]->nameTmva.Data(), numberNew(sig_figs, _cuts[i]) );
    else 
      printf("  %30s < %g\n", Vars::variables[i]->nameTmva.Data(), numberNew(sig_figs, _cuts[i]) );
    
  
  }

}


  
