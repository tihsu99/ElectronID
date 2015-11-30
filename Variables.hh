#ifndef CUTVARIABLES_H
#define CUTVARIABLES_H

#include "TString.h"

namespace Vars {

  const int nVariables = 9;
  const int nSpectatorVariables = 2;
  
  struct Variables{
    TString name;
    TString nameTmva; // for addition of abs
    char type;        // "F" for float, "I" for int
    bool symmetric;       // cuts symmetric or one-sided 
    // Constructor
    Variables(TString nameIn, TString nameTmvaIn, char typeIn, bool symIn):
      name(nameIn), nameTmva(nameTmvaIn), type(typeIn), symmetric(symIn){};
  };
  
  Variables *variables [nVariables] = {
    new Variables("full5x5_sigmaIetaIeta","full5x5_sigmaIetaIeta",'F',false),
    new Variables("dEtaIn"               ,"abs(dEtaIn)"          ,'F',true),
    new Variables("dPhiIn"               ,"abs(dPhiIn)"          ,'F',true),
    new Variables("hOverE"               ,"hOverE"               ,'F',false),
    new Variables("relIsoWithEA"      ,"relIsoWithEA"      ,'F',false),
    new Variables("ooEmooP"              ,"ooEmooP"              ,'F',false),
    new Variables("d0"                   ,"abs(d0)"              ,'F',true),
    new Variables("dz"                   ,"abs(dz)"              ,'F',true),
    new Variables("expectedMissingInnerHits","expectedMissingInnerHits",'I',false)
  };
  
  Variables *spectatorVariables [nSpectatorVariables] = {
    new Variables("pt"   , "pt"   , 'F',   false),
    new Variables("etaSC", "etaSC", 'F',   false)
  };
  
};

#endif
