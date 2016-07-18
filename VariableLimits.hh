#ifndef VARIABLELIMITS_HH
#define VARIABLELIMITS_HH

#include "TString.h"

namespace VarLims{

  const int nVarLimits = 8;
  
  struct VariableLimits {
    TString name;
    float   max;
    VariableLimits(TString nameIn, float maxIn):
      name(nameIn), max(maxIn){};
  };

  // For the veto working point, do not restrict anything
  VariableLimits * limitsNoRestrictions[nVarLimits] = {
    new VariableLimits("full5x5_sigmaIetaIeta"   ,1e30),
    new VariableLimits("dEtaSeed"                  ,1e30),
    new VariableLimits("dPhiIn"                  ,1e30),
    new VariableLimits("hOverE"                  ,1e30),
    new VariableLimits("relIsoWithEA"            ,1e30),
    new VariableLimits("ooEmooP"                 ,1e30),
    new VariableLimits("d0"                      ,1e30),
    new VariableLimits("dz"                      ,1e30)
    // new VariableLimits("expectedMissingInnerHits",1e30) // Removed from optimization
  };

  // For Loose, Medium, Tight working points:
  //   restrict d0 (at a very loose value)
  //   restrict missing hits, to Run 1-motivated value
  //      (to make sure ID discriminates well against 
  //       hard photons faking electrons, which may not be
  //       plentiful in the sample on which we are optimizing).
  //   restrict H/E, dPhiIn: all are very loose common sense
  //       restrictions.

  VariableLimits * limitsWPAnyV1[nVarLimits] = {
    new VariableLimits("full5x5_sigmaIetaIeta"   ,1e30),
    new VariableLimits("dEtaSeed"                ,0.03), // Endcap 99.9+%
    new VariableLimits("dPhiIn"                  ,0.27), // Endcap 99.9+%
    new VariableLimits("hOverE"                  ,0.37), // Barrel 99.9%
    new VariableLimits("relIsoWithEA"            ,1e30),
    new VariableLimits("ooEmooP"                 ,1e30),
    new VariableLimits("d0"                      ,0.40), // Endcap 99.9%
    new VariableLimits("dz"                      ,1e30)
    // new VariableLimits("expectedMissingInnerHits",1e30)  //NO Endcap 99.9%  // Removed from optimization
  };

  // For Medium and Tight working point, impose HLT-safe restrictions
  // See more details in:
  //   https://indico.cern.ch/event/491536/contributions/2202104/attachments/1288921/1918602/talk_electron_ID_spring16.pdf

  VariableLimits * limitsHLTSafeBarrel[nVarLimits] = {
    new VariableLimits("full5x5_sigmaIetaIeta"   ,0.011), // HLT WPLoose emulation
    new VariableLimits("dEtaSeed"                ,0.004), // HLT WPLoose emulation
    new VariableLimits("dPhiIn"                  ,1.000), // DISABLED HLT WPLoose emulation (default 0.020) because it won't be applied before ICHEP16   
    new VariableLimits("hOverE"                  ,0.06),  // HLT WPLoose emulation
    new VariableLimits("relIsoWithEA"            ,0.100), // UCCOM application
    new VariableLimits("ooEmooP"                 ,0.013), // HLT WPLoose emulation
    new VariableLimits("d0"                      ,0.40),  // Endcap 99.9%
    new VariableLimits("dz"                      ,1e30)  // No restriction
    // new VariableLimits("expectedMissingInnerHits",1e30)   // No restriction  // Removed from optimization
  };


  VariableLimits * limitsHLTSafeEndcap [nVarLimits] = {
    new VariableLimits("full5x5_sigmaIetaIeta"   ,0.031), // HLT WPLoose emulation
    new VariableLimits("dEtaSeed"                ,0.007), // HLT WPLoose emulation
    new VariableLimits("dPhiIn"                  ,1.000), // DISABLED HLT WPLoose emulation (default 0.020) because it won't be applied before ICHEP16
    new VariableLimits("hOverE"                  ,0.065), // HLT WPLoose emulation
    new VariableLimits("relIsoWithEA"            ,0.100), // UCCOM application
    new VariableLimits("ooEmooP"                 ,0.013),  // HLT WPLoose emulation
    new VariableLimits("d0"                      ,0.40),  // Endcap 99.9%
    new VariableLimits("dz"                      ,1e30)   // No restriction
    // new VariableLimits("expectedMissingInnerHits",2.05)   // HLT WPLoose emulation  // Removed from optimization
  };

}

#endif
