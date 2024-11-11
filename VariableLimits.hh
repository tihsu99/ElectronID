#ifndef VARIABLELIMITS_HH
#define VARIABLELIMITS_HH

#include "TString.h"

namespace VarLims{

  const int nVarLimits = 6;
  
  struct VariableLimits {
    TString name;
    float   max;
    VariableLimits(TString nameIn, float maxIn):
      name(nameIn), max(maxIn){};
  };

  // For the veto working point, do not restrict anything
  VariableLimits * limitsNoRestrictions[nVarLimits] = {
    new VariableLimits("full5x5_sigmaIetaIeta"   ,1e30),
    new VariableLimits("dEtaSeed"                ,1e30),
    new VariableLimits("dPhiIn"                  ,1e30),
    new VariableLimits("hOverEscaled"            ,1e30),
    new VariableLimits("relIsoWithEA"            ,1e30),
    new VariableLimits("ooEmooP"                 ,1e30)
    //new VariableLimits("d0"                      ,1e30), // Removed from optimization
    //new VariableLimits("dz"                      ,1e30) // Removed from optimization
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
    new VariableLimits("dEtaSeed"                ,0.03), // Barrel/Endcap 99.9+%
    new VariableLimits("dPhiIn"                  ,0.25), // Barrel/Endcap 99.9+%
    new VariableLimits("hOverEscaled"                  ,0.35), // Barrel/Endcap 99,9%
    new VariableLimits("relIsoWithEA"            ,1.00), // Barrel/Endcap 99.9%
    new VariableLimits("ooEmooP"                 ,0.25)  // Barrel/Endcap 99.9%
    // new VariableLimits("d0"                      ,0.40), // Endcap 99.9%  // Removed from optimization
    // new VariableLimits("dz"                      ,1e30)  // Removed from optimization
    // new VariableLimits("expectedMissingInnerHits",1e30)  //NO Endcap 99.9%  // Removed from optimization
  };


  // Allow for possibility of HLT safe cuts for barrel and
  // endcap

  VariableLimits * limitsHLTSafeBarrel[nVarLimits] = {
    new VariableLimits("full5x5_sigmaIetaIeta"   ,1e30), // Disabled
    new VariableLimits("dEtaSeed"                ,1e30), // Disabled
    new VariableLimits("dPhiIn"                  ,1e30), // Disabled
    new VariableLimits("hOverEscaled"                  ,1e30), // Disabled
    new VariableLimits("relIsoWithEA"            ,1e30), // Disabled
    new VariableLimits("ooEmooP"                 ,1e30)  // Disabled
    // new VariableLimits("d0"                      ,1e30),  // Disabled
    // new VariableLimits("dz"                      ,1e30),  // Disabled
    // new VariableLimits("expectedMissingInnerHits",1e30)   // Disabled
  };


  VariableLimits * limitsHLTSafeEndcap [nVarLimits] = {
    new VariableLimits("full5x5_sigmaIetaIeta"   ,1e30), // 
    new VariableLimits("dEtaSeed"                ,1e30), // 
    new VariableLimits("dPhiIn"                  ,1e30), // 
    new VariableLimits("hOverEscaled"                  ,1e30), // 
    new VariableLimits("relIsoWithEA"            ,1e30), // 
    new VariableLimits("ooEmooP"                 ,1e30)  // 
    // new VariableLimits("d0"                      ,1e30),  // 
    // new VariableLimits("dz"                      ,1e30)   // 
    // new VariableLimits("expectedMissingInnerHits",1e30)   // 
  };

}

#endif
