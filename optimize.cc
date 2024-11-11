#include "TSystem.h"

#include "optimize.hh"

// Manipulations with this is a guess, not documented in TMVA 
// at this point.
const TString datasetname = "dataset";

//
// Main method
//
void optimize(TString cutMaxFileName, TString cutsOutFileNameBase, TString trainingDataOutputBase, VarLims::VariableLimits **userDefinedCutLimits, int region){

  // region: 0->barrel ; 1->endcap; 2(else)->extend;

  TString fnameSignal     = Opt::fnameSignalBarrel;
  if (region == 0)     fnameSignal = Opt::fnameSignalBarrel;
  else if(region == 1) fnameSignal = Opt::fnameSignalEndcap;
  else                 fnameSignal = Opt::fnameSignalExtend;

  TString fnameBackground = Opt::fnameBackgroundBarrel;
  if (region == 0)     fnameBackground = Opt::fnameBackgroundBarrel;
  else if(region == 1) fnameBackground = Opt::fnameBackgroundEndcap;
  else                 fnameBackground = Opt::fnameBackgroundExtend;

  printf("\n Take true electrons from %s tree %s\n\n",       fnameSignal.Data(),     Opt::signalTreeName.Data());
  printf("\n Take background electrons from %s tree %s\n\n", fnameBackground.Data(), Opt::backgroundTreeName.Data());

  TTree *signalTree     = getTreeFromFile(fnameSignal,     Opt::signalTreeName,     &Opt::fileSignal);
  TTree *backgroundTree = getTreeFromFile(fnameBackground, Opt::backgroundTreeName, &Opt::fileBackground);

  // Configure output details
  TString trainingOutputDir = TString("trainingData/") + trainingDataOutputBase;
  printf("The directory where the xml results of the training is:\n");
  printf("         %s\n", trainingOutputDir.Data());

  FileStat_t buf;
  if( gSystem->GetPathInfo(trainingOutputDir.Data(), buf) ){
    printf("     this directory does not exist, creating it.\n");
    gSystem->MakeDirectory(trainingOutputDir.Data());
  }
  TMVA::gConfig().GetIONames().fWeightFileDir = trainingOutputDir;
  
  TString outfileName = trainingOutputDir + TString("/") + TString("TMVA_") + trainingDataOutputBase + TString(".root");
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
  printf("The ROOT file with train/test distributions from TMVA:\n");
  printf("         %s\n", outfileName.Data());

  // Create the factory object. Later you can choose the methods.
  // The factory is the only TMVA object you have to interact with
  //
  // The first argument is the base of the name of all the
  // weightfiles in the directory weight/
  //
  // The second argument is the output file for the training results
  // All TMVA output can be suppressed by removing the "!" (not) in
  // front of the "Silent" argument in the option string
  TString factoryOptions = "!V:!Silent:Color:DrawProgressBar:Transformations=I";
  TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile, factoryOptions);

  // Data loader handles trees, variables, etc
  TMVA::DataLoader *dataloader=new TMVA::DataLoader(datasetname);

  // Define the input variables that shall be used in the optimization
  configureVariables(dataloader);

  // Define per-tree weights and add trees to the data loader
  Double_t signalWeight     = 1.0;
  Double_t backgroundWeight = 1.0;
  dataloader->AddSignalTree    ( signalTree,     signalWeight     );
  dataloader->AddBackgroundTree( backgroundTree, backgroundWeight );
    
  // Set individual event weights (the variables must exist in the original TTree)
  // -  for signal    : `dataloader->SetSignalWeightExpression    ("weight1*weight2");`
  // -  for background: `dataloader->SetBackgroundWeightExpression("weight1*weight2");`
  dataloader->SetSignalWeightExpression("abs(genWeight*kinWeight)");
  dataloader->SetBackgroundWeightExpression("abs(genWeight*kinWeight)");

  // Configure training and test trees  
  TString trainAndTestOptions = getTrainAndTestOptions(region);

  // Apply additional cuts on the signal and background samples (can be different)
  TCut signalCuts = "";
  TCut backgroundCuts = "";
  configureCuts(signalCuts, backgroundCuts, region);

  // Tell the dataloader how to use the training and testing events
  dataloader->PrepareTrainingAndTestTree(signalCuts, backgroundCuts, trainAndTestOptions );
  
  // Book the Cuts method with the factory
  TString methodName = "Cuts";
  TString methodOptions = getMethodOptions(cutMaxFileName, userDefinedCutLimits);
  factory->BookMethod( dataloader, TMVA::Types::kCuts, methodName,methodOptions);
  
  // Do the work: optimization, testing, and evaluation
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();
  
  // Save working points into files.
  writeWorkingPoints(factory, cutsOutFileNameBase, region);

  // Clean up

  outputFile->Close();
  // When the factory is deleted, there appears to be a problem.
  // Things might crash or root does not exit. Commented it out until
  // better understanding.
  // delete factory;

  if(Opt::fileSignal != 0)     Opt::fileSignal->Close();
  if(Opt::fileBackground != 0) Opt::fileBackground->Close();

  delete factory;
  delete dataloader;

  return;
}

// Get a given tree from a given file name.
// Note: the **file is the way to return a pointer to a file
// back into the calling method.
TTree *getTreeFromFile(TString fname, TString tname, TFile **file){
  *file = new TFile( fname );
  TTree *tree = (TTree*) (*file)->Get(tname);
  return tree;
}

void configureVariables(TMVA::DataLoader *dataloader){

  // Define the input variables that shall be used for the MVA training
  // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
  // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]

  // Variables for cut optimization
  printf("Configure data loader variables for optimization\n");
  for(int i=0; i<Vars::nVariables; i++){
    TString varName = Vars::variables[i]->nameTmva;
    char varType = Vars::variables[i]->type;
    printf("    add variable %s of the type %c\n", varName.Data(), varType);
    dataloader->AddVariable( varName, varType );
  }
  
  // Spectator variables
  printf("Configure data loader spectator variables\n");
  for(int i=0; i<Vars::nSpectatorVariables; i++){
    TString varName = Vars::spectatorVariables[i]->nameTmva;
    char varType = Vars::spectatorVariables[i]->type;
    printf("    add spectator variable %s of the type %c\n", varName.Data(), varType);
    dataloader->AddSpectator( varName, varType );
  }
}

void configureCuts(TCut &signalCuts, TCut &backgroundCuts, int region){

  // Define all cuts 
 
  TCut etaCut = "";
  if(region == 0){
    printf("\n\nTraining for BARREL electrons\n\n");
    etaCut = Opt::etaCutBarrel;
  }
  else if(region == 1){
    printf("\n\nTraining for ENDCAP electrons\n\n");
    etaCut = Opt::etaCutEndcap;
  }
  else{
    printf("\n\nTraining for EXTEND electrons\n\n");
    etaCut = Opt::etaCutExtend;
  }
  TCut kinematicCuts = Opt::ptCut && etaCut;

  TCut preselectionCuts = kinematicCuts && Opt::otherPreselectionCuts;
  
  signalCuts = preselectionCuts && Opt::trueEleCut;
  backgroundCuts = preselectionCuts && Opt::fakeEleCut;  

}

TString getTrainAndTestOptions(int region){

  TString options = "SplitMode=Random:!V";
  options += ":nTrain_Signal=";
  //options += "2000";
  if(region == 0)      options += Opt::nTrain_SignalBarrel;
  else if(region == 1) options += Opt::nTrain_SignalEndcap;
  else                 options += Opt::nTrain_SignalExtend;
  options += ":nTrain_Background=";
  //options += "2000";
  if(region == 0)      options += Opt::nTrain_BackgroundBarrel;
  else if(region == 1) options += Opt::nTrain_BackgroundEndcap;
  else                 options += Opt::nTrain_BackgroundExtend;
  options += ":nTest_Signal=";
  //options += "2000";
  if(region == 0)      options += Opt::nTest_SignalBarrel;
  else if(region == 1) options += Opt::nTest_SignalEndcap;
  else                 options += Opt::nTest_SignalExtend;
  options += ":nTest_Background=";
  //options += "2000";
  if(region == 0)      options += Opt::nTest_BackgroundBarrel;
  else if(region == 1) options += Opt::nTest_BackgroundEndcap;
  else                 options += Opt::nTest_BackgroundExtend;
 
  printf("INFO: training and test options: %s\n", options.Data());
  return options;
}

TString getMethodOptions(TString cutMaxFileName, VarLims::VariableLimits **userDefinedCutLimits){

  TString methodOptions = Opt::methodCutsBaseOptions;

  // Next, put together cut-specific options
  TString cutsFileName = Opt::cutRepositoryDir;
  cutsFileName += "/";
  cutsFileName += cutMaxFileName;

  TFile *cutsFile = new TFile(cutsFileName);
  if( !cutsFile ) assert(0);
  VarCut *cutMax = (VarCut*)cutsFile->Get("cuts");
  if( !cutMax ) assert(0);

  if(!userDefinedCutLimits) assert(0);
  // Make sure the user defined cut limits array is consistent with the optimization
  // variables set
  bool checkPassed = true;
  if( Vars::nVariables != VarLims::nVarLimits ) checkPassed = false;
  for(int i=0; i<Vars::nVariables; i++){
    if( Vars::variables[i]->name != userDefinedCutLimits[i]->name )
      checkPassed = false;
  }
  if( !checkPassed ){
    printf("ERROR: the list of optimization variables is not consistent with the list\n");
    printf("       of the variables with user defined cut limits.\n");
    assert(0);
  }

  // As all cuts are upper cuts, we set the lower cut to -inf
  // Note: we do not have any negative vars, the vars that can be negative
  // are symmetric and enter as abs(XXX).
  for(int i=0; i<Vars::nVariables; i++){
    methodOptions += TString::Format(":VarProp[%d]=FMin",i);
  }
  // Add all cut ranges:
  for(int i=0; i<Vars::nVariables; i++){
    float max = cutMax->getCutValue(Vars::variables[i]->name);
    if( max > userDefinedCutLimits[i]->max ) max = userDefinedCutLimits[i]->max;
    methodOptions += TString::Format(":CutRangeMax[%d]=%.6f", i, max);
  }
  
  printf("\nMethod configuration: method options are\n");
  printf("%s\n", methodOptions.Data());
  printf("\n");

  return methodOptions;
}

void writeWorkingPoints(const TMVA::Factory *factory, TString cutsOutFileNameBase, int region){

  TString cutsFileName = Opt::cutRepositoryDir;
  cutsFileName += "/";
  cutsFileName += cutsOutFileNameBase;

  // Loop over four working points
  printf("The working points being saved:\n");
  for(int iwp=0; iwp<Opt::nWP; iwp++){
    TString cutsFileNameWP = cutsFileName;
    cutsFileNameWP += "_";
    cutsFileNameWP += Opt::wpNames[iwp];
    cutsFileNameWP += ".root";
    TFile *cutsFile = new TFile(cutsFileNameWP, "recreate");
    if( !cutsFile ) assert(0);
    VarCut *cutMax = new VarCut();
    
    const TMVA::MethodCuts *method = dynamic_cast<TMVA::MethodCuts*> (factory->GetMethod(datasetname,"Cuts"));
    if( method == 0 ) assert(0);

    std::vector <double> cutLo;
    std::vector <double> cutHi;
    if(region == 0)      method->GetCuts(Opt::effBarrel[iwp], cutLo, cutHi);
    else if(region == 1) method->GetCuts(Opt::effEndcap[iwp], cutLo, cutHi);
    else                 method->GetCuts(Opt::effExtend[iwp], cutLo, cutHi);
    // NOTE: this relies on filling the factory with AddVarilables
    // in exactly the same order (using the same loop) 
    // Start with a sanity check:
    if( Vars::nVariables != cutHi.size()) assert(0);
    // Now, fill the cut values into the storable object.
    for(uint ivar=0; ivar<cutHi.size(); ivar++){
      cutMax->setCutValue(Vars::variables[ivar]->name, cutHi.at(ivar));
    }
    printf("   working point %s\n", Opt::wpNames[iwp].Data());
    cutMax->printCuts();
    cutMax->Write("cuts");
    cutsFile->Close();
  }

}
