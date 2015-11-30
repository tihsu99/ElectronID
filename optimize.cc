#include "TSystem.h"

#include "optimize.hh"

//
// Main method
//
void optimize(TString cutMaxFileName, TString cutsOutFileNameBase,
	      TString trainingDataOutputBase,
	      VarLims::VariableLimits **userDefinedCutLimits){

  // Input signal tree
  printf("\n Take true electrons from %s   tree %s\n\n", 
	 Opt::fnameSignal.Data(), Opt::signalTreeName.Data());
  TTree *signalTree = getTreeFromFile( Opt::fnameSignal, Opt::signalTreeName, 
				       &Opt::fileSignal );

  // Input background tree  
  printf("\n Take background electrons from %s   tree %s\n\n", 
	 Opt::fnameBackground.Data(), Opt::backgroundTreeName.Data());
  TTree *backgroundTree = getTreeFromFile( Opt::fnameBackground, 
					   Opt::backgroundTreeName,
					   &Opt::fileBackground);
  
  // Configure output details
  TString trainingOutputDir = TString("trainingData/")
    + trainingDataOutputBase;
  printf("The directory where the xml results of the training is:\n");
  printf("         %s\n", trainingOutputDir.Data());
  FileStat_t buf;
  if( gSystem->GetPathInfo(trainingOutputDir.Data(), buf) ){
    printf("     this directory does not exist, creating it.\n");
    gSystem->MakeDirectory(trainingOutputDir.Data());
  }
  TMVA::gConfig().GetIONames().fWeightFileDir = trainingOutputDir;
  
  TString outfileName = trainingOutputDir + TString("/")
    + TString("TMVA_") + trainingDataOutputBase
    + TString(".root");
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
  printf("The ROOT file with train/test distributions from TMVA:\n");
  printf("         %s\n", outfileName.Data());

  // Create factory
  TString factoryOptions = "!V:!Silent:Color:DrawProgressBar:Transformations=I";
  TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification",
					      outputFile, factoryOptions);
  configureFactoryVariables(factory);

  // Define weights and add trees to the factory
  Double_t signalWeight     = 1.0;
  Double_t backgroundWeight = 1.0;
  factory->AddSignalTree    ( signalTree,     signalWeight     );
  factory->AddBackgroundTree( backgroundTree, backgroundWeight );

  factory->SetBackgroundWeightExpression("genWeight");
  factory->SetSignalWeightExpression("genWeight");

  // Configure training and test trees  
  TString trainAndTestOptions = getTrainAndTestOptions();

  TCut signalCuts = "";
  TCut backgroundCuts = "";
  configureCuts(signalCuts, backgroundCuts);

  factory->PrepareTrainingAndTestTree( signalCuts, backgroundCuts,
				       trainAndTestOptions );
  
  // Book the Cuts method with the factory
  TString methodName = "Cuts";
  TString methodOptions = getMethodOptions(cutMaxFileName, userDefinedCutLimits);
  factory->BookMethod( TMVA::Types::kCuts, methodName,methodOptions);
  
  // Do the work: optimization, testing, and evaluation
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();
  
  // Save working points into files.
  writeWorkingPoints(factory, cutsOutFileNameBase);

  // Clean up

  outputFile->Close();
  // When the factory is deleted, there appears to be a problem.
  // Things might crash or root does not exit. Commented it out until
  // better understanding.
  // delete factory;

  if( Opt::fileSignal != 0 ){
    Opt::fileSignal->Close();
  }
  if( Opt::fileBackground != 0 ){
    Opt::fileBackground->Close();
  }

  return;
}

// Get a given tree from a given file name.
// Note: the **file is the way to return a pointer to a file
// back into the calling method.
TTree *getTreeFromFile(TString fname, TString tname, TFile **file){

  *file = new TFile( fname );
  TTree *tree     = (TTree*) (*file)->Get(tname);
  
  return tree;
}

void configureFactoryVariables(TMVA::Factory *factory){

  // Variables for cut optimization
  printf("Configure factory variables for optimization\n");
  for(int i=0; i<Vars::nVariables; i++){
    TString varName = Vars::variables[i]->nameTmva;
    char varType = Vars::variables[i]->type;
    printf("    add variable %s of the type %c\n", varName.Data(), varType);
    factory->AddVariable( varName, varType );
  }
  
  // Spectator variables
  printf("Configure factory spectator variables\n");
  for(int i=0; i<Vars::nSpectatorVariables; i++){
    TString varName = Vars::spectatorVariables[i]->nameTmva;
    char varType = Vars::spectatorVariables[i]->type;
    printf("    add spectator variable %s of the type %c\n", varName.Data(), varType);
    factory->AddSpectator( varName, varType );
  }
  
}

void configureCuts(TCut &signalCuts, TCut &backgroundCuts){

  // Define all cuts 
 
  TCut etaCut = "";
  if( Opt::useBarrel ){
    printf("\n\nTraining for BARREL electrons\n\n");
    etaCut = Opt::etaCutBarrel;
  }else{
    printf("\n\nTraining for ENDCAP electrons\n\n");
    etaCut = Opt::etaCutEndcap;
  }
  TCut kinematicCuts = Opt::ptCut && etaCut;

  TCut preselectionCuts = kinematicCuts && Opt::otherPreselectionCuts;
  
  signalCuts = preselectionCuts && Opt::trueEleCut;
  backgroundCuts = preselectionCuts && Opt::fakeEleCut;  

}

TString getTrainAndTestOptions(){

  TString options = "SplitMode=Random:!V";
  options += ":nTrain_Signal=";
  options += Opt::nTrain_Signal;
  options += ":nTrain_Background=";
  options += Opt::nTrain_Background;
  options += ":nTest_Signal=";
  options += Opt::nTest_Signal;
  options += ":nTest_Background=";
  options += Opt::nTest_Background;
 
  printf("INFO: training and test options: %s\n", options.Data());
  return options;
}

TString getMethodOptions(TString cutMaxFileName, 
			 VarLims::VariableLimits **userDefinedCutLimits){

  TString methodOptions = Opt::methodCutsBaseOptions;

  // Next, put together cut-specific options
  TString cutsFileName = Opt::cutRepositoryDir;
  cutsFileName += "/";
  cutsFileName += cutMaxFileName;

  TFile *cutsFile = new TFile(cutsFileName);
  if( !cutsFile )
    assert(0);
  VarCut *cutMax = (VarCut*)cutsFile->Get("cuts");
  if( !cutMax )
    assert(0);

  if(!userDefinedCutLimits)
    assert(0);
  // Make sure the user defined cut limits array is consistent with the optimization
  // variables set
  bool checkPassed = true;
  if( Vars::nVariables != VarLims::nVarLimits )
    checkPassed = false;
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
    if( max > userDefinedCutLimits[i]->max )
      max = userDefinedCutLimits[i]->max;
    methodOptions += TString::Format(":CutRangeMax[%d]=%.6f", i, max);
  }
  
  printf("\nMethod configuration: method options are\n");
  printf("%s\n", methodOptions.Data());
  printf("\n");

  return methodOptions;
}

void writeWorkingPoints(const TMVA::Factory *factory, TString cutsOutFileNameBase){

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
    if( !cutsFile )
      assert(0);
    VarCut *cutMax = new VarCut();
    
    const TMVA::MethodCuts *method = dynamic_cast<TMVA::MethodCuts*> (factory->GetMethod("Cuts"));
    if( method == 0 )
      assert(0);

    std::vector <double> cutLo;
    std::vector <double> cutHi;
    method->GetCuts(Opt::eff[iwp], cutLo, cutHi);
    // NOTE: this relies on filling the factory with AddVarilables
    // in exactly the same order (using the same loop) 
    // Start with a sanity check:
    if( Vars::nVariables != cutHi.size() )
      assert(0);
    // Now, fill the cut values into the storable object.
    for(uint ivar=0; ivar<cutHi.size(); ivar++){
      cutMax->setCutValue(Vars::variables[ivar]->name, 
			  cutHi.at(ivar));
    }
    printf("   working point %s\n", Opt::wpNames[iwp].Data());
    cutMax->print();
    cutMax->Write("cuts");
    cutsFile->Close();
  }

}

