#! /usr/bin/env python

import ROOT,os,numpy
os.system('root -b -q VarCut.cc+ &> /dev/null')
ROOT.gSystem.Load('VarCut_cc')
ROOT.gROOT.SetBatch(True)

#
# Fill cuts into a VarCut object
#
def fillCuts(filename, dict):
  cutsFile   = ROOT.TFile('cut_repository/' + filename + '.root', 'recreate')
  cutsObject = ROOT.VarCut()
  for name, cut in dict.iteritems():
    cutsObject.setCutValue(name, cut) 
  cutsObject.Write("cuts")
  cutsFile.Close()

fillCuts('cuts_barrel_2016_WP_Veto',  {'full5x5_sigmaIetaIeta': 0.0115,  'dEtaSeed': 0.00749, 'dPhiIn': 0.228,  'hOverE': 0.356,  'relIsoWithEA': 0.175,  'ooEmooP': 0.299})
fillCuts('cuts_barrel_2016_WP_Veto',  {'full5x5_sigmaIetaIeta': 0.0115,  'dEtaSeed': 0.00749, 'dPhiIn': 0.228,  'hOverE': 0.356,  'relIsoWithEA': 0.175,  'ooEmooP': 0.299})
fillCuts('cuts_barrel_2016_WP_Loose', {'full5x5_sigmaIetaIeta': 0.011,   'dEtaSeed': 0.00477, 'dPhiIn': 0.222,  'hOverE': 0.298,  'relIsoWithEA': 0.0994, 'ooEmooP': 0.241})
fillCuts('cuts_barrel_2016_WP_Medium',{'full5x5_sigmaIetaIeta': 0.00998, 'dEtaSeed': 0.00311, 'dPhiIn': 0.103,  'hOverE': 0.253,  'relIsoWithEA': 0.0695, 'ooEmooP': 0.134})
fillCuts('cuts_barrel_2016_WP_Tight', {'full5x5_sigmaIetaIeta': 0.00998, 'dEtaSeed': 0.00308, 'dPhiIn': 0.0816, 'hOverE': 0.0414, 'relIsoWithEA': 0.0588, 'ooEmooP': 0.0129})
fillCuts('cuts_endcap_2016_WP_Veto',  {'full5x5_sigmaIetaIeta': 0.037,   'dEtaSeed': 0.00895, 'dPhiIn': 0.213,  'hOverE': 0.211,  'relIsoWithEA': 0.159,  'ooEmooP': 0.15})
fillCuts('cuts_endcap_2016_WP_Loose', {'full5x5_sigmaIetaIeta': 0.0314,  'dEtaSeed': 0.00868, 'dPhiIn': 0.213,  'hOverE': 0.101,  'relIsoWithEA': 0.107,  'ooEmooP': 0.14})
fillCuts('cuts_endcap_2016_WP_Medium',{'full5x5_sigmaIetaIeta': 0.0298,  'dEtaSeed': 0.00609, 'dPhiIn': 0.045,  'hOverE': 0.0878, 'relIsoWithEA': 0.0821, 'ooEmooP': 0.13})
fillCuts('cuts_endcap_2016_WP_Tight', {'full5x5_sigmaIetaIeta': 0.0292,  'dEtaSeed': 0.00605, 'dPhiIn': 0.0394, 'hOverE': 0.0641, 'relIsoWithEA': 0.0571, 'ooEmooP': 0.0129})
