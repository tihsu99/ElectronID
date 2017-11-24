#! /usr/bin/env python

import ROOT,os
from common import loadClasses, workingPoints, getTreeFromFile, drawFromTree
loadClasses('VarCut.cc', 'OptimizationConstants.hh')


def findSignalEfficiency(signalTree, cuts, wp, barrel, missingHits, tuneC0=None, C_E=None, C_rho=None):
  preselectionCuts = ROOT.TCut(ROOT.Opt.ptCut) + (ROOT.Opt.etaCutBarrel if barrel else ROOT.Opt.etaCutEndcap) + ROOT.Opt.otherPreselectionCuts;
  signalCuts       = preselectionCuts + ROOT.Opt.trueEleCut;
  
  if tuneC0:
    selectionCuts  = ROOT.TCut('&&'.join(['(' + var + '<' + str(cuts.getCutValue(var)) + ')' for var in ["full5x5_sigmaIetaIeta", "dEtaSeed", "dPhiIn", "relIsoWithEA", "ooEmooP"]]))
    selectionCuts += ROOT.TCut('hOverE<' + str(tuneC0) + '+' + str(C_E) + '/eSC+' + str(C_rho) + '*rho/eSC')
  else:
    selectionCuts  = cuts.getCut()

  if missingHits:
    selectionCuts += ROOT.TCut('expectedMissingInnerHits<='+str(wp.missingHitsBarrel if barrel else wp.missingHitsEndcap))

  hS_num = drawFromTree(signalTree, signalCuts + selectionCuts, 'pt', 'hS_num', [0.,10000.], wp)
  hS_den = drawFromTree(signalTree, signalCuts,                 'pt', 'hS_den', [0.,10000.], wp)
  return hS_num.GetSumOfWeights()/hS_den.GetSumOfWeights()



def tuneC0(region, tag):
  C_E   = 1.63   if region=='barrel' else 2.65
  C_rho = 0.0368 if region=='barrel' else 0.201

  signalTree = getTreeFromFile('2017-11-16/DYJetsToLL_flat_ntuple_true_' + region + '_full.root', ROOT.Opt.signalTreeName)

  for wp in workingPoints[tag]:
    fileName = os.path.join('cut_repository', wp.cutsFileBarrel if region=='barrel' else wp.cutsFileEndcap)
    file     = ROOT.TFile(fileName)
    cuts     = file.Get('cuts')

    effSignal = findSignalEfficiency(signalTree, cuts, wp, region=='barrel', True)
    print wp.name + ' --> tuning for ' + str(effSignal)

    C_0  = 0.05
    step = 0.02
    while True:
      effTemp     = findSignalEfficiency(signalTree, cuts, wp, region=='barrel', True, C_0, C_E, C_rho)
      print 'Using  ' + str(C_0) + ' --> eff: ' + str(effTemp)
      if   abs(effTemp - effSignal) < 0.001: break
      if effTemp < effSignal:
        C_0  += step
        step *= 2./3.
      C_0 -= step

    cuts.setCutValue('hOverE', min(C_0,0.05))
    cuts.printCuts()

    file = ROOT.TFile(fileName.replace('WP', 'retuned_WP'), 'recreate')
    cuts.Write("cuts")
    file.Close()

for region in ['barrel','endcap']:
  tuneC0(region, 'default')
