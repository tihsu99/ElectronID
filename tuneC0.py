#! /usr/bin/env python

import ROOT,os
from common import loadClasses, workingPoints, getTreeFromFile, drawFromTree
loadClasses('VarCut.cc', 'OptimizationConstants.hh')


def findSignalEfficiency(signalTree, cuts, wp, barrel, missingHits, trueEle = True, range=[0.,10000.]):
  preselectionCuts = ROOT.TCut(ROOT.Opt.ptCut) + (ROOT.Opt.etaCutBarrel if barrel else ROOT.Opt.etaCutEndcap) + ROOT.Opt.otherPreselectionCuts;
  signalCuts       = (preselectionCuts + ROOT.Opt.trueEleCut) if trueEle else preselectionCuts;
  selectionCuts    = cuts.getCut()

  if missingHits:
    selectionCuts += ROOT.TCut('expectedMissingInnerHits<='+str(wp.missingHitsBarrel if barrel else wp.missingHitsEndcap))

  hS_num = drawFromTree(signalTree, signalCuts + selectionCuts, 'pt', 'hS_num', [0.,10000.], wp)
  hS_den = drawFromTree(signalTree, signalCuts,                 'pt', 'hS_den', [0.,10000.], wp)
  return hS_num.GetSumOfWeights()/hS_den.GetSumOfWeights()



def tuneC0(region, tag):
  C_E   = 1.16   if region=='barrel' else 2.54
  C_rho = 0.0324 if region=='barrel' else 0.183
  C_pt  = 0.506  if region=='barrel' else 0.963

  signalTree = getTreeFromFile('2018-03-18/DYJetsToLL_flat_ntuple_true_' + region + '_full.root', ROOT.Opt.signalTreeName)
  highPtTree = getTreeFromFile('2018-03-18/DoubleEleFlat_flat_ntuple_trueAndFake_alleta_full.root', ROOT.Opt.signalTreeName)

  for wp in workingPoints[tag]:
    fileName = os.path.join('cut_repository', (wp.cutsFileBarrel if region=='barrel' else wp.cutsFileEndcap) + '.root')
    file     = ROOT.TFile(fileName)
    cuts     = file.Get('cuts')

    effSignal = findSignalEfficiency(signalTree, cuts, wp, region=='barrel', True)
    print wp.name + ' --> tuning for ' + str(effSignal)

    # HoverE tuning
    cuts.setConstantValue('C_E',   C_E)
    cuts.setConstantValue('C_rho', C_rho)
    C_0  = 0.05
    step = 0.02
    while True:
      cuts.setCutValue('hOverE', C_0)
      effTemp   = findSignalEfficiency(signalTree, cuts, wp, region=='barrel', True)
      effHighPt = findSignalEfficiency(highPtTree, cuts, wp, region=='barrel', True, trueEle=False, range=[1000.,2000.])
      print 'hOverE tuning with C_0=' + str(C_0) + ' --> eff: ' + str(effTemp)
      if abs(effTemp - effSignal) < 0.0005: break
      if effTemp > effSignal and effHighPt - 0.10 < effSignal:
        print 'eff at 1-2 TeV dropped 0.10 below target efficiency, stop here to avoid degrading efficiency at high pt'
        break # make sure that the efficiency does not drop too much at high pt
      if C_0 > 0.06: break
      if effTemp < effSignal:
        C_0  += step
        step *= 2./3.
      C_0 -= step

    cuts.setCutValue('hOverE', min(C_0,0.05))

    # relIso tuning
    cuts.setConstantValue('C_pt', C_pt)
    C_0  = cuts.getCutValue('relIsoWithEA')
    step = 0.02
    effTemp = effSignal
    while True:
      cuts.setCutValue('relIsoWithEA', C_0)
      effLast = effTemp
      effTemp = findSignalEfficiency(signalTree, cuts, wp, region=='barrel', True)
      print 'relIso tuning with C_0=' + str(C_0) + ' --> eff: ' + str(effTemp)
      if abs(effTemp - effSignal) < 0.0005: break
      if abs(effTemp - effLast) < 0.000001: 
        print 'To few difference, giving up'
        break
      if effTemp < effSignal:
        C_0  += step
        step *= 2./3.
      C_0 -= step

    cuts.printCuts()

    tag = 'retuned_WP3'
    file = ROOT.TFile(fileName.replace('WP', tag), 'recreate')
    cuts.Write("cuts")
    file.Close()

for region in ['barrel','endcap']:
  tuneC0(region, 'training94')
