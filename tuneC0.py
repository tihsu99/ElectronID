#! /usr/bin/env python

import ROOT,os
from common import loadClasses, workingPoints, getTreeFromFile, drawFromTree
loadClasses('VarCut.cc', 'OptimizationConstants.hh')


def findSignalEfficiency(signalTree, cuts, wp, region, missingHits, trueEle = True, range=[0.,10000.]):

  preselectionCuts = ROOT.TCut(ROOT.Opt.ptCut)
  if(region == 'barrel'):
    preselectionCuts = ROOT.TCut(ROOT.Opt.ptCut) + ROOT.Opt.etaCutBarrel + ROOT.Opt.otherPreselectionCuts;
  elif(region == 'endcap'):
    preselectionCuts = ROOT.TCut(ROOT.Opt.ptCut) + ROOT.Opt.etaCutEndcap + ROOT.Opt.otherPreselectionCuts;
  else:
    preselectionCuts = ROOT.TCut(ROOT.Opt.ptCut) + ROOT.Opt.etaCutExtend + ROOT.Opt.otherPreselectionCuts;

  signalCuts       = (preselectionCuts + ROOT.Opt.trueEleCut) if trueEle else preselectionCuts;
  selectionCuts    = cuts.getCut()

  if missingHits:
    if(region == 'barrel'):
      selectionCuts += ROOT.TCut('expectedMissingInnerHits<='+str(wp.missingHitsBarrel))
    elif(region == 'endcap'):
      selectionCuts += ROOT.TCut('expectedMissingInnerHits<='+str(wp.missingHitsEndcap))
    else:
      selectionCuts += ROOT.TCut('expectedMissingInnerHits<='+str(wp.missingHitsExtend))

  hS_num = drawFromTree(signalTree, signalCuts + selectionCuts, 'pt', 'hS_num', [0.,10000.], wp)
  hS_den = drawFromTree(signalTree, signalCuts,                 'pt', 'hS_den', [0.,10000.], wp)
  return hS_num.GetSumOfWeights()/hS_den.GetSumOfWeights()



def tuneC0(region, tag):
  C_E = 0.0
  C_rho = 0.0
  C_pt = 0.0
  if(region == 'barrel'):
    C_E = 1.39
    C_rho = 0.0422
    C_pt = 1.33
  elif(region == 'endcap'):
    C_E = 2.25
    C_rho = 0.266
    C_pt = 1.19
  else:
    C_E = 1.44
    C_rho = 0.697
    C_pt = 0.956
  signalTree = getTreeFromFile('2019-08-23/DY_ext_flat_ntuple_true_' + region + '_full.root', ROOT.Opt.signalTreeName)
  highPtTree = getTreeFromFile('2019-08-23/DoubleEleFlat_flat_ntuple_trueAndFake_alleta_full.root', ROOT.Opt.signalTreeName)

  for wp in workingPoints[tag]:

    if(region == 'barrel'):
      fileName = os.path.join('cut_repository', wp.cutsFileBarrel + '.root')
    elif(region == 'endcap'):
      fileName = os.path.join('cut_repository', wp.cutsFileEndcap + '.root')
    else:
      fileName = os.path.join('cut_repository', wp.cutsFileExtend + '.root')

    file     = ROOT.TFile(fileName)
    cuts     = file.Get('cuts')

    effSignal = findSignalEfficiency(signalTree, cuts, wp, region, True)
    print wp.name + ' --> tuning for ' + str(effSignal)

    # HoverE tuning
    cuts.setConstantValue('C_E',   C_E)
    cuts.setConstantValue('C_rho', C_rho)
    C_0  = 0.05
    step = 0.02
    while True:
      cuts.setCutValue('hOverE', C_0)
      effTemp   = findSignalEfficiency(signalTree, cuts, wp, region, True)
#      effTemp_100To200 = findSignalEfficiency(highPtTree, cuts, wp, region, True, trueEle=False, range=[100.,200.])
      effHighPt = findSignalEfficiency(highPtTree, cuts, wp, region, True, trueEle=False, range=[1000.,2000.])

      print 'hOverE tuning with C_0=' + str(C_0) + ' --> eff: ' + str(effTemp)
      if ((abs(effTemp - effSignal) < 0.0005) and not ((effHighPt < effTemp - 0.05 and C_E > 0.5))): break

      if (((abs(effTemp - effSignal) < 0.0005) or (effTemp > effSignal)) and ((effHighPt < effTemp - 0.05) and C_E > 0.5)): # extend even can not achieve this at Veto WP
        C_E -= 0.01
        C_E = max(C_E,0.5)
        print 'eff at 1-2 TeV dropped 0.05 below peak efficiency, tune CE to retrain again.'
        print 'eff at 1-2 TeV=' + str(effHighPt) + ' eff(target)=' + str(effTemp)
        print 'C_E=' + str(C_E)
        cuts.setConstantValue('C_E',C_E)
        C_0  = 0.05
        step = 0.02
        continue
#        break # make sure that the efficiency does not drop too much at high pt
      if C_0 > 0.06: break
      if effTemp < effSignal:
        C_0  += step
        step *= 2./3.
      C_0 -= step

    cuts.setCutValue('hOverE', max(min(C_0,0.05),0.02))

    # relIso tuning
    cuts.setConstantValue('C_pt', C_pt)
    C_0  = cuts.getCutValue('relIsoWithEA')
    step = 0.02
    effTemp = effSignal
    while True:
      cuts.setCutValue('relIsoWithEA', C_0)
      effLast = effTemp
      effTemp = findSignalEfficiency(signalTree, cuts, wp, region, True)
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

    tag = 'training106_tuneC0CE_v4'
    file2 = ROOT.TFile(fileName.replace('WP', tag), 'recreate')
    file2.cd()
    cuts.Write("cuts")
    file2.Close()

for region in ['barrel','endcap','extend']:
  tuneC0(region, 'training106')
