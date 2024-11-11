#! /usr/bin/env python

import ROOT,os
from common import loadClasses, workingPoints, getTreeFromFile, drawFromTree, makeSubDirs, setColors, HistoryData_barrel, HistoryData_endcap
loadClasses('VarCut.cc', 'OptimizationConstants.hh')

def tmvaFileName(wp, region, tag):
  if 'Veto'   in wp.name: wpPass = 'pass1'
  if 'Loose'  in wp.name: wpPass = 'pass2'
  if 'Medium' in wp.name: wpPass = 'pass3'
  if 'Tight'  in wp.name: wpPass = 'pass4'
  if   tag=='default':     name = 'training_results_' + region + '_' + wpPass + '_2017-11-07'
  elif tag=='retuneMVA':   name = 'training_results_' + region + '_' + wpPass + '_2017-11-16'
  elif tag=='training94':  name = 'training_results_' + region + '_' + wpPass + '_2018-03-18'
  elif tag=='training106': name = 'training_results_' + region + '_' + wpPass + '_2019-08-23'
  elif tag=='training106_tuneC0': name = 'training_results_' + region + '_' + wpPass + '_2019-08-23'
  elif tag=='training106_tuneC0CE_v3': name = 'training_results_' + region + '_' + wpPass + '_2019-08-23'
  elif tag=='training106_tuneCE_vs': name = 'training_results_' + region + '_' + wpPass + '_2019-08-23'
  elif tag=='training106_tuneC0CE_v5': name = 'training_results_' + region + '_' + wpPass + '_2019-08-23'
  return './trainingData/' + name + '/TMVA_' + name + '.root'

def findEfficiencies(signalTree, backgroundTree, cuts, wp, region, missingHits):
  if(region == 'barrel'):
    preselectionCuts = ROOT.TCut(ROOT.Opt.ptCut) + ROOT.Opt.etaCutBarrel + ROOT.Opt.otherPreselectionCuts;
  elif(region == 'endcap'):
    preselectionCuts = ROOT.TCut(ROOT.Opt.ptCut) + ROOT.Opt.etaCutEndcap + ROOT.Opt.otherPreselectionCuts;
  else:
    preselectionCuts = ROOT.TCut(ROOT.Opt.ptCut) + ROOT.Opt.etaCutExtend + ROOT.Opt.otherPreselectionCuts;

  signalCuts       = preselectionCuts + ROOT.Opt.trueEleCut;
  backgroundCuts   = preselectionCuts + ROOT.Opt.fakeEleCut;  
  selectionCuts    = ROOT.TCut(cuts.getCut())
 
  if missingHits:
    if(region == 'barrel'):
      selectionCuts += ROOT.TCut('expectedMissingInnerHits<='+str(wp.missingHitsBarrel))
    elif(region == 'endcap'):
      selectionCuts += ROOT.TCut('expectedMissingInnerHits<='+str(wp.missingHitsEndcap))
    else:
      selectionCuts += ROOT.TCut('expectedMissingInnerHits<='+str(wp.missingHitsExtend))

  cuts.printCuts()
  hS_num  = drawFromTree(signalTree,     signalCuts + selectionCuts,     'pt', 'hS_num',  [0.,10000.], wp)
  hS_den  = drawFromTree(signalTree,     signalCuts,                     'pt', 'hS_den',  [0.,10000.], wp)
  hBG_num = drawFromTree(backgroundTree, backgroundCuts + selectionCuts, 'pt', 'hBG_num', [0.,10000.], wp)
  hBG_den = drawFromTree(backgroundTree, backgroundCuts,                 'pt', 'hBG_den', [0.,10000.], wp)

  return hS_num.GetSumOfWeights()/hS_den.GetSumOfWeights(), hBG_num.GetSumOfWeights()/hBG_den.GetSumOfWeights()



def drawROCandWP(region, missingHits, tag):
  c1 = ROOT.TCanvas("c1","",10,10,600,600)

  hROC     = ROOT.TH1F("hROC", "", 100, 0, 1)
  hROC.min = 1

  if(region == 'barrel'):
    signalTree     = getTreeFromFile(ROOT.Opt.fnameSignalBarrel,     ROOT.Opt.signalTreeName)
    backgroundTree = getTreeFromFile(ROOT.Opt.fnameBackgroundBarrel, ROOT.Opt.backgroundTreeName)
  elif(region == 'endcap'):
    signalTree     = getTreeFromFile(ROOT.Opt.fnameSignalEndcap,     ROOT.Opt.signalTreeName)
    backgroundTree = getTreeFromFile(ROOT.Opt.fnameBackgroundEndcap, ROOT.Opt.backgroundTreeName)
  else:
    signalTree     = getTreeFromFile(ROOT.Opt.fnameSignalExtend,     ROOT.Opt.signalTreeName)
    backgroundTree = getTreeFromFile(ROOT.Opt.fnameBackgroundExtend, ROOT.Opt.backgroundTreeName)

  leg = ROOT.TLegend(0.15, 0.45, 0.5, 0.7)
  leg.SetFillStyle(0)
  leg.SetBorderSize(0)

  markers = {}
  markers_historical = {}
  setColors(workingPoints[tag])
  for wp in reversed(workingPoints[tag]):
    if(region == 'barrel'):
      file = ROOT.TFile('cut_repository/' + (wp.cutsFileBarrel) + '.root')
    elif(region == 'endcap'):
      file = ROOT.TFile('cut_repository/' + (wp.cutsFileEndcap) + '.root')
    else:
      file = ROOT.TFile('cut_repository/' + (wp.cutsFileExtend) + '.root')
    cuts = file.Get('cuts')

    effSignal, effBackground = findEfficiencies(signalTree, backgroundTree, cuts, wp, region, missingHits);
    
    c1.cd()
    markers[wp] = ROOT.TMarker(effSignal, 1.0-effBackground, 20)
    markers[wp].SetMarkerSize(2)
    markers[wp].SetMarkerColor(wp.sColor)

    leg.AddEntry(markers[wp], wp.name.replace('106','122'), "p")

    if(region == 'barrel'):
      history_data = HistoryData_barrel[wp.name.replace('106','94')]
      markers_historical[wp] = ROOT.TMarker(history_data[0], history_data[1], 20)
      markers_historical[wp].SetMarkerSize(2)
      markers_historical[wp].SetMarkerColor(history_data[2])
      leg.AddEntry(markers_historical[wp], wp.name.replace('106','94'), "p")
    elif(region == 'endcap'):
      history_data = HistoryData_endcap[wp.name.replace('106','94')]
      markers_historical[wp] = ROOT.TMarker(history_data[0], history_data[1], 20)
      markers_historical[wp].SetMarkerSize(2)
      markers_historical[wp].SetMarkerColor(history_data[2])
      leg.AddEntry(markers_historical[wp], wp.name.replace('106','94'), "p")


    tmvaFile = ROOT.TFile(tmvaFileName(wp, region, tag))
    wpROC    = tmvaFile.Get("dataset/Method_Cuts/Cuts/MVA_Cuts_rejBvsS")
    assert tmvaFile, hROC

    for i in range(hROC.min, 101):
      if wpROC.GetBinContent(i) < 0.5:
        hROC.min = i
        break
      hROC.SetBinContent(i, wpROC.GetBinContent(i))

  hROC.SetStats(0)
  hROC.SetLineWidth(2)
  hROC.SetTitle("")
  hROC.GetXaxis().SetTitle("signal efficiency")
  hROC.GetYaxis().SetTitle("background rejection")
  hROC.GetYaxis().SetTitleOffset(1.4)
  hROC.GetXaxis().SetRangeUser(0.6, 1.0)
  hROC.GetYaxis().SetRangeUser(0.931 if region=='barrel' else 0.8, 1.0)

  c1.cd()
  hROC.Draw("L")
  for m in markers.values(): m.Draw("same")
  for m in markers_historical.values(): m.Draw("same")
  leg.Draw()

  comment = ROOT.TLatex(0.2, 0.2, region + ' electrons')
  comment.SetNDC(True)
  comment.Draw()

  dirName  = os.path.join('figures', tag if tag!='default' else '')
  fileName = os.path.join(dirName, 'plot_ROCandWP_' + region + ('' if missingHits else '_noMissingHits') + '.png')
  c1.Print(makeSubDirs(fileName))


for tag in ['training106_tuneC0CE_v5']:
  for region in ['endcap','barrel','extend']:
    for missingHits in [True, False]:
      drawROCandWP(region, missingHits, tag)
