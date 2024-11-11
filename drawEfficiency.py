#! /usr/bin/env python

import ROOT,os,numpy,shutil
from common import loadClasses, workingPoints, getTreeFromFile, drawFromTree, getCuts, makeSubDirs, setColors
loadClasses('VarCut.cc', 'OptimizationConstants.hh')

dateTag = "2019-08-23"


#
# Helper functions for trees and binning
#
def getPtBins(longPtRange):
  bins = [10.]
  for i in range(30): bins.append(bins[-1] + 1)
  for i in range(30): bins.append(bins[-1] + 2)
  for i in range(10): bins.append(bins[-1] + 5)
  for i in range(5):  bins.append(bins[-1] + 10)
  if(longPtRange):
    for i in range(20):  bins.append(bins[-1] + 20)
    for i in range(10):  bins.append(bins[-1] + 40)
    for i in range(10):  bins.append(bins[-1] + 100)
  return bins;

def getEtaBins():
  nBins = 60
  xmin  = -3.0
  xmax  = 3.0
  delta = (xmax - xmin)/nBins

  bins = []
  for i in range(0, nBins+1): bins.append(xmin + i*delta);
  return bins

def getNvtxBins():
  return [0., 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 27, 30, 35, 40, 45, 50, 60, 70, 90]

#
# Calculating the efficiency histogram with correct errors
#
def calculateEffAndErrors(numHist, denHist, eName, binning):
  ehist = ROOT.TH1F(eName, eName, len(binning)-1, numpy.array(binning))
  ehist.Sumw2()

  for i in range(1, len(binning)):
    eff = 0.0;
    effErr2 = 0.0;
    effErr = 0.0;
    num    = numHist.GetBinContent(i);
    numErr = numHist.GetBinError(i);
    den    = denHist.GetBinContent(i);
    denErr = denHist.GetBinError(i);
    tot    = num+den;

    eff     = num/den if den > 0 else 0
    effErr2 = (num*num*den*den)/(tot*tot*tot*tot)*(numErr*numErr/(num*num)+denErr*denErr/(den*den)) if den > 0 and num > 0 else 0
    effErr = numpy.sqrt(effErr2 if effErr2 > 0 else 0)

    ehist.SetBinContent(i, eff)
    ehist.SetBinError(i, effErr)

  return ehist;

#
# Histogram style
#
def setHistogram(hist, wp, isSignal):
  hist.SetLineColor(wp.sColor if isSignal else wp.bgColor)
  hist.SetLineWidth(2)

  hist.SetMarkerStyle(20)
  hist.SetMarkerSize(0.8)
  hist.SetMarkerColor(wp.sColor if isSignal else wp.bgColor)

#
# Main function
#
def drawEfficiency(mode, tag, region, selectVar):
  if '2TeV' in mode:
#    return # not available
    signalFileName = dateTag + '/' + "DoubleEleFlat_flat_ntuple_trueAndFake_alleta_full.root"
  elif '1to500GeV' in mode:
    signalFileName = dateTag + '/' + "DoubleEle1to500_flat_ntuple_trueAndFake_alleta_full.root"
  elif '500to1000GeV' in mode:
    signalFileName = dateTag + '/' + "DoubleEle500to1000_flat_ntuple_trueAndFake_alleta_full.root"
  else:              signalFileName = dateTag + '/' + "DY_ext_flat_ntuple_true_alleta_full.root"
  signalTree     = getTreeFromFile(signalFileName,                                                ROOT.Opt.signalTreeName)
  backgroundTree = getTreeFromFile(dateTag + '/' + "TT_flat_ntuple_trueAndFake_alleta_full.root", ROOT.Opt.backgroundTreeName) #if '2TeV' not in mode else None

  if('genPt' in mode):  binning = getPtBins(mode.count("eV"))
  elif('pt' in mode):   binning = getPtBins(mode.count("eV"))
  elif('eta' in mode):  binning = getEtaBins()
  elif('nvtx' in mode): binning = getNvtxBins()

  if('genPt' in mode):  axisLabel = "generator p_{T} [GeV]"
  elif('pt' in mode):   axisLabel = "p_{T} [GeV]"
  elif('eta' in mode):  axisLabel = "#eta_{SC}"
  elif('nvtx' in mode): axisLabel = "Nvtx"

  if(mode == "eta"): comment = '';
  else:              comment = region + ' electrons'

  preselectionCuts = ROOT.Opt.ptCut + ROOT.Opt.otherPreselectionCuts
  if(mode == "eta"):        preselectionCuts += ROOT.TCut('abs(etaSC)<3.0')
  elif(region == "barrel"): preselectionCuts += ROOT.Opt.etaCutBarrel
  elif(region == "endcap"): preselectionCuts += ROOT.Opt.etaCutEndcap
  else:                     preselectionCuts += ROOT.Opt.etaCutExtend

  signalCuts     = preselectionCuts + (ROOT.Opt.trueEleCut if not 'DoubleEle' in signalFileName else ROOT.TCut())
  backgroundCuts = preselectionCuts + ROOT.Opt.fakeEleCut

  c1 = ROOT.TCanvas("c1","c1",10,10,600,600)
  ROOT.gStyle.SetOptStat(0)
  c1.cd()

  dummy = ROOT.TH2D("dummy","", 100, binning[0], binning[-1], 100, 0, 1);
  dummy.GetXaxis().SetTitle(axisLabel);
  dummy.GetXaxis().SetTitleOffset(1.2);
  dummy.GetYaxis().SetTitle("efficiency");
  dummy.GetYaxis().SetTitleOffset(1.2);
  dummy.Draw();

  sigEff = {}
  bgEff = {}
  setColors(workingPoints[tag])
  for wp in workingPoints[tag]:
    is2016 = '2016' in wp.name
    wpCutsBarrel = getCuts(wp, 'barrel', selectVar)
    wpCutsEndcap = getCuts(wp, 'endcap', selectVar)
    wpCutsExtend = getCuts(wp, 'extend', selectVar)

    selectionCutsBarrel = ROOT.TCut(wpCutsBarrel) + ROOT.Opt.etaCutBarrel
    selectionCutsEndcap = ROOT.TCut(wpCutsEndcap) + ROOT.Opt.etaCutEndcap
    selectionCutsExtend = ROOT.TCut(wpCutsExtend) + ROOT.Opt.etaCutExtend
    selectionCuts       = ROOT.TCut('(' + selectionCutsBarrel.GetTitle() + ')||(' + selectionCutsEndcap.GetTitle() + ')||(' + selectionCutsExtend.GetTitle() + ')')

    if('genPt' in mode):  varName = "genPt"
    elif('pt' in mode):   varName = "pt"
    elif('eta' in mode):  varName = "etaSC"
    elif('nvtx' in mode): varName = "nPV"

    sigNum = drawFromTree(signalTree,     (signalCuts + selectionCuts),     varName, 'sigNum'+ wp.name.split()[0], binning, wp)
    sigDen = drawFromTree(signalTree,     (signalCuts),                     varName, 'sigDen'+ wp.name.split()[0], binning, wp)
    bgNum  = drawFromTree(backgroundTree, (backgroundCuts + selectionCuts), varName, 'bgNum' + wp.name.split()[0], binning, wp) if backgroundTree else 1
    bgDen  = drawFromTree(backgroundTree, (backgroundCuts),                 varName, 'bgDen' + wp.name.split()[0], binning, wp) if backgroundTree else 1

    sigEff[wp] = calculateEffAndErrors(sigNum, sigDen, 'sigEff' + wp.name, binning)
    bgEff[wp]  = calculateEffAndErrors(bgNum,  bgDen,  'bgEff'  + wp.name, binning) if backgroundTree else None

    setHistogram(sigEff[wp], wp, True)
    sigEff[wp].Draw("same,pe")

    if backgroundTree:
      setHistogram( bgEff[wp], wp, False)
      bgEff[wp].Draw("same,pe")
    c1.Update()

  if('pt' in mode or 'Pt' in mode): leg = ROOT.TLegend(0.2, 0.3, 0.6, 0.7)
  elif('eta' in mode):              leg = ROOT.TLegend(0.38, 0.2, 0.75, 0.6)
  elif('nvtx' in mode):             leg = ROOT.TLegend(0.2, 0.2, 0.6, 0.6)

  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(0, "Signal:", "")
  for wp in workingPoints[tag]: leg.AddEntry(sigEff[wp], wp.name.replace('106X','run3'), "pl")

  if backgroundTree:
    leg.AddEntry(0, "", "")
    leg.AddEntry(0, "Background:", "")
    for wp in workingPoints[tag]: leg.AddEntry(bgEff[wp], wp.name.replace('106X','run3'), "pl")
  leg.Draw("same");

  lat = ROOT.TLatex(0.53, 0.93, comment);
  lat.SetNDC(True);
  lat.Draw("same");

  cmsText = ROOT.TLatex()
  cmsText.SetNDC();
  cmsText.SetTextAlign(11);
  cmsText.SetTextFont(61);
  cmsText.SetTextSize(0.04);
  cmsText.DrawLatex(0.12, 0.93, "CMS");

  extraText = ROOT.TLatex()
  extraText.SetNDC();
  extraText.SetTextAlign(11);
  extraText.SetTextFont(52);
  extraText.SetTextSize(0.04);
  extraText.DrawLatex(0.21, 0.93, "Preliminary");

  dirName  = os.path.join('figures', 'efficiencies', tag, selectVar)
  fileName = os.path.join(dirName, "eff_" + ((region + '_') if mode != 'eta' else '') + mode + '.png')
  c1.Print(makeSubDirs(fileName))
  fileName = os.path.join(dirName, "eff_" + ((region + '_') if mode != 'eta' else '') + mode + '.pdf')
  c1.Print(makeSubDirs(fileName))


import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--tag',       action='store',      default='prelim2017')
argParser.add_argument('--mode',      action='store',      default=None)
argParser.add_argument('--region',    action='store',      default=None)
argParser.add_argument('--selectVar', action='store',      default='all')
argParser.add_argument('--runLocal',  action='store_true')
argParser.add_argument('--subJob',    action='store_true')
args = argParser.parse_args()


#
# Run with above arguments, if None for an argument submit jobs on qsub using all its variations
#
if args.mode and args.region and args.selectVar != 'all':
  drawEfficiency(args.mode, args.tag, args.region, args.selectVar)
elif not args.subJob:
  for mode in ([args.mode] if args.mode else ['eta','pt', 'genPt', 'pt_2TeV', 'genPt_2TeV', 'nvtx', 'pt_1to500GeV','pt_500to1000GeV']):
    for region in ([args.region] if args.region else (['barrel','endcap','extend'] if mode != 'eta' else ['full'])):
      for selectVar in [args.selectVar] if args.selectVar!='all' else ["", "expectedMissingInnerHits", "full5x5_sigmaIetaIeta", "dEtaSeed", "dPhiIn", "hOverE", "relIsoWithEA", "ooEmooP"]:
        if args.runLocal: drawEfficiency(mode, args.tag, region, selectVar)
        else:
          command = './drawEfficiency.py --subJob --mode=%s --tag=%s --region=%s --selectVar=%s' % (mode, args.tag, region, selectVar)
          logFile = 'log/%s-%s-%s-%s.log' % (mode, args.tag, region, selectVar)
          os.system('mkdir -p log')
          os.system("qsub -v dir=" + os.getcwd() + ",command=\"" + command + "\" -q localgrid@cream02 -o " + logFile + " -e " + logFile + " -l walltime=2:00:00 $CMSSW_BASE/src/ElectronID/runOnCream02.sh")
