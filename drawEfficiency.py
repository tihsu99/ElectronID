#! /usr/bin/env python

import ROOT,os,numpy,shutil
os.system('root -b -q VarCut.cc+ &> /dev/null')
os.system('root -b -q OptimizationConstants.hh+ &> /dev/null')
ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load('VarCut_cc')
ROOT.gSystem.Load('OptimizationConstants_hh')

dateTag = "2017-11-16"

#
# Listing workingpoints to compare
#
from collections import namedtuple
workingPoint = namedtuple('workingpoint', 'name cutsFileBarrel cutsFileEndcap sColor bgColor missingHitsBarrel missingHitsEndcap')

workingPoints = {
  'default':     [workingPoint('Veto',          'cuts_barrel_2017-11-07_WP_Veto',           'cuts_endcap_2017-11-07_WP_Veto',           2, 9,  2, 3),
                  workingPoint('Loose',         'cuts_barrel_2017-11-07_WP_Loose',          'cuts_endcap_2017-11-07_WP_Loose',          4, 46, 1, 1),
                  workingPoint('Medium',        'cuts_barrel_2017-11-07_WP_Medium',         'cuts_endcap_2017-11-07_WP_Medium',         6, 30, 1, 1),
                  workingPoint('Tight',         'cuts_barrel_2017-11-07_WP_Tight',          'cuts_endcap_2017-11-07_WP_Tight',          8, 42, 1, 1)],
  'retuneMVA':   [workingPoint('Veto',          'cuts_barrel_2017-11-16_WP_Veto',           'cuts_endcap_2017-11-16_WP_Veto',           2, 9,  2, 3),
                  workingPoint('Loose',         'cuts_barrel_2017-11-16_WP_Loose',          'cuts_endcap_2017-11-16_WP_Loose',          4, 46, 1, 1),
                  workingPoint('Medium',        'cuts_barrel_2017-11-16_WP_Medium',         'cuts_endcap_2017-11-16_WP_Medium',         6, 30, 1, 1),
                  workingPoint('Tight',         'cuts_barrel_2017-11-16_WP_Tight',          'cuts_endcap_2017-11-16_WP_Tight',          8, 42, 1, 1)],
  'retuneEff':   [workingPoint('Veto',          'cuts_barrel_2017-11-07_retuned_WP_Veto',   'cuts_endcap_2017-11-07_retuned_WP_Veto',   2, 9,  2, 3),
                  workingPoint('Loose',         'cuts_barrel_2017-11-07_retuned_WP_Loose',  'cuts_endcap_2017-11-07_retuned_WP_Loose',  4, 46, 1, 1),
                  workingPoint('Medium',        'cuts_barrel_2017-11-07_retuned_WP_Medium', 'cuts_endcap_2017-11-07_retuned_WP_Medium', 6, 30, 1, 1),
                  workingPoint('Tight',         'cuts_barrel_2017-11-07_retuned_WP_Tight',  'cuts_endcap_2017-11-07_retuned_WP_Tight',  8, 42, 1, 1)],
  'compare2016': [workingPoint('Tight (2016)',  'cuts_barrel_2016_WP_Tight',                'cuts_endcap_2016_WP_Tight',                6, 30, 1, 1),
                  workingPoint('Tight',         'cuts_barrel_2017-11-07_WP_Tight',          'cuts_endcap_2017-11-07_WP_Tight',          8, 42, 1, 1),
                  workingPoint('Loose (2016)',  'cuts_barrel_2016_WP_Loose',                'cuts_endcap_2016_WP_Loose',                2, 9,  1, 1),
                  workingPoint('Loose',         'cuts_barrel_2017-11-07_WP_Loose',          'cuts_endcap_2017-11-07_WP_Loose',          4, 46, 1, 1)],
}


#
# Helper functions for trees and binning
#

def getTreeFromFile(fname, tname):
  chain = ROOT.TChain(tname)
  chain.Add(fname)
  return chain

def getPtBins(longPtRange):
  bins = [20.]
  for i in range(20): bins.append(bins[-1] + 1)
  for i in range(30): bins.append(bins[-1] + 2)
  for i in range(10): bins.append(bins[-1] + 5)
  for i in range(5):  bins.append(bins[-1] + 10)
  if(longPtRange):
    for i in range(20):  bins.append(bins[-1] + 20)
    for i in range(10):  bins.append(bins[-1] + 40)
    for i in range(10):  bins.append(bins[-1] + 100)
  return bins;

def getEtaBins():
  nBins = 50
  xmin  = -2.5
  xmax  = 2.5
  delta = (xmax - xmin)/nBins
  
  bins = []
  for i in range(0, nBins+1): bins.append(xmin + i*delta);
  return bins

def getNvtxBins():
  return [0., 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 27, 30, 35, 40, 45, 50]

def drawFromTree(tree, selectionCuts, varName, histName, binning, wp):
  hist = ROOT.TH1F(histName, histName, len(binning)-1, numpy.array(binning))
  hist.Sumw2()
  cutString  = "genWeight*kinWeight*(%s)" % selectionCuts.GetTitle()
  command    = varName + '>>' + histName
  if '2017-11-16' in wp.cutsFileBarrel or 'retuned' in wp.cutsFileBarrel: 
    cutString = cutString.replace('hOverE','hOverEscaled')
    command   = command.replace('hOverE','hOverEscaled')
  tree.Draw(command, cutString, "goff")
  return hist

def getCuts(wp, barrel, selectVar):
  file           = ROOT.TFile('cut_repository/' + (wp.cutsFileBarrel if barrel else wp.cutsFileEndcap) + '.root')
  cuts           = file.Get('cuts')
  selectionCuts  = ROOT.TCut(cuts.getCut(selectVar))
  selectionCuts += ROOT.TCut('expectedMissingInnerHits<='+str(wp.missingHitsBarrel if barrel else wp.missingHitsEndcap))
  return selectionCuts

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
  if '2TeV' in mode: signalFileName = dateTag + '/' + "DoubleEleFlat_flat_ntuple_trueAndFake_alleta_full.root"
  else:              signalFileName = dateTag + '/' + "DYJetsToLL_flat_ntuple_true_alleta_full.root" 
  signalTree     = getTreeFromFile(signalFileName,                                                    ROOT.Opt.signalTreeName)
  backgroundTree = getTreeFromFile(dateTag + '/' + "TTJets_flat_ntuple_trueAndFake_alleta_full.root", ROOT.Opt.backgroundTreeName)

  if('pt' in mode):     binning = getPtBins(mode=="pt_2TeV")
  elif('eta' in mode):  binning = getEtaBins()
  elif('nvtx' in mode): binning = getNvtxBins()

  if('pt' in mode):     axisLabel = "p_{T} [GeV]"
  elif('eta' in mode):  axisLabel = "#eta_{SC}"
  elif('nvtx' in mode): axisLabel = "Nvtx"

  if(mode == "eta"): comment = '';
  else:              comment = region + ' electrons'

  preselectionCuts = ROOT.Opt.ptCut + ROOT.Opt.otherPreselectionCuts
  if(mode == "eta"):        preselectionCuts += ROOT.TCut('abs(etaSC)<2.5')
  elif(region == "barrel"): preselectionCuts += ROOT.Opt.etaCutBarrel
  else:                     preselectionCuts += ROOT.Opt.etaCutEndcap

  signalCuts     = preselectionCuts + (ROOT.Opt.trueEleCut if not mode == 'pt_2TeV' else ROOT.TCut())
  backgroundCuts = preselectionCuts + ROOT.Opt.fakeEleCut;

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
  for wp in workingPoints[tag]:
    is2016 = '2016' in wp.name
    wpCutsBarrel = getCuts(wp, True,  selectVar)
    wpCutsEndcap = getCuts(wp, False, selectVar)

    selectionCutsBarrel = wpCutsBarrel + ROOT.Opt.etaCutBarrel
    selectionCutsEndcap = wpCutsEndcap + ROOT.Opt.etaCutEndcap
    selectionCuts       = ROOT.TCut('(' + selectionCutsBarrel.GetTitle() + ')||(' + selectionCutsEndcap.GetTitle() + ')')

    if('pt' in mode):     varName = "pt"
    elif('eta' in mode):  varName = "etaSC"
    elif('nvtx' in mode): varName = "nPV"

    sigNum = drawFromTree(signalTree,     (signalCuts + selectionCuts),     varName, 'sigNum'+ wp.name.split()[0], binning, wp)
    sigDen = drawFromTree(signalTree,     (signalCuts),                     varName, 'sigDen'+ wp.name.split()[0], binning, wp)
    bgNum  = drawFromTree(backgroundTree, (backgroundCuts + selectionCuts), varName, 'bgNum' + wp.name.split()[0], binning, wp)
    bgDen  = drawFromTree(backgroundTree, (backgroundCuts),                 varName, 'bgDen' + wp.name.split()[0], binning, wp)

    sigEff[wp] = calculateEffAndErrors(sigNum, sigDen, 'sigEff' + wp.name, binning)
    bgEff[wp]  = calculateEffAndErrors(bgNum,  bgDen,  'bgEff'  + wp.name, binning)

    setHistogram(sigEff[wp], wp, True)
    setHistogram( bgEff[wp], wp, False)

    sigEff[wp].Draw("same,pe")
    if('2TeV' not in mode): bgEff[wp].Draw("same,pe")
    c1.Update()
  
  if('pt' in mode):     leg = ROOT.TLegend(0.2, 0.3, 0.6, 0.7)
  elif('eta' in mode):  leg = ROOT.TLegend(0.38, 0.2, 0.75, 0.6)
  elif('nvtx' in mode): leg = ROOT.TLegend(0.2, 0.2, 0.6, 0.6)

  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(0, "Signal:", "")
  for wp in workingPoints[tag]: leg.AddEntry(sigEff[wp], wp.name, "pl")

  if('2TeV' not in mode):
    leg.AddEntry(0, "", "")
    leg.AddEntry(0, "Background:", "")
    for wp in workingPoints[tag]: leg.AddEntry(bgEff[wp], wp.name, "pl")
  leg.Draw("same");

  lat = ROOT.TLatex(0.5, 0.95, comment);
  lat.SetNDC(True);
  lat.Draw("same");

  def createSubDirs(listOfSubDirs):
    dirName = os.path.join(*listOfSubDirs)
    try:    os.makedirs(dirName)
    except: pass
    for i in range(2,len(listOfSubDirs)+1): shutil.copy("figures/index.php", os.path.join(*listOfSubDirs[:i]))
    return dirName

  subDirs = ['figures','efficiencies']
  if(selectVar != ""): subDirs.append(selectVar)
  if(tag!='default'):  subDirs.append(tag)
  dirName = createSubDirs(subDirs)

  fileName = os.path.join(dirName, "eff_" + ((region + '_') if mode != 'eta' else '') + mode + '.png')
  if '2017-11-16' in wp.cutsFileBarrel or 'retuned' in wp.cutsFileBarrel: fileName.replace('hOverE','hOverEscaled')
  c1.Print(fileName) 


import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--tag',       action='store',      default='default')
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
  for mode in ([args.mode] if args.mode else ['eta','pt', 'pt_2TeV', 'nvtx']):
    for region in ([args.region] if args.region else (['endcap','barrel'] if mode != 'eta' else ['full'])):
      for selectVar in [args.selectVar] if args.selectVar!='all' else ["", "expectedMissingInnerHits", "full5x5_sigmaIetaIeta", "dEtaSeed", "dPhiIn", "hOverE", "relIsoWithEA", "ooEmooP"]:
        if args.runLocal: drawEfficiency(mode, args.tag, region, selectVar)
        else:             
          command = './drawEfficiency.py --subJob --mode=%s --tag=%s --region=%s --selectVar=%s' % (mode, args.tag, region, selectVar)
          logFile = 'log/%s-%s-%s-%s.log' % (mode, args.tag, region, selectVar)
          os.system('mkdir -p log')
          os.system("qsub -v dir=" + os.getcwd() + ",command=\"" + command + "\" -q localgrid@cream02 -o " + logFile + " -e " + logFile + " -l walltime=1:00:00 $CMSSW_BASE/src/ElectronID/runOnCream02.sh")


