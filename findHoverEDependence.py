#! /usr/bin/env python

import ROOT, numpy
from common import getTreeFromFile, makeSubDirs
from math import sqrt

dateTag = "2018-03-18"

def mergeBinning(hist):
  varBinLimits = [0]
  varBinIndex  = [1]
  startBinOrig = 0

  def addBins(bins, step):
    for i in range(0, bins):
      varBinIndex.append(varBinIndex[-1] + step)
      varBinLimits.append(hist.GetXaxis().GetBinUpEdge(varBinIndex[-1]))

  addBins(50, 2)
  addBins(20, 5)
  addBins(20, 10)
  addBins(4,  50)
  addBins(4,  100)

  varBinWeightedCenters = []
  varBinRangePos = []
  varBinRangeNeg = []

  # Find weighted centers
  for i in range(len(varBinLimits)-1):
    histTmp = hist.ProjectionX().Clone("histTmp")
    for j in range(histTmp.GetNbinsX()+2):
      if j < varBinIndex[i] or j > varBinIndex[i+1]:
        histTmp.SetBinContent(j, 0)
        histTmp.SetBinError(j, 0)

    mean = histTmp.GetMean()
    varBinWeightedCenters.append(mean)
    varBinRangeNeg.append(mean - varBinLimits[i])
    varBinRangePos.append(varBinLimits[i+1] - mean)

  return varBinIndex, varBinWeightedCenters, varBinRangeNeg, varBinRangePos


def interpolate(x1, x2, y1, y2, y):
  return x1 + (y-y1)*(x2-x1)/(y2-y1)

def findCutoff(hist, cutOffFraction):
  if hist.GetEntries() < 10: return (0,0,0)

  if not cutOffFraction: # use mean based 
    return (hist.GetMean(), hist.GetRMS()/hist.GetSumOfWeights(), hist.GetRMS()/hist.GetSumOfWeights())

  total = hist.Integral()
  grEff = ROOT.TGraphErrors(0);
  for i in range(1, hist.GetNbinsX()+1):
    eff           = sum(hist.GetBinContent(j) for j in range(1, i+1))/total if total > 0 else 0
    effErr        = sqrt(eff*(1-eff)/total)                                 if total > 0 else 1
    cutVal        = hist.GetXaxis().GetBinUpEdge(i)
    cutValErr     = hist.GetXaxis().GetBinWidth(i)/2.0
    grEff.SetPoint(i, cutVal, eff);
    grEff.SetPointError(i, cutValErr, effErr);

  xArray       = grEff.GetX()
  effArray     = grEff.GetY()
  effArrayUp   = [effArray[i]+grEff.GetEY()[i] for i in range(len(effArray))]
  effArrayDown = [effArray[i]-grEff.GetEY()[i] for i in range(len(effArray))]

  def cutOff(array):
    for i in range(1, len(array)):
      if(array[i] > cutOffFraction):
        return interpolate(xArray[i-1], xArray[i], array[i-1], array[i], cutOffFraction)

  xcutoff     = cutOff(effArray)
  xcutoffUp   = cutOff(effArrayUp)-xcutoff
  xcutoffDown = xcutoff-cutOff(effArrayDown)

  return xcutoff, xcutoffUp, xcutoffDown

def findDependence(region, option): 
  signalTree = getTreeFromFile('~/eleIdTuning/tuples/DYJetsToLL_cutID_tuning_94X_v3.root', 'ntupler/ElectronTree')

  cut = 'isTrue==1 && pt>10 && hOverE>0'
  if region=='barrel': cut += '&& abs(etaSC)<1.4442'
  else:                cut += '&& abs(etaSC)>1.566 && abs(etaSC)<2.5'

  if option=='hoeVsE':
    xAxisTitle     = 'E_{SC}'
    yAxisTitle     = 'H/E'
    binning        = (1000, 0, 1000, 1000, 0, .5)
    specialBinning = True
    command        = 'hOverE:eSC'
    fitFunc        = '[0]+[1]/x'
    funcText       = 'f(E)=A+B/E, B=%.2e #pm %.2e' 
    fitMin, fitMax = 20, 1000
    cutOffFraction = None
  elif option=='hVsRho':
    yAxisTitle     = 'HCAL energy [GeV]'
    xAxisTitle     = '#rho'
    binning        = (51, -0.5, 50.5, 1000, 0, 1000)
    specialBinning = False
    command        = 'hOverE*eSC:rho'
    fitFunc        = '[0]+[1]*x'
    funcText       = 'slope %.2e #pm %.2e'
    fitMin, fitMax = 5, 30
    cutOffFraction = 0.90
  elif option=='hoeCorrVsE':
    Crho           = (0.0368 if region == 'barrel' else 0.201)
    xAxisTitle     = 'E_{SC}'
    yAxisTitle     = 'H/E-%.4f#rho/E' % Crho 
    binning        = (1000, 0, 1000, 1000, 0, .5)
    specialBinning = True
    command        = 'hOverE-%.4f*rho/eSC:eSC' % Crho
    fitFunc        = '[0]+[1]/x'
    funcText       = 'f(E)=A+B/E, B=%.2e #pm %.2e' 
    fitMin, fitMax = 20, 1000
    cutOffFraction = None

  hist2D = ROOT.TH2F(option, '', *binning)
  signalTree.Draw(command + '>>' + option, cut, 'colz')

  if specialBinning:
    varBinIndex, varBinWeightedCenters, varBinRangeNeg, varBinRangePos = mergeBinning(hist2D)

  # Loop over all slices along X and determine the cut-off
  graph = ROOT.TGraphAsymmErrors()
  yhigh = 0
  for i in range(len(varBinWeightedCenters) if specialBinning else hist2D.GetNbinsX()):
    hist1D = hist2D.ProjectionY("_py%d" % i, varBinIndex[i] if specialBinning else i+1, varBinIndex[i+1] if specialBinning else i+1)
    val, errPos, errNeg = findCutoff(hist1D, cutOffFraction)
    xCurrent            = varBinWeightedCenters[i] if specialBinning else hist2D.GetXaxis().GetBinCenter(i+1)
    xRangeLow           = varBinRangeNeg[i]        if specialBinning else hist2D.GetXaxis().GetBinWidth(i+1)/2.
    xRangeHigh          = varBinRangePos[i]        if specialBinning else hist2D.GetXaxis().GetBinWidth(i+1)/2.
    graph.SetPoint(i, xCurrent, val)
    graph.SetPointError(i, xRangeLow, xRangeHigh, errPos, errNeg)
    if val > yhigh: yhigh = val

  func = ROOT.TF1("func", fitFunc, fitMin, fitMax)
  func.SetParLimits(0,0,100)
  func.SetParLimits(1,0,10)
  graph.Fit("func","WQRN");
  graph.Fit("func","RM+");

  c1 = ROOT.TCanvas("c1", "c1", 10, 10, 800, 800)
  c1.cd()
  ROOT.gStyle.SetOptStat(0)

  dummy = ROOT.TH2F("dummy","",100, 0, hist2D.GetXaxis().GetBinUpEdge(hist2D.GetNbinsX()), 100, 0, yhigh*1.5)
  dummy.GetXaxis().SetTitle(xAxisTitle)
  dummy.GetYaxis().SetTitle(yAxisTitle)
  dummy.GetYaxis().SetTitleOffset(1.4)
  dummy.Draw()

  graph.SetMarkerStyle(20)
  graph.SetMarkerSize(1.0)
  graph.Draw("P,same")

  lat1 = ROOT.TLatex(0.15, 0.8, region)
  lat2 = ROOT.TLatex(0.15, 0.75, "contours at %.0f" % (100*cutOffFraction) if cutOffFraction else "markers: mean of H/E in E_{SC} slices")
  lat3 = ROOT.TLatex(0.15, 0.7,  funcText % (func.GetParameter(1), func.GetParError(1)))

  for l in [lat2,lat3]:        l.SetTextSize(0.03)
  for l in [lat1, lat2, lat3]: l.SetNDC(True)
  for l in [lat1, lat2, lat3]: l.Draw()

  c1.Print(makeSubDirs('figures/hOverEdependence/' + option + '_' + region + '.png'))

for region in ['barrel', 'endcap']:
  findDependence(region, 'hVsRho')
  findDependence(region, 'hoeVsE')
  findDependence(region, 'hoeCorrVsE')
