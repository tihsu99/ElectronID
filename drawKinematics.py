#! /usr/bin/env python

import ROOT
from common import loadClasses, getTreeFromFile
loadClasses('OptimizationConstants.hh')

dateTag = "2019-08-23"


def drawKinematics(var):
  signalTree     = getTreeFromFile(dateTag + '/' + "DY_ext_flat_ntuple_true_alleta_full.root", ROOT.Opt.signalTreeName)
  backgroundTree = getTreeFromFile(dateTag + '/' + "TT_flat_ntuple_trueAndFake_alleta_full.root", ROOT.Opt.backgroundTreeName)

  c1 = ROOT.TCanvas("c1","c1",10,10,600,600)
  ROOT.gStyle.SetOptStat(0)
  c1.cd()

  def drawFromTree(tree, selectionCuts, varName, histName, binning, kinWeights):
    hist = ROOT.TH1F(histName, '', *binning)
    hist.Sumw2()
    cutString  = ("kinWeight*(%s)" if kinWeights else "(%s)") % selectionCuts.GetTitle()
    command    = varName + '>>' + histName
    tree.Draw(command, cutString, "goff")
    return hist

  if var=='eta': varName = 'etaSC'
  else:          varName = var

  if   var=='pt':  binning = (100, 0., 200.)
  elif var=='eta': binning = (50, -2.5, 2.5)

  preselectionCuts = ROOT.Opt.ptCut + ROOT.Opt.otherPreselectionCuts
  sig  = drawFromTree(signalTree,     (preselectionCuts + ROOT.Opt.trueEleCut), varName, 'sig',  binning, False)
  sigW = drawFromTree(signalTree,     (preselectionCuts + ROOT.Opt.trueEleCut), varName, 'sigW', binning, True)
  bg   = drawFromTree(backgroundTree, (preselectionCuts + ROOT.Opt.fakeEleCut), varName, 'bg',   binning, False)

  for h in [sig, sigW, bg]: h.Scale(1./h.GetSumOfWeights())
  for h in [sig, sigW, bg]: h.SetLineWidth(2)

  bg.GetXaxis().SetTitle('p_{T} [GeV]' if var=='pt' else '#eta_{SC}')
  bg.GetYaxis().SetTitle('1/N dN/dx')
  bg.GetYaxis().SetTitleOffset(1.5)
  bg.Draw('hist')

  sig.SetLineColor(ROOT.kRed)
  sig.Draw('hist same')

  sigW.SetLineColor(ROOT.kBlue)
  sigW.SetMarkerColor(ROOT.kBlue)
  sigW.SetMarkerStyle(20)
  sigW.SetMarkerSize(1)
  sigW.Draw('same,pe')

  if var=='pt':    leg = ROOT.TLegend(0.5, 0.5, 0.8, 0.8)
  elif var=='eta': leg = ROOT.TLegend(0.3, 0.15, 0.7, 0.5)
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.AddEntry(sig,  'signal',            'l')
  leg.AddEntry(bg,   'background',        'l')
  leg.AddEntry(sigW, 'signal reweighted', 'pl')
  leg.Draw('same')

  c1.Print('figures/plot_kinematics_' + var + '.png')

drawKinematics('pt')
drawKinematics('eta')
