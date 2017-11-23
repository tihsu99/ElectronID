import ROOT,os
ROOT.gROOT.SetBatch(True)

#
# Load custom classes
#
def loadClasses(*args):
  for c in args:
    os.system('root -b -q ' + c + '+ &> /dev/null')
    ROOT.gSystem.Load(c.replace('.','_'))

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
# Helper functions
#
def getTreeFromFile(fname, tname):
  chain = ROOT.TChain(tname)
  chain.Add(fname)
  return chain

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
