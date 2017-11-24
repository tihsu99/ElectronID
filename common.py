import ROOT,os,shutil,numpy
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
  'default':     [workingPoint('Veto',                'cuts_barrel_2017-11-07_WP_Veto',           'cuts_endcap_2017-11-07_WP_Veto',           2, 9,  2, 3),
                  workingPoint('Loose',               'cuts_barrel_2017-11-07_WP_Loose',          'cuts_endcap_2017-11-07_WP_Loose',          4, 46, 1, 1),
                  workingPoint('Medium',              'cuts_barrel_2017-11-07_WP_Medium',         'cuts_endcap_2017-11-07_WP_Medium',         6, 30, 1, 1),
                  workingPoint('Tight',               'cuts_barrel_2017-11-07_WP_Tight',          'cuts_endcap_2017-11-07_WP_Tight',          8, 42, 1, 1)],
  'retuneMVA':   [workingPoint('Veto (retune MVA)',   'cuts_barrel_2017-11-16_WP_Veto',           'cuts_endcap_2017-11-16_WP_Veto',           2, 9,  2, 3),
                  workingPoint('Loose (retune MVA)',  'cuts_barrel_2017-11-16_WP_Loose',          'cuts_endcap_2017-11-16_WP_Loose',          4, 46, 1, 1),
                  workingPoint('Medium (retune MVA)', 'cuts_barrel_2017-11-16_WP_Medium',         'cuts_endcap_2017-11-16_WP_Medium',         6, 30, 1, 1),
                  workingPoint('Tight (retune MVA)',  'cuts_barrel_2017-11-16_WP_Tight',          'cuts_endcap_2017-11-16_WP_Tight',          8, 42, 1, 1)],
  'retuneEff':   [workingPoint('Veto (tune C0)',      'cuts_barrel_2017-11-07_retuned_WP_Veto',   'cuts_endcap_2017-11-07_retuned_WP_Veto',   2, 9,  2, 3),
                  workingPoint('Loose (tune C0)',     'cuts_barrel_2017-11-07_retuned_WP_Loose',  'cuts_endcap_2017-11-07_retuned_WP_Loose',  4, 46, 1, 1),
                  workingPoint('Medium (tune C0)',    'cuts_barrel_2017-11-07_retuned_WP_Medium', 'cuts_endcap_2017-11-07_retuned_WP_Medium', 6, 30, 1, 1),
                  workingPoint('Tight (tune C0)',     'cuts_barrel_2017-11-07_retuned_WP_Tight',  'cuts_endcap_2017-11-07_retuned_WP_Tight',  8, 42, 1, 1)],
  '2016':        [workingPoint('Veto (2016)',         'cuts_barrel_2016_WP_Tight',                'cuts_endcap_2016_WP_Tight',                2, 9,  2, 3),
                  workingPoint('Loose (2016)',        'cuts_barrel_2016_WP_Loose',                'cuts_endcap_2016_WP_Loose',                4, 46, 1, 1),
                  workingPoint('Medium (2016)',       'cuts_barrel_2016_WP_Loose',                'cuts_endcap_2016_WP_Loose',                6, 30, 1, 1),
                  workingPoint('Tight (2016)',        'cuts_barrel_2016_WP_Loose',                'cuts_endcap_2016_WP_Loose',                8, 42, 1, 1)],
}

from itertools import chain
def compareWP(set1, set2, selectPoints):
  zippedList = list(chain(*zip(workingPoints[set1], workingPoints[set2])))
  return [wp for wp in zippedList if wp.name.split()[0] in selectPoints]

workingPoints['compare2016']      = compareWP('default',  '2016',      ['Loose','Tight'])
workingPoints['compareRetuneMVA'] = compareWP('default',  'retuneMVA', ['Loose','Tight'])
workingPoints['compareRetuneEff'] = compareWP('retuneEff','retuneMVA', ['Loose','Tight'])

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

def makeSubDirs(fileName):
  listOfSubDirs = fileName.split('/')[:-1]
  dirName = os.path.join(*listOfSubDirs)
  try:    os.makedirs(dirName)
  except: pass
  for i in range(2,len(listOfSubDirs)+1): shutil.copy("figures/index.php", os.path.join(*listOfSubDirs[:i]))
  return fileName
