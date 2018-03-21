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
class workingPoint:
  def __init__(self, name, cutsFileBarrel, cutsFileEndcap, missingHitsBarrel, missingHitsEndcap):
    self.name              = name
    self.cutsFileBarrel    = cutsFileBarrel
    self.cutsFileEndcap    = cutsFileEndcap
    self.missingHitsBarrel = missingHitsBarrel
    self.missingHitsEndcap = missingHitsEndcap

workingPoints = {
  '2016':        [workingPoint('Veto (2016)',              'cuts_barrel_2016_WP_Tight',                 'cuts_endcap_2016_WP_Tight',                 2, 3),
                  workingPoint('Loose (2016)',             'cuts_barrel_2016_WP_Loose',                 'cuts_endcap_2016_WP_Loose',                 1, 1),
                  workingPoint('Medium (2016)',            'cuts_barrel_2016_WP_Loose',                 'cuts_endcap_2016_WP_Loose',                 1, 1),
                  workingPoint('Tight (2016)',             'cuts_barrel_2016_WP_Loose',                 'cuts_endcap_2016_WP_Loose',                 1, 1)],

  'default':     [workingPoint('Veto',                     'cuts_barrel_2017-11-07_WP_Veto',            'cuts_endcap_2017-11-07_WP_Veto',            2, 3),
                  workingPoint('Loose',                    'cuts_barrel_2017-11-07_WP_Loose',           'cuts_endcap_2017-11-07_WP_Loose',           1, 1),
                  workingPoint('Medium',                   'cuts_barrel_2017-11-07_WP_Medium',          'cuts_endcap_2017-11-07_WP_Medium',          1, 1),
                  workingPoint('Tight',                    'cuts_barrel_2017-11-07_WP_Tight',           'cuts_endcap_2017-11-07_WP_Tight',           1, 1)],
  'retuneMVA':   [workingPoint('Veto (retune MVA)',        'cuts_barrel_2017-11-16_WP_Veto',            'cuts_endcap_2017-11-16_WP_Veto',            2, 3),
                  workingPoint('Loose (retune MVA)',       'cuts_barrel_2017-11-16_WP_Loose',           'cuts_endcap_2017-11-16_WP_Loose',           1, 1),
                  workingPoint('Medium (retune MVA)',      'cuts_barrel_2017-11-16_WP_Medium',          'cuts_endcap_2017-11-16_WP_Medium',          1, 1),
                  workingPoint('Tight (retune MVA)',       'cuts_barrel_2017-11-16_WP_Tight',           'cuts_endcap_2017-11-16_WP_Tight',           1, 1)],
  'retuneEff':   [workingPoint('Veto (tune C0)',           'cuts_barrel_2017-11-07_retuned_WP_Veto',    'cuts_endcap_2017-11-07_retuned_WP_Veto',    2, 3),
                  workingPoint('Loose (tune C0)',          'cuts_barrel_2017-11-07_retuned_WP_Loose',   'cuts_endcap_2017-11-07_retuned_WP_Loose',   1, 1),
                  workingPoint('Medium (tune C0)',         'cuts_barrel_2017-11-07_retuned_WP_Medium',  'cuts_endcap_2017-11-07_retuned_WP_Medium',  1, 1),
                  workingPoint('Tight (tune C0)',          'cuts_barrel_2017-11-07_retuned_WP_Tight',   'cuts_endcap_2017-11-07_retuned_WP_Tight',   1, 1)],
  'retuneEff2':  [workingPoint('Veto (tune C0, 2)',        'cuts_barrel_2017-11-07_retuned2_WP_Veto',   'cuts_endcap_2017-11-07_retuned2_WP_Veto',   2, 3),
                  workingPoint('Loose (tune C0, 2)',       'cuts_barrel_2017-11-07_retuned2_WP_Loose',  'cuts_endcap_2017-11-07_retuned2_WP_Loose',  1, 1),
                  workingPoint('Medium (tune C0, 2)',      'cuts_barrel_2017-11-07_retuned2_WP_Medium', 'cuts_endcap_2017-11-07_retuned2_WP_Medium', 1, 1),
                  workingPoint('Tight (tune C0, 2)',       'cuts_barrel_2017-11-07_retuned2_WP_Tight',  'cuts_endcap_2017-11-07_retuned2_WP_Tight',  1, 1)],
  'retuneEff3':  [workingPoint('Veto (tune C0, 3)',        'cuts_barrel_2017-11-07_retuned3_WP_Veto',   'cuts_endcap_2017-11-07_retuned3_WP_Veto',   2, 3),
                  workingPoint('Loose (tune C0, 3)',       'cuts_barrel_2017-11-07_retuned3_WP_Loose',  'cuts_endcap_2017-11-07_retuned3_WP_Loose',  1, 1),
                  workingPoint('Medium (tune C0, 3)',      'cuts_barrel_2017-11-07_retuned3_WP_Medium', 'cuts_endcap_2017-11-07_retuned3_WP_Medium', 1, 1),
                  workingPoint('Tight (tune C0, 3)',       'cuts_barrel_2017-11-07_retuned3_WP_Tight',  'cuts_endcap_2017-11-07_retuned3_WP_Tight',  1, 1)],
  'retuneEff4':  [workingPoint('Veto (tune C0, CE=0.5)',   'cuts_barrel_2017-11-07_retuned4_WP_Veto',   'cuts_endcap_2017-11-07_retuned4_WP_Veto',   2, 3),
                  workingPoint('Loose (tune C0, CE=0.5)',  'cuts_barrel_2017-11-07_retuned4_WP_Loose',  'cuts_endcap_2017-11-07_retuned4_WP_Loose',  1, 1),
                  workingPoint('Medium (tune C0, CE=0.5)', 'cuts_barrel_2017-11-07_retuned4_WP_Medium', 'cuts_endcap_2017-11-07_retuned4_WP_Medium', 1, 1),
                  workingPoint('Tight (tune C0, CE=0.5)',  'cuts_barrel_2017-11-07_retuned4_WP_Tight',  'cuts_endcap_2017-11-07_retuned4_WP_Tight',  1, 1)],
  'prelim2017':  [workingPoint('Veto',                     'cuts_barrel_2017-11-07_retuned5_WP_Veto',   'cuts_endcap_2017-11-07_retuned5_WP_Veto',   2, 3),
                  workingPoint('Loose',                    'cuts_barrel_2017-11-07_retuned5_WP_Loose',  'cuts_endcap_2017-11-07_retuned5_WP_Loose',  1, 1),
                  workingPoint('Medium',                   'cuts_barrel_2017-11-07_retuned5_WP_Medium', 'cuts_endcap_2017-11-07_retuned5_WP_Medium', 1, 1),
                  workingPoint('Tight',                    'cuts_barrel_2017-11-07_retuned5_WP_Tight',  'cuts_endcap_2017-11-07_retuned5_WP_Tight',  1, 1)],

  'training94':  [workingPoint('Veto',                     'cuts_barrel_2018-03-18_WP_Veto',            'cuts_endcap_2018-03-18_WP_Veto',            2, 3),
                  workingPoint('Loose',                    'cuts_barrel_2018-03-18_WP_Loose',           'cuts_endcap_2018-03-18_WP_Loose',           1, 1),
                  workingPoint('Medium',                   'cuts_barrel_2018-03-18_WP_Medium',          'cuts_endcap_2018-03-18_WP_Medium',          1, 1),
                  workingPoint('Tight',                    'cuts_barrel_2018-03-18_WP_Tight',           'cuts_endcap_2018-03-18_WP_Tight',           1, 1)],
}

from itertools import chain
def compareWP(set1, set2, selectPoints):
  zippedList = list(chain(*zip(workingPoints[set1], workingPoints[set2])))
  return [wp for wp in zippedList if wp.name.split()[0] in selectPoints]

workingPoints['compare2016']      = compareWP('default',  '2016',       ['Loose','Tight'])
workingPoints['compareRetuneMVA'] = compareWP('default',  'retuneMVA',  ['Loose','Tight'])
workingPoints['compareRetuneEff'] = compareWP('default',  'retuneEff',  ['Loose','Tight'])
workingPoints['compareRetune2']   = compareWP('retuneEff','retuneEff2', ['Medium','Tight'])
workingPoints['compareRetune4']   = compareWP('retuneEff','retuneEff4', ['Medium','Tight'])
workingPoints['compareRetune4a']  = compareWP('retuneEff','retuneEff4', ['Veto','Loose'])

def setColors(set):
  sColors  = [2, 4, 6, 8]
  bgColors = [9, 46, 30, 42]
  for i, wp in enumerate(set):
    wp.sColor  = sColors[i]
    wp.bgColor = bgColors[i]

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
