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
  def __init__(self, name, cutsFile, missingHitsBarrel, missingHitsEndcap):
    self.name              = name
    self.cutsFileBarrel    = cutsFile.replace('region','barrel')
    self.cutsFileEndcap    = cutsFile.replace('region','endcap')
    self.missingHitsBarrel = missingHitsBarrel
    self.missingHitsEndcap = missingHitsEndcap

workingPoints = {
  '2016':        [workingPoint('Veto (2016)',              'cuts_region_2016_WP_Tight',                 2, 3),
                  workingPoint('Loose (2016)',             'cuts_region_2016_WP_Loose',                 1, 1),
                  workingPoint('Medium (2016)',            'cuts_region_2016_WP_Loose',                 1, 1),
                  workingPoint('Tight (2016)',             'cuts_region_2016_WP_Loose',                 1, 1)],

  'training92':  [workingPoint('Veto (92X)',               'cuts_region_2017-11-07_WP_Veto',            2, 3),
                  workingPoint('Loose (92X)',              'cuts_region_2017-11-07_WP_Loose',           1, 1),
                  workingPoint('Medium (92X)',             'cuts_region_2017-11-07_WP_Medium',          1, 1),
                  workingPoint('Tight (92X)',              'cuts_region_2017-11-07_WP_Tight',           1, 1)],
  'retuneMVA':   [workingPoint('Veto (retune MVA)',        'cuts_region_2017-11-16_WP_Veto',            2, 3),
                  workingPoint('Loose (retune MVA)',       'cuts_region_2017-11-16_WP_Loose',           1, 1),
                  workingPoint('Medium (retune MVA)',      'cuts_region_2017-11-16_WP_Medium',          1, 1),
                  workingPoint('Tight (retune MVA)',       'cuts_region_2017-11-16_WP_Tight',           1, 1)],
  'retuneEff':   [workingPoint('Veto (tune C0)',           'cuts_region_2017-11-07_retuned_WP_Veto',    2, 3),
                  workingPoint('Loose (tune C0)',          'cuts_region_2017-11-07_retuned_WP_Loose',   1, 1),
                  workingPoint('Medium (tune C0)',         'cuts_region_2017-11-07_retuned_WP_Medium',  1, 1),
                  workingPoint('Tight (tune C0)',          'cuts_region_2017-11-07_retuned_WP_Tight',   1, 1)],
  'retuneEff2':  [workingPoint('Veto (tune C0, 2)',        'cuts_region_2017-11-07_retuned2_WP_Veto',   2, 3),
                  workingPoint('Loose (tune C0, 2)',       'cuts_region_2017-11-07_retuned2_WP_Loose',  1, 1),
                  workingPoint('Medium (tune C0, 2)',      'cuts_region_2017-11-07_retuned2_WP_Medium', 1, 1),
                  workingPoint('Tight (tune C0, 2)',       'cuts_region_2017-11-07_retuned2_WP_Tight',  1, 1)],
  'retuneEff3':  [workingPoint('Veto (tune C0, 3)',        'cuts_region_2017-11-07_retuned3_WP_Veto',   2, 3),
                  workingPoint('Loose (tune C0, 3)',       'cuts_region_2017-11-07_retuned3_WP_Loose',  1, 1),
                  workingPoint('Medium (tune C0, 3)',      'cuts_region_2017-11-07_retuned3_WP_Medium', 1, 1),
                  workingPoint('Tight (tune C0, 3)',       'cuts_region_2017-11-07_retuned3_WP_Tight',  1, 1)],
  'retuneEff4':  [workingPoint('Veto (tune C0, CE=0.5)',   'cuts_region_2017-11-07_retuned4_WP_Veto',   2, 3),
                  workingPoint('Loose (tune C0, CE=0.5)',  'cuts_region_2017-11-07_retuned4_WP_Loose',  1, 1),
                  workingPoint('Medium (tune C0, CE=0.5)', 'cuts_region_2017-11-07_retuned4_WP_Medium', 1, 1),
                  workingPoint('Tight (tune C0, CE=0.5)',  'cuts_region_2017-11-07_retuned4_WP_Tight',  1, 1)],
  'prelim2017':  [workingPoint('Veto (92X)',               'cuts_region_2017-11-07_retuned5b_WP_Veto',  2, 3),
                  workingPoint('Loose (92X)',              'cuts_region_2017-11-07_retuned5b_WP_Loose', 1, 1),
                  workingPoint('Medium (92X)',             'cuts_region_2017-11-07_retuned5b_WP_Medium',1, 1),
                  workingPoint('Tight (92X)',              'cuts_region_2017-11-07_retuned5b_WP_Tight', 1, 1)],

  'training94':  [workingPoint('Veto (94X)',               'cuts_region_2018-03-18_WP_Veto',            2, 3),
                  workingPoint('Loose (94X)',              'cuts_region_2018-03-18_WP_Loose',           1, 1),
                  workingPoint('Medium (94X)',             'cuts_region_2018-03-18_WP_Medium',          1, 1),
                  workingPoint('Tight (94X)',              'cuts_region_2018-03-18_WP_Tight',           1, 1)],
  'retuned94':   [workingPoint('Veto (94X)',               'cuts_region_2018-03-18_retuned_WP_Veto',    2, 3),
                  workingPoint('Loose (94X)',              'cuts_region_2018-03-18_retuned_WP_Loose',   1, 1),
                  workingPoint('Medium (94X)',             'cuts_region_2018-03-18_retuned_WP_Medium',  1, 1),
                  workingPoint('Tight (94X)',              'cuts_region_2018-03-18_retuned_WP_Tight',   1, 1)],
  'retuned94b':  [workingPoint('Veto (94X)',               'cuts_region_2018-03-18_retuned_WP2_Veto',   2, 3),
                  workingPoint('Loose (94X)',              'cuts_region_2018-03-18_retuned_WP2_Loose',  1, 1),
                  workingPoint('Medium (94X)',             'cuts_region_2018-03-18_retuned_WP2_Medium', 1, 1),
                  workingPoint('Tight (94X)',              'cuts_region_2018-03-18_retuned_WP2_Tight',  1, 1)],
  'retuned94c':  [workingPoint('Veto (94X)',               'cuts_region_2018-03-18_retuned_WP3_Veto',   2, 3),
                  workingPoint('Loose (94X)',              'cuts_region_2018-03-18_retuned_WP3_Loose',  1, 1),
                  workingPoint('Medium (94X)',             'cuts_region_2018-03-18_retuned_WP3_Medium', 1, 1),
                  workingPoint('Tight (94X)',              'cuts_region_2018-03-18_retuned_WP3_Tight',  1, 1)],
  'retuned94d':  [workingPoint('Veto (94X)',               'cuts_region_2018-03-18_retuned_WP4_Veto',   2, 3),
                  workingPoint('Loose (94X)',              'cuts_region_2018-03-18_retuned_WP4_Loose',  1, 1),
                  workingPoint('Medium (94X)',             'cuts_region_2018-03-18_retuned_WP4_Medium', 1, 1),
                  workingPoint('Tight (94X)',              'cuts_region_2018-03-18_retuned_WP4_Tight',  1, 1)],
  'retuned94e':  [workingPoint('Veto (94X)',               'cuts_region_2018-03-18_retuned_WP5_Veto',   2, 3),
                  workingPoint('Loose (94X)',              'cuts_region_2018-03-18_retuned_WP5_Loose',  1, 1),
                  workingPoint('Medium (94X)',             'cuts_region_2018-03-18_retuned_WP5_Medium', 1, 1),
                  workingPoint('Tight (94X)',              'cuts_region_2018-03-18_retuned_WP5_Tight',  1, 1)],

  'training106': [workingPoint('Veto (106X)',               'cuts_region_2019-08-23_WP_Veto',           2, 3),
                  workingPoint('Loose (106X)',              'cuts_region_2019-08-23_WP_Loose',          1, 1),
                  workingPoint('Medium (106X)',             'cuts_region_2019-08-23_WP_Medium',         1, 1),
                  workingPoint('Tight (106X)',              'cuts_region_2019-08-23_WP_Tight',          1, 1)],
}

from itertools import chain
def compareWP(set1, set2, selectPoints):
  zippedList = list(chain(*zip(workingPoints[set1], workingPoints[set2])))
  return [wp for wp in zippedList if wp.name.split()[0] in selectPoints]

workingPoints['compareTrainingVM']  = compareWP('training92','training94', ['Veto','Medium'])
workingPoints['compareTrainingLT']  = compareWP('training92','training94', ['Loose','Tight'])
workingPoints['compare94VM']        = compareWP('prelim2017','retuned94',  ['Veto','Medium'])
workingPoints['compare94LT']        = compareWP('prelim2017','retuned94',  ['Loose','Tight'])
workingPoints['compare94bVM']       = compareWP('prelim2017','retuned94b', ['Veto','Medium'])
workingPoints['compare94bLT']       = compareWP('prelim2017','retuned94b', ['Loose','Tight'])
workingPoints['compare94cVM']       = compareWP('prelim2017','retuned94c', ['Veto','Medium'])
workingPoints['compare94cLT']       = compareWP('prelim2017','retuned94c', ['Loose','Tight'])
workingPoints['compare94dVM']       = compareWP('prelim2017','retuned94d', ['Veto','Medium'])
workingPoints['compare94dLT']       = compareWP('prelim2017','retuned94d', ['Loose','Tight'])
workingPoints['compare94eVM']       = compareWP('prelim2017','retuned94e', ['Veto','Medium'])
workingPoints['compare94eLT']       = compareWP('prelim2017','retuned94e', ['Loose','Tight'])
workingPoints['compare106VM']       = compareWP('training106','training94', ['Veto','Medium'])
workingPoints['compare106LT']       = compareWP('training106','training94', ['Loose','Tight'])

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
