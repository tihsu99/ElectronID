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
    self.cutsFileExtend    = cutsFile.replace('region','extend')
    self.missingHitsBarrel = missingHitsBarrel
    self.missingHitsEndcap = missingHitsEndcap
    self.missingHitsExtend = missingHitsEndcap # Set it to be the same with endcap in current stage

HistoryData_barrel = {
  'Tight (94X)': [0.696,0.9968, 16],
  'Medium (94X)': [0.804, 0.993, 14],
  'Loose (94X)' : [0.892, 0.984, 12],
  'Veto (94X)' : [0.942, 0.964, 11]
}
HistoryData_endcap = {
  'Veto (94X)': [0.941, 0.894, 11],
  'Loose (94X)': [0.89, 0.945, 12],
  'Medium (94X)' : [0.799, 0.974, 14],
  'Tight (94X)': [0.689, 0.987, 16]
}

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
  'training106_tuneC0': [workingPoint('Veto (106X)',        'cuts_region_2019-08-23_training106_tuneC0_Veto',           2, 3),
                  workingPoint('Loose (106X)',              'cuts_region_2019-08-23_training106_tuneC0_Loose',          1, 1),
                  workingPoint('Medium (106X)',             'cuts_region_2019-08-23_training106_tuneC0_Medium',         1, 1),
                  workingPoint('Tight (106X)',              'cuts_region_2019-08-23_training106_tuneC0_Tight',          1, 1)],
  'training106_tuneC0CE_v3': [workingPoint('Veto (106X)',        'cuts_region_2019-08-23_training106_tuneC0CE_v3_Veto',           2, 3),
                  workingPoint('Loose (106X)',              'cuts_region_2019-08-23_training106_tuneC0CE_v3_Loose',          1, 1),
                  workingPoint('Medium (106X)',             'cuts_region_2019-08-23_training106_tuneC0CE_v3_Medium',         1, 1),
                  workingPoint('Tight (106X)',              'cuts_region_2019-08-23_training106_tuneC0CE_v3_Tight',          1, 1)],
 'training106_tuneC0CE_v4': [workingPoint('Veto (106X)',        'cuts_region_2019-08-23_training106_tuneC0CE_v4_Veto',           2, 3),
                  workingPoint('Loose (106X)',              'cuts_region_2019-08-23_training106_tuneC0CE_v4_Loose',          1, 1),
                  workingPoint('Medium (106X)',             'cuts_region_2019-08-23_training106_tuneC0CE_v4_Medium',         1, 1),
                  workingPoint('Tight (106X)',              'cuts_region_2019-08-23_training106_tuneC0CE_v4_Tight',          1, 1)],
 'training106_tuneC0CE_v5': [workingPoint('Veto (106X)',        'cuts_region_2019-08-23_training106_tuneC0CE_v5_Veto',           2, 3),
                  workingPoint('Loose (106X)',              'cuts_region_2019-08-23_training106_tuneC0CE_v5_Loose',          1, 1),
                  workingPoint('Medium (106X)',             'cuts_region_2019-08-23_training106_tuneC0CE_v5_Medium',         1, 1),
                  workingPoint('Tight (106X)',              'cuts_region_2019-08-23_training106_tuneC0CE_v5_Tight',          1, 1)],
   'training106_tuneCE_vs': [workingPoint('Veto (106X) (fix CE)',        'cuts_region_2019-08-23_training106_tuneC0_Veto',           2, 3),
                  workingPoint('Loose (106X) (fix CE)',              'cuts_region_2019-08-23_training106_tuneC0_Loose',          1, 1),
                  workingPoint('Medium (106X) (fix CE)',             'cuts_region_2019-08-23_training106_tuneC0_Medium',         1, 1),
                  workingPoint('Tight (106X) (fix CE)',              'cuts_region_2019-08-23_training106_tuneC0_Tight',          1, 1),
                  workingPoint('Veto (106X) (tune CE)',        'cuts_region_2019-08-23_training106_tuneC0CE_v3_Veto',           2, 3),
                  workingPoint('Loose (106X) (tune CE)',              'cuts_region_2019-08-23_training106_tuneC0CE_v3_Loose',          1, 1),
                  workingPoint('Medium (106X) (tune CE)',             'cuts_region_2019-08-23_training106_tuneC0CE_v3_Medium',         1, 1),
                  workingPoint('Tight (106X) (tune CE)',              'cuts_region_2019-08-23_training106_tuneC0CE_v3_Tight',          1, 1)

  ],
  'training106_tuneC0_norel': [workingPoint('Veto (106X)',               'cuts_region_2019-08-23_training106_tuneC0_norel_Veto',           2, 3),
                  workingPoint('Loose (106X)',              'cuts_region_2019-08-23_training106_tuneC0_norel_Loose',          1, 1),
                  workingPoint('Medium (106X)',             'cuts_region_2019-08-23_training106_tuneC0_norel_Medium',         1, 1),
                  workingPoint('Tight (106X)',              'cuts_region_2019-08-23_training106_tuneC0_norel_Tight',          1, 1)],
  'training94_official': [workingPoint('Veto (94X)',       'cuts_region_2019-08-23_94X_Veto',           2, 3),
                  workingPoint('Loose (94X)',              'cuts_region_2019-08-23_94X_Loose',          1, 1),
                  workingPoint('Medium (94X)',             'cuts_region_2019-08-23_94X_Medium',         1, 1),
                  workingPoint('Tight (94X)',              'cuts_region_2019-08-23_94X_Tight',          1, 1)],
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
workingPoints['compare106v94MT']    = compareWP('training106_tuneC0','training94_official',['Medium','Tight'])

def setColors(set):
  sColors  = [2, 4, 6, 8, 11 ,12, 14 ,16]
  bgColors = [9, 46, 30, 42, 44, 48, 52, 56]
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

def drawFromTree2D(tree, selectionCuts, varNameY, varNameX, histName, binningY, binningX, wp):
  hist = ROOT.TH2F(histName, histName, len(binningX)-1, numpy.array(binningX), len(binningY)-1, numpy.array(binningY))
  cutString  = "genWeight*kinWeight*(%s)" % selectionCuts.GetTitle()
  command    = varNameY + ':' + varNameX + '>>' + histName
  tree.Draw(command, cutString, "goff")
  return hist

def getCuts(wp, region, selectVar):

  file           = ROOT.TFile('cut_repository/' + wp.cutsFileBarrel + '.root')
  if(region == 'barrel'):
    file = ROOT.TFile('cut_repository/' + wp.cutsFileBarrel + '.root')
  elif(region == 'endcap'):
    file = ROOT.TFile('cut_repository/' + wp.cutsFileEndcap + '.root')
  else:
    file = ROOT.TFile('cut_repository/' + wp.cutsFileExtend + '.root')

  cuts           = file.Get('cuts')
  if(region == 'barrel'):
    print('cut_repository/' + wp.cutsFileBarrel + '.root')
  elif(region == 'endcap'):
    print('cut_repository/' + wp.cutsFileEndcap + '.root')
  else:
    print('cut_repository/' + wp.cutsFileExtend + '.root')
  cuts.printCuts()

  selectionCuts  = ROOT.TCut(cuts.getCut(selectVar))

  if(region == 'barrel'):
    selectionCuts += ROOT.TCut('expectedMissingInnerHits<='+str(wp.missingHitsBarrel))
  elif(region == 'endcap'):
    selectionCuts += ROOT.TCut('expectedMissingInnerHits<='+str(wp.missingHitsEndcap))
  else:
    selectionCuts += ROOT.TCut('expectedMissingInnerHits<='+str(wp.missingHitsExtend))
  return selectionCuts

def makeSubDirs(fileName):
  listOfSubDirs = fileName.split('/')[:-1]
  dirName = os.path.join(*listOfSubDirs)
  try:    os.makedirs(dirName)
  except: pass
  for i in range(2,len(listOfSubDirs)+1): shutil.copy("figures/index.php", os.path.join(*listOfSubDirs[:i]))
  return fileName
