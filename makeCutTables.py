#! /usr/bin/env python

import ROOT,os,glob,shutil
from common      import loadClasses, workingPoints, compareWP
from collections import OrderedDict
#loadClasses('VarCut.cc', 'OptimizationConstants.hh')

latexVars = {'full5x5_sigmaIetaIeta' : 'full $5\\times5 \\sigma_{i\\eta i\\eta} <$',
             'dEtaSeed'              : '$|$dEtaInSeed$|$ $<$',
             'dPhiIn'                : '$|$dPhiIn$|$ $<$',
             'hOverE'                : '$H/E <$',
             'hOverEscaled'          : '$C_0 <$',
             'relIsoWithEA'          : 'relIsoWithEA $<$',
             'ooEmooP'               : '$|1/E-1/p| <$',
             'missingHits'           : 'expected missing inner hits $<=$',
             'conversionVeto'        : 'pass conversion veto'}

twikiVars = {'full5x5_sigmaIetaIeta' : 'full5x5_sigmaIetaIeta <',
             'dEtaSeed'              : 'abs(dEtaSeed) <',
             'dPhiIn'                : 'abs(dPhiIn) <',
             'hOverE'                : 'H/E <',
             'hOverEscaled'          : 'C_0 <',
             'relIsoWithEA'          : 'relIsoWithEA <',
             'ooEmooP'               : 'abs(1/E-1/p) <',
             'missingHits'           : 'expected missing inner hits <=',
             'conversionVeto'        : 'pass conversion veto'}

def writeTwikiLine(f, *args): f.write(('|  %-40s  ' + '|  %20s  '*(len(args)-1) + '| \n')     % args) 
def writeLatexLine(f, *args): f.write(('   %-40s  ' + '&  %20s  '*(len(args)-1) + ' \\\\ \n') % args)

def writeTwikiTable(f, cutValues):
  writeTwikiLine(f, '', *(wp.name for wp in cutValues.keys()))
  for var in cutValues[cutValues.keys()[0]]:
    writeTwikiLine(f, twikiVars[var], *(cutValues[wp][var] for wp in cutValues.keys()))

def writeLatexTable(f, cutValues):
  f.write('\\begin{tabular}{l|' + 'ab'*int(len(cutValues)/2) + 'a'*(len(cutValues)%2) +  '} \n')
  writeLatexLine(f, '', *(wp.name for wp in cutValues.keys()))
  f.write('\\hline \n')
  for var in cutValues[cutValues.keys()[0]]:
    writeLatexLine(f, latexVars[var], *(cutValues[wp][var] for wp in cutValues.keys()))
  f.write('\\end{tabular}')

def makeTables(outFileName, wps):
  for barrel in [False, True]:
    cutValues = OrderedDict()
    for wp in wps:
      file = ROOT.TFile('cut_repository/' + (wp.cutsFileBarrel if barrel else wp.cutsFileEndcap) + '.root')
      cuts = file.Get('cuts')
      
      cutValues[wp] = OrderedDict()
      for var in ["full5x5_sigmaIetaIeta", "dEtaSeed", "dPhiIn", "hOverEscaled" if outFileName.count('etune') else "hOverE", "relIsoWithEA", "ooEmooP"]:
        cutValues[wp][var] = '%.3g' % cuts.getCutValue(var.replace('scaled',''))
      cutValues[wp]['missingHits'] = str(wp.missingHitsBarrel if barrel else wp.missingHitsEndcap)
      cutValues[wp]['conversionVeto'] = 'yes'

    with open('tables/' + outFileName + '-' + ('barrel' if barrel else 'endcap') + '.txt', 'w') as f: writeTwikiTable(f, cutValues)
    with open('tables/' + outFileName + '-' + ('barrel' if barrel else 'endcap') + '.tex', 'w') as f: writeLatexTable(f, cutValues)

try:    os.makedirs('tables')
except: pass
#for tag in ['default', 'retuneMVA', 'retuneEff','2016', 'compare2016','compareRetuneEff']:
#  makeTables(tag, workingPoints[tag])


def compileLatex():
  for f in glob.glob('tables/*.tex'):
    os.system("pdflatex --jobname=temp '\\input{latexDefinitions}\\begin{document}\\resizebox{0.6\\textwidth}{!}{\\input{" + f + "}}\\end{document}' &> /dev/null")
    os.system("pdfcrop temp.pdf temp.pdf &> /dev/null")
    os.system("convert -density 1000 -resize 2000 temp.pdf temp.png")
    shutil.copy('temp.pdf', f.replace('.tex','.pdf'))
    shutil.copy('temp.png', f.replace('.tex','.png'))
    os.system("rm temp*")

shutil.copy('figures/index.php', 'tables')
compileLatex()
