import os
import glob
pypath = 'F:\Schmithorst\CLINICAL\ClinicalfMRI\Python'
funcpath = 'F:\Schmithorst\CLINICAL\ClinicalfMRI\ABIGEAL_LABAHN'
os.chdir(pypath)
import fmriproc
os.chdir(funcpath)

timedir = '.\TIME_FILES\\'

fmriproc.processall('SYNONYMS5MINMB',timedir+'5MinutesMB.time')
fmriproc.processall('VISUAL5MINMB',timedir+'5MinutesMB.time')

fmriproc.processmotor('MOTOR6MINMB',timedir+'LeftHandMB.time',timedir+'RightHandMB.time')

logfile = glob.glob("*Verbs*.log")
logfile = logfile[0]
fmriproc.verbslogtotime(logfile,'Verbs.time')
fmriproc.hushprocessall('VERBSSILENTGRADIENTMB','Verbs.time',4)

fmriproc.hushprocessall('STORIESSILENTGRADIENTMB',timedir+'StoriesMB.time',3)


