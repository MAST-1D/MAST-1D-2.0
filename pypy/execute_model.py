#  By Katie
#  This adds the MAST-1D and Output folders to the system path (so pypy knows to look/write there).
#  It then asks the user which model setup to run.
#  It then executes 'Inputs_mainpage' for the selected model via pypy.

import sys
import os

filelist = []
riverlist = []
riverchoice = {}

modelfolder = str(os.path.join(os.pardir, "Input_pages"))
outputfolder = str(os.path.join(os.pardir, "Output"))

files = os.listdir(modelfolder)
for f in files:
	if 'Inputs_mainpage_' in f:
		filelist.append(f)
		
for f in filelist:
	rivername = f[16:-3]
	riverlist.append([rivername,f])

print '\n'	
print "Which model would you like to run:"
print '\n'

i = 0
for river in riverlist:
	print "for " + river[0] + ' press ' + str(i)
	riverchoice[str(i)] = river[1]
	i = i + 1

print '\n'
	
runfile = raw_input()
modelfile = str(os.path.join(os.pardir, "Input_pages", riverchoice[runfile]))
sys.path = [modelfolder,outputfolder] + sys.path

print '\n'
print 'Running ' + riverchoice[runfile]
print '\n'

execfile(modelfile)
