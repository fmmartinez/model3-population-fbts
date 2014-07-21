#!/usr/bin/python
import os
import shutil
from sys import argv

#argument provided is the name of the folder, which in turn has all the
#information needed for the script
if len(argv) == 1:
	print "you must provide name of simulation (folder name) as an argument"
	quit()

dataname = argv[1]

nproc = int(dataname[8:11])

dirs = []

if nproc > 9:
	if nproc > 99:
		if nproc > 999:
			print "error: too many processsors will be used"
		else:
			for i in range(0,10):
				name = 'part-00' + str(i)
				dirs.append(name)
			for i in range(10,100):
				name = 'part-0'  + str(i)
				dirs.append(name)
			for i in range(100,nproc):
				name = 'part-'  + str(i)
				dirs.append(name)
	else:
		for i in range(0,10):
			name = 'part-00' + str(i)
			dirs.append(name)
		for i in range(10,nproc):
			name = 'part-0'  + str(i)
			dirs.append(name)
else:
	for i in range(0,nproc):
		name = 'part-00'  + str(i)
		dirs.append(name)

currentpath = os.getcwd()
for i in range(0,nproc):
	workingd = currentpath + '/' + dataname + '/'  + dirs[i]
	#print workingd
	os.chdir(workingd)
	os.system('qsub submit_cluster.pbs')
