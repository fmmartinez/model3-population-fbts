#!/usr/bin/python
import os
import shutil

#Machine stuff

#modules to load for execution in certain clusters, leave blank if none
#following is used for jasper
module = 'compiler/intel/13.0.1'
#module = ''

#TORQUE will use the following email to inform details about your jobs
email = 'fmmartin@ualberta.ca'

#number of processors to use
nproc = 100

#trajectories per processor
tpp = 5000

timestep = '5d-5'

#Wall time
walltime = '24:00:00'

#Model stuff

#energies
eg = 0
eb = 240
ed = 240

#laser thingies
e0 = 0.9
time = 0.15
tau = 0.045
omega = 260

#coupled vibration frequency
vibomega = 37.7

#basis functions used
g = 1  
b = 1  
d = 1  

#delta
delta = 1

#end of basic input variables
#edit below if you are sure what are you doing
#############################################################################

#not treating more than 100 basis functions per center, this can change in the future
if (g >= 100 or b >= 100 or d >= 100):
	print "error: too many basis functions used per center"
	quit()

#folder tree generation with consistent length
if g < 10:
	gpart = '0' + str(g)
else:
	gpart = str(g)

if b < 10:
	bpart = '0' + str(b)
else:
	bpart = str(b)

if d < 10:
	dpart = '0' + str(d)
else:
	dpart = str(d)

if nproc > 9:
	if nproc > 99:
		if nproc > 999:
			print "error: too many processsors will be used"
		else:
			cpus = str(nproc)
	else:
		cpus = '0' + str(nproc)
else:
	cpus = '00' + str(nproc)

gendirname = 'p' + gpart + bpart + dpart + '-' + cpus + '/'

#genenrate global folder
if not(os.path.exists(gendirname)):
	os.mkdir(gendirname)

dirs = []

if nproc > 9:
	if nproc > 99:
		if nproc > 999:
			print "error: too many processsors will be used"
		else:
			for i in range(0,10):
				name = gendirname + 'part-00' + str(i)
				dirs.append(name)
			for i in range(10,100):
				name = gendirname + 'part-0'  + str(i)
				dirs.append(name)
			for i in range(100,nproc):
				name = gendirname + 'part-'  + str(i)
				dirs.append(name)
	else:
		for i in range(0,10):
			name = gendirname + 'part-00' + str(i)
			dirs.append(name)
		for i in range(10,nproc):
			name = gendirname + 'part-0'  + str(i)
			dirs.append(name)
else:
	for i in range(0,nproc):
		name = gendirname + 'part-00'  + str(i)
		dirs.append(name)

#generate folders inside global folder
for i in range(0,nproc):
	if not(os.path.exists(dirs[i])):
			os.mkdir(dirs[i])

# input files generation
l = []
m = []

#skeleton of md.in file:
l.append('Np\tDELTA\tNOSC\tOME_MAX\t\r\n')
l.append('1\t' + str(delta) + 'd0\t20\t50\r\n')
l.append('NMCS\tNMDS\tseed\tDT\tLAMDA_D\r\n')
l.append('')
l.append('Eg\tEb\tEd\tmu\tE0\tbeta\tvib_omega\r\n')
l.append(str(eg) + '\t' + str(eb) + '\t' + str(ed) + '\t1\t' + str(e0) + '\t0.24\t' + str(vibomega) + '\r\n')
l.append('TIME_J\tTAU_J\tOMEGA_J\r\n')
l.append(str(time) + '\t' + str(tau) + '\t' + str(omega) + '\r\n')
l.append('BATH(0:B EQ 1:T EQ)\tINIT\r\n')
l.append('0\t3\r\n')
l.append('mapinG\tmapinB\tmapinD\r\n')
l.append(str(g) + '\t' + str(b)+ '\t' + str(d) + '\r\n')

#creating md.in files
for i in range(0,nproc):
	mdfile = open('./' + dirs[i] + '/md.in','w')
	
	seed = i*9 + 5

	l[3] = str(tpp) + '\t4000\t' + str(seed) + '\t' + timestep + '\t10d0\r\n'

	for j in range(0,12):
		mdfile.write(l[j])

	mdfile.close()

#skeleton of pbs fie
m.append('#!/bin/bash -l\n')
m.append('#PBS -r n\n')
m.append('')
m.append('#PBS -l walltime=' + walltime + '\n')
m.append('#PBS -l procs=1\n')
m.append('#PBS -l mem=256mb\n')
m.append('')
m.append('#PBS -M' + email + '\n')
m.append('cd $PBS_O_WORKDIR\n')
m.append('module load ' + module + '\n')
m.append('time ./a.out < md.in\n')

#creating pbs files
for i in range(0,nproc):
	pbsfile = open('./' + dirs[i] + '/submit_cluster.pbs','w')
	
	m[2] = '#PBS -N p' + str(g) + str(b) + str(d) + '-' + str(i) + '\n'
	
	if ((i == 0)|(i == nproc-1)):
		m[6] = '#PBS -m bea\n'
	else:
		m[6] = '#PBS -m n\n'

	for j in range(0,11):
		pbsfile.write(m[j])

	pbsfile.close()


#copy executables
shutil.copy2('merge.out',gendirname)

for i in range(0,nproc):
	shutil.copy2('a.out',dirs[i])
#	shutil.copy2('m_vib.mod',dirs[i])
