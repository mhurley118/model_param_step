from numpy import *
import numpy as np
import numpy.random as npr
import argparse
import os
import model_builder as mdb
import matplotlib
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('--brs',type=int,help='number of samples in the file',default=1000)
parser.add_argument('--path',default='data_files')
parser.add_argument('--step',default=None,type=int)
args = parser.parse_args()

cwd = os.getcwd()
os.chdir(args.path)

print 'Loading file'
steps = np.loadtxt('param_steps_'+str(args.brs)+'.dat')
print steps.shape
steps.shape = (args.brs+1,-1)
(num_samp,num_param) = steps.shape
print num_samp
print num_param
os.chdir(cwd)

####################CALCULATION####################
stdstep = []
avgstep = []
print 'Calculateing std and avg'
for i in xrange(num_param):
    svalue = std(steps[:,i])
    avalue = mean(steps[:,i])
    stdstep.append(svalue)
    avgstep.append(avalue)

print 'Verifying all parameters have been calculated, SHOULD BE %s INDICES LONG'%(num_param)
print len(stdstep)
print len(avgstep)
if len(stdstep) != num_param:
    print 'Error: Not all parameters were calculated for std'
    exit()
elif len(avgstep) != num_param:
    print 'Error: Not all parameters were calculated for avg'
    exit()
else:
    pass

try:
   os.mkdir('data_files')
except:
   pass
os.chdir(args.path)
print 'saving lists as data files'
np.savetxt('['+str(args.brs)+']std_model_param_step.dat',stdstep)
np.savetxt('['+str(args.brs)+']avg_model_param_step.dat',avgstep)
os.chdir(cwd)
####################CALCULATION####################

####################PLOTS####################
try:
    os.mkdir('plots')
except:
    pass
os.chdir('plots')

print 'plotting figures, std and avg'
xval = range(num_param)
fig = plt.figure()
sx = fig.add_subplot(211,title='Standard Deviations of Model Param Steps')
ax = fig.add_subplot(212,title='Averages of Model Param Steps')
sx.scatter(xval,stdstep)
ax.scatter(xval,avgstep)
fig.savefig('['+str(args.brs)+']plots_of_avg_and_std.png')

print 'plotting figures, std as percent difference of avg'
percents = []
for i in xrange(num_param):
    per = abs(stdstep[i]/avgstep[i])
    percents.append(per)
pfig = plt.figure()
px = pfig.add_subplot(111,title='Std/Avg Model Param Step')
px.hist(percents, bins=100)
pfig.savefig('['+str(args.brs)+']plot_of_std-avg.png')

print 'plotting histogram of just the std'
stdfig = plt.figure()
stdx = stdfig.add_subplot(111,title='Hist of Std')
stdx.hist(stdstep,bins=100)
stdfig.savefig('['+str(args.brs)+']hist_std.png')

print 'plotting histogram of just the averages'
avgfig = plt.figure()
avgx = avgfig.add_subplot(111,title='Hist of Avg')
avgx.hist(avgstep,bins=100)
avgfig.savefig('['+str(args.brs)+']hist_avg.png')
####################PLOTS####################

####################STEP####################
if args.step:
    os.chdir('%s/%s'%(cwd,args.path))
    print 'loading file'
    steps.shape = (args.brs+1,-1)
    print 'adjusting shape'
    (num_samples, num_param) = steps.shape
    print steps.shape

    stdx = stdstep[args.step]
    sample = steps[:,args.step]
    print 'Should be %s values'%(num_samples)
    print sample.size

    print 'plotting histogram'
    plt.title('Histogram of Step Values for Param #%s'%(args.step),fontsize=22)
    plt.ylabel('Count',fontsize=20)
    plt.xlabel('Change',fontsize=20)
    plt.hist(sample,bins=100,normed=False,color='blue')
    os.chdir('%s/plots'%cwd)
    plt.savefig('hist_param_%s[%s].png'%(args.step,args.brs))

#    print 'plotting scatter plot'
#    sfig = plt.figure()
#    sx = sfig.add_subplot(111,title='Scatter Plot of Step Values for Param %s'%(args.step))
#    xval = range(sample.size)
#    for i in xrange(sample.size):
#        sx.scatter(xval[i],sample[i])
#    sfig.savefig('scatter_param_%s[%s].png'%(args.step,args.brs))


print 'plotting figures, std as percent difference of avg FILTERED'
percents = []
for i in xrange(num_param):
    per = abs(stdstep[i]/avgstep[i])
    percents.append(per)
for i in xrange(len(percents)):
    if percents[i] > 200:
        percents[i] = 201
pfig = plt.figure()
px = pfig.add_subplot(111,title='Std/Avg Model Param Step FILTERED 200')
px.hist(percents, bins=100)
pfig.savefig('['+str(args.brs)+']plot_of_std-avgFILTERED.png')




print 'COMPLETE'
exit()
