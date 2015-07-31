from numpy import *
import numpy as np
import numpy.random as npr
import argparse
import os
import model_builder as mdb
import matplotlib
import matplotlib.pyplot as plt

####################CALCULATION####################
def calc_avg_std(steps):
    cwd = args.cwd
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
    try:
        np.savetxt('['+str(args.brs)+']std_model_param_step'+xp+'.dat',stdstep)
        np.savetxt('['+str(args.brs)+']avg_model_param_step'+xp+'.dat',avgstep)
    except:
        np.savetxt('['+str(args.brs)+']std_model_param_step.dat',stdstep)
        np.savetxt('['+str(args.brs)+']avg_model_param_step.dat',avgstep)

    os.chdir(cwd)

    return avgstep,stdstep
####################CALCULATION####################

####################PLOTS####################
def plot(avgstep,stdstep,xp):
    cwd = args.cwd
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
    try:
        fig.savefig('['+str(args.brs)+']plots_of_avg_and_std'+xp+'.png')
    except:
        fig.savefig('['+str(args.brs)+']plots_of_avg_and_std.png')

    print 'plotting figures, std as percent difference of avg'
    percents = []
    for i in xrange(num_param):
        per = abs(stdstep[i]/avgstep[i])
        percents.append(per)
    pfig = plt.figure()
    px = pfig.add_subplot(111,title='Std/Avg Model Param Step')
    px.hist(percents, bins=100)
    try:
        pfig.savefig('['+str(args.brs)+']plot_of_std-avg'+xp+'.png')
    except:
        pfig.savefig('['+str(args.brs)+']plot_of_std-avg.png')

    print 'plotting figures, std as percent difference of avg FILTERED'
    percents = []
    for i in xrange(num_param):
        per = abs(stdstep[i]/avgstep[i])
        if per > 200:
            per = 201
            percents.append(per)
        else:
            percents.append(per)
    pfig = plt.figure()
    px = pfig.add_subplot(111,title='Std/Avg Model Param Step')
    px.hist(percents, bins=100)
    try:
        pfig.savefig('['+str(args.brs)+']plot_of_std-avg'+xp+'FILTERED.png')
    except:
        pfig.savefig('['+str(args.brs)+']plot_of_std-avgFILTERED.png')

    print 'plotting histogram of just the std'
    stdfig = plt.figure()
    stdx = stdfig.add_subplot(111,title='Hist of Std')
    stdx.hist(stdstep,bins=100)
    try:
        stdfig.savefig('['+str(args.brs)+']hist_std'+xp+'.png')
    except:
        stdfig.savefig('['+str(args.brs)+']hist_std.png')

    print 'plotting histogram of just the averages'
    avgfig = plt.figure()
    avgx = avgfig.add_subplot(111,title='Hist of Avg')
    avgx.hist(avgstep,bins=100)
    try:
        avgfig.savefig('['+str(args.brs)+']hist_avg'+xp+'.png')
    except:
        avgfig.savefig('['+str(args.brs)+']hist_avg.png')

    os.chdir(cwd)
####################PLOTS####################

####################STEP####################
def step():
    cwd = args.cwd
    print 'loading file'
    steps = np.loadtxt('param_steps_'+str(args.brs)+'.dat')
    steps.shape = (args.brs,-1)
    print 'adjusting shape'
    (num_samples, num_param) = steps.shape
    print steps.shape

    sample = steps[:,args.step]
    print 'Should be %s values'%(num_samples)
    print sample.size

    print 'plotting histogram'
    hfig = plt.figure()
    hx = hfig.add_subplot(111,title='Histogram of Step Values for Param %s'%(args.step))
    hx.hist(sample,bins=100)
    hfig.savefig('hist_param_%s['+str(args.brs)+'].png'%(args.step))

    print 'plotting scatter plot'
    sfig = plt.figure()
    sx = sfig.add_subplot(111,title='Scatter Plot of Step Values for Param %s'%(args.step))
    xval = range(sample.size)
    for i in xrange(sample.size):
        sx.scatter(xval[i],sample[i])
    sfig.savefig('scatter_param_%s['+str(args.brs)+'].png'%(args.step))

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

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--brs','--resamplings',type=int,help='number of samples in the file',default=1000)
    parser.add_argument('--path',default=os.getcwd())
    parser.add_argument('--savelocation',default=os.getcwd())
    parser.add_argument('--xp',type=int,help='number of xp files',default=1)
    parser.add_argument('--cwd',default=os.getcwd())
    args = parser.parse_args()

    return args



if __name__ == '__main__':
    args = get_args()

    print 'Loading file'
    for i in xrange(args.xp):
        os.chdir(args.path)

        print os.getcwd()
        if args.xp ==1:
            xp = None
        else:
            xp = ',xp%s'%(i+1)
        print xp
        
        try:
            steps = np.loadtxt('param_steps_'+str(args.brs)+xp+'.dat')
        except:
            steps = np.loadtxt('param_steps_'+str(args.brs)+'.dat')

        print steps.shape
        steps.shape = (args.brs,-1)
        (num_samp,num_param) = steps.shape
        print num_samp
        print num_param

        os.chdir(args.cwd)
        average,standard = calc_avg_std(steps)
        plot(average,standard,xp)

    print 'COMPLETE'
exit()
