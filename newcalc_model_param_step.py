import project_tools.parameter_fitting.newton_solver.Truncated_SVD as truncated
import project_tools.parameter_fitting.FRET.compute_Jacobian as compute
import project_tools.parameter_fitting.util.util as util
import model_builder as mdb
from numpy import *
import numpy as np
import numpy.random as npr
import argparse
import os
import mdtraj as md

parser = argparse.ArgumentParser()
parser.add_argument('--brs',default=1000,type=int,help='number of resamplings')
args = parser.parse_args()

cwd = os.getcwd()
model, fitopts = mdb.inputs.load_model('1PB7',dry_run=True)
global GAS_CONSTANT_KJ_MOL
GAS_CONSTANT_KJ_MOL = 0.0083144621
def_FRET_pairs = [[219, 371]]
defspacing = 0.1

print 'TRUNCATE VALUE IS:'
if 'truncate_value' in fitopts:
    truncate_value = fitopts['truncate_value']
else:
    truncate_value = 0.01
print truncate_value

cwd = os.getcwd()
subdir = model.name
iteration = fitopts["iteration"]
fit_temp = fitopts["t_fit"]
sub = "%s/%s/iteration_%d" % (cwd,subdir,iteration)
subtemp = "%s/%d_0" % (sub,fit_temp)
os.chdir(subtemp)
traj = md.load("traj.xtc",top="Native.pdb")
print 'computing original FRETr distances'
oFRETr = md.compute_distances(traj,def_FRET_pairs, periodic=False)
os.chdir(cwd)

#############################COMPUTE_JACOBIAN##############################
def get_sim_params(model,fitopts,FRETr):
    fit_temp = fitopts["t_fit"]
    cwd = os.getcwd()
    subdir = model.name
    iteration = fitopts["iteration"]
    
    sub = "%s/%s/iteration_%d" % (cwd,subdir,iteration)
    subdirec = "%s/fitting_%d" % (sub,iteration)
    parmfile = "%s/simf-params%d.dat" % (subdirec,fit_temp)
    if not os.path.isfile(parmfile):
        find_sim_bins(subdirec, FRETr[:,0], fit_temp, residues=residues, spacing=spacing, weights=weights) 
    parms = np.loadtxt(parmfile)
    num_bins, ran_size, spacing = compute.analyze_sim_params(parms)
    return num_bins, ran_size, spacing

def get_target_feature(model,fitopts,FRETr):
    fit_temp = fitopts["t_fit"]
    cwd = os.getcwd()
    subdir = model.name
    iteration = fitopts["iteration"]
    
    sub = "%s/%s/iteration_%d" % (cwd,subdir,iteration)
    subdirec = "%s/fitting_%d" % (sub,iteration)
    simfile = "%s/simf_centers%d.dat" % (subdirec,fit_temp)
    FRETfile = "%s/FRET_hist.dat" % subdirec
    
    bin_centers = compute.get_sim_centers(model, fitopts)
    bin_size, ran_size, spacing = get_sim_params(model, fitopts,FRETr)
    ##Re-round to utilize the same range for ran_size as used in the simparams, in case floating point is round-off errored
    ran_size = (int(round(ran_size[0]/spacing))*spacing, int(round(ran_size[1]/spacing))*spacing)
    if not os.path.isfile(FRETfile):
        compute.fret_hist_calc(model, fitopts, bin_size, ran_size, spacing)  
    FRETdata = np.loadtxt(FRETfile)  
    
    print "initial FRET data and bin_centers"

    print FRETdata[:,0]
    print bin_centers
    
    if compute.check_exp_data(FRETdata[:,0], bin_centers):
        print "Mismatched experimental data and simulated data. Attempting Re-binning"
        compute.fret_hist_calc(model, fitopts, bin_size, ran_size, spacing)
        FRETdata = np.loadtxt(FRETfile)
        if compute.check_exp_data(FRETdata[:,0], bin_centers):
            compute.add_error_log("Catastrophic miscalculation of FRET bins", fit_temp)
            print "Found the FRETdata and bin_centers to be:"
            print FRETdata[:,0]
            print bin_centers
    
    target = FRETdata[:,1]
    target_err = target**0.5 ##for lack of a better way, take sqrt of bins for error estimate
    
    return target, target_err

def calculate_average_Jacobian(model,fitopts,FRETr,FRET_pairs=def_FRET_pairs,spacing=defspacing):
    print "Working on calculating model's trajectory and contact info"
    if "t_fit" in fitopts:
        fit_temp = fitopts["t_fit"]
    else:
        raise IOError("Missing the fit_temperature, please specify in .ini file")
    if "fret_pairs" in fitopts:
        fret_pairs = fitopts["fret_pairs"]
        FRET_pairs = np.array(fret_pairs) - 1
        print "The FRET pairs are:"
        print FRET_pairs
    if "y_shift" in fitopts:
        y_shift = fitopts["y_shift"]   
    else:
        y_shift = 0.0
        fitopts["y_shift"] = 0.0        
    if "spacing" in fitopts:
        spacing = fitopts["spacing"]

    cwd = os.getcwd()
    subdir = model.name
    iteration = fitopts["iteration"]
    sub = "%s/%s/iteration_%d" % (cwd,subdir,iteration)
    traj_location = "%s/%d_0" % (sub, fit_temp)
    sim_location = "%s/fitting_%d" % (sub,iteration)
    ##define location of logical files
    os.chdir(traj_location)
    traj,rij,qij = util.get_rij_Vp(model)

    sim_feature, sim_slices = compute.find_sim_bins(sim_location, FRETr[:,0], fit_temp, residues=FRET_pairs, spacing=spacing, weights=None)
    
    beta = 1.0 / (GAS_CONSTANT_KJ_MOL*float(fit_temp))
    os.chdir(cwd)
    
    print "Computing Jacobian and Simparams for the temperature %d, with spacing %f" % (fit_temp, spacing)
    Jacobian = compute.compute_Jacobian_basic(qij,sim_feature*spacing, sim_slices, beta)
    Jacobian /= spacing        
    sim_feature_err = sim_feature ** 0.5
    Jacobian_err = np.zeros(np.shape(Jacobian))
    
    return sim_feature, sim_feature_err, Jacobian, Jacobian_err
#############################COMPUTE JACOBIAN##############################

#############################FUNCTIONS##############################
def calc_Jacobian(FRETr):
    sim_feature, sim_feature_err, Jacobian, Jacobian_err = calculate_average_Jacobian(model,fitopts,FRETr,FRET_pairs=def_FRET_pairs,spacing=defspacing)
    target_feature, target_feature_err = get_target_feature(model,fitopts,FRETr)

    np.savetxt('target_feature.dat',target_feature)
    np.savetxt('sim_feature.dat',sim_feature)
    np.savetxt('Jacobian.dat', Jacobian)

def calc_step():
    print 'Calcualting steps'
    truncated.find_solutions(model,scaling=False)

def find_step():
    print'Loading lambdas'
    lambdas = np.loadtxt('lambdas.dat')
    
    real_trunc_value = min(lambdas, key=lambda x:abs(x-truncate_value))
    print 'REAL TRUNCATED VALUE IS:'
    print real_trunc_value
    print 'MODEL PARAM STEPS IN FILE:'
    for i in xrange(lambdas.size):
        if lambdas[i] == real_trunc_value:
            xp_num = i
    print 'xp_%s.dat'%(xp_num)
    xp_file = np.loadtxt('xp_%s.dat'%(xp_num))
    return xp_file.tolist()

def save_step(xp_file):
    os.chdir('data_files')
    print 'Saving data to file'
    stepfile = open('param_steps_'+str(args.brs)+'.dat','a')
    for i in xrange(len(xp_file)):
        stepfile.write('%s '%(xp_file[i]))
    os.chdir(cwd)
#############################FUNCTIONS##############################

#############################FIRST_RUN##############################
print 'Creating initial Jacobian files'
calc_Jacobian(oFRETr)

print 'Calculating param step'
calc_step()
xp_file = find_step()
print 'Saving data to file'
try:
    os.mkdir('data_files')
except:
    pass
os.chdir('data_files')
stepfile = open('param_steps_'+str(args.brs)+'.dat','a')
for i in xrange(len(xp_file)):
    stepfile.write('%s '%(xp_file[i]))
os.chdir(cwd)
#############################FIRST_RUN##############################

#############################BOOTSTRAPPING##############################
print 'Bootstrapping'
(x,y) = oFRETr.shape
oFRETr = oFRETr.ravel()
failedruns = []
for i in xrange(args.brs):
    try:
        FRETr = npr.choice(oFRETr,oFRETr.size)
        FRETr.shape = (x,y)
        calc_Jacobian(FRETr)

        calc_step()
        xp_file = find_step()
        save_step(xp_file)

        print 'NUMBER OF MODEL PARAM STEPS CALCULATED'
        print str(i+2)
    except:
        pass
#############################BOOTSTRAPPING##############################
