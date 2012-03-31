"""
Module of routines to take the DART results from a GeoClaw simulation 
and compare to the DART observations, shifting by an optimal t0 and
computing the RMS.

The GeoClaw code must first be run and the gauge output converted to a
suitable form.
"""

import pylab
import numpy as np
from scipy import interpolate 

#gext = 'eps'  # graphics extension
gext = 'png'  # graphics extension

# source models:
models = """GCMT Hayes UCSB3 Ammon Caltech Fujii Saito Gusman GusmanAU pmelWei 
          """.split()

modelnum = {}
modelnum["GCMT"] = '1'
modelnum["Hayes"] = '2'
modelnum["UCSB3"] = '3'
modelnum["Ammon"] = '4'
modelnum["Caltech"] = '5'
modelnum["Fujii"] = '6'
modelnum["Saito"] = '7'
modelnum["Gusman"] = '8a'
modelnum["GusmanAU"] = '8b'
modelnum["pmelWei"] = '9'

#models = ["pmelWei"]

gaugenos = [21401, 21413, 21418, 21419]

pTd = 0  # p in weighted norm for minimizing RMS in finding Td
p0 = 0   # p in weighted norm for computing RMS to plot
         # possibly both equal to 0.5?

def make_interp_fcn(t,eta):
    """
    Create and return a piecewise linear interpolating function
    based on the discrete data t,eta.
    """
    if t[0] > t[-1]:
        t = np.flipud(t)
        eta = np.flipud(eta)
    f = interpolate.interp1d(t,eta,bounds_error=False,fill_value=0.)
    return f

def get_data(gaugeno, model):

    #fname = 'DART/%s_notide.txt' % gaugeno
    fname = '%s_notide.txt' % gaugeno
    t_dart,eta_dart = np.loadtxt(fname, unpack=True)
    eta_dart_fcn = make_interp_fcn(t_dart,eta_dart)

    #fname = 'Simulations/%s%s.txt' % (model,gaugeno)
    fname = 'simulation_results/%s_%s.txt' % (model,gaugeno)
    gaugedata = np.loadtxt(fname)
    #t_sim = gaugedata[:,2]
    #eta_sim = gaugedata[:,6]
    t_sim = gaugedata[:,0]
    eta_sim = gaugedata[:,1]
    eta_sim_fcn = make_interp_fcn(t_sim,eta_sim)
    return t_dart,eta_dart,eta_dart_fcn, t_sim,eta_sim,eta_sim_fcn
    
def error(t,eta_dart_fcn,eta_sim_fcn,t0,t1,t2,p):
    eta_dart = eta_dart_fcn(t)
    eta_sim = eta_sim_fcn(t-t0)
    err = eta_sim - eta_dart
    weighted_err = err * abs(eta_dart)**p
    rms = np.sqrt(sum(weighted_err[(t>=t1) & (t<=t2)]**2))
    if 0:
        pylab.figure(2)
        pylab.plot([t0],[rms],'g^')
        pylab.figure(10)
        pylab.clf()
        pylab.plot(t,err,'b')
        pylab.plot(t,weighted_err,'r')
        pylab.title('p = %s, t0 = %s,  rms = %s' % (p,t0,rms))
        pylab.draw()
        #raw_input("Hit return...")
    return rms
   
def plot_rms(t,eta_dart_fcn,eta_sim_fcn,t0a,t0b,t1,t2,p):
    t0_vals = np.linspace(t0a,t0b,20)
    rms = []
    for t0 in t0_vals:
        rms.append(error(t,eta_dart_fcn,eta_sim_fcn,t0,t1,t2,p))
    rms = np.array(rms)
    pylab.figure(2)
    pylab.clf()
    pylab.plot(t0_vals,rms,'ob')
    return rms

def minimize_rms(t,eta_dart_fcn,eta_sim_fcn,t0_initial,t1,t2,p):
    from scipy.optimize import fsolve
    def error1(t0):
        return error(t,eta_dart_fcn,eta_sim_fcn,t0,t1,t2,p)
    t0_min = fsolve(error1, t0_initial)
    if len(t0_min) > 1:
        print "*** Warning, more than one root found: ",t0_min
    t0_min = t0_min[0]
    rms_min = error(t,eta_dart_fcn,eta_sim_fcn,t0_min,t1,t2,p)
    #print "+++ p= ",p,t0_min,t1,t2,len(t)
    #print "+++ rms_min= ",rms_min
    t0a = t0_min - 120
    t0b = t0_min + 120
    plot_rms(t,eta_dart_fcn,eta_sim_fcn,t0a,t0b,t1,t2,p)
    pylab.plot([t0_min],[rms_min],'or')

    pylab.figure(3)
    pylab.clf()
    pylab.plot(t,eta_dart_fcn(t),'k')
    pylab.plot(t,eta_sim_fcn(t),'b')
    pylab.plot(t,eta_sim_fcn(t-t0_min),'r')
    pylab.legend(['DART','S(t)','S(t-t0)'])
    #print "+++ rms_min= ",rms_min
    return t0_min, rms_min

    
def solve_and_plot(gaugeno,model,t1,t2):
    t_dart,eta_dart,eta_dart_fcn, t_sim,eta_sim,eta_sim_fcn = \
         get_data(gaugeno,model)
    t0_initial = 0.
    p = pTd
    t = np.arange(t1-500,t2+500,15)
    t0_min, rms_min = minimize_rms(t,eta_dart_fcn,eta_sim_fcn,t0_initial,t1,t2,p)
    #print "+++ rms_min= ",rms_min

    # Replot over tspan hours from time t1:
    tspan = 8
    pylab.figure(3)
    pylab.clf()
    tplot = np.arange(t1,t1+tspan*3600.,60)
    pylab.plot(tplot,eta_dart_fcn(tplot),'k',linewidth=2)
    pylab.plot(tplot,eta_sim_fcn(tplot),'b',linewidth=2)
    pylab.plot(tplot,eta_sim_fcn(tplot-t0_min),'r',linewidth=2)
    pylab.legend(['DART','S(t)','S(t-t0)'])
    pylab.title("DART %s with source %s:  t0 = %7.2f" \
              % (gaugeno,modelnum[model],t0_min))
    pylab.xlim([t1,t2])
    pylab.xlabel("Seconds")
    pylab.ylabel("Elevation (m)")

    # Recompute RMS over tspan hours with p=0, both shifted by t0 and
    # unshifted:
    #tspan = 1.5*3600.
    #tspan = 2.0*3600.
    #rms_min = error(t,eta_dart_fcn,eta_sim_fcn,t0_min,t1,t1+tspan,p0)
    #print "+++ p= ",p0,t0_min,t1,t1+tspan,len(t)
    #print "+++ new rms_min= ",rms_min

    rms_shift = rms_min
    rms_noshift = error(t,eta_dart_fcn,eta_sim_fcn,0,t1,t2,p0)

    #print "+++ rms_noshift = ",rms_noshift
    #print "tmax = ", t.max()

    # compute RMS for null solution:
    def eta_nosim_fcn(t): 
        return 0
    rms_nosim = error(t,eta_dart_fcn,eta_nosim_fcn,0.,t1,t2,p0)
    #print "RMS with 0 function for simulation: ",rms_nosim

    return t0_min, rms_shift, rms_noshift, rms_nosim

def run_all():
    t0_dict = {}
    rms_shift_dict = {}
    rms_noshift_dict = {}
    rms_nosim_dict = {}

    for gaugeno in gaugenos:
        #if gaugeno==21401: t1 = 3000; t2 = 6000;
        #if gaugeno==21413: t1 = 4000; t2 = 7000;
        #if gaugeno==21418: t1 = 1400; t2 = 4400;
        #if gaugeno==21419: t1 = 4000; t2 = 8000;
        if gaugeno==21401: t1 = 3000; t2 = t1+2.*3600;
        if gaugeno==21413: t1 = 4000; t2 = t1+2.*3600;
        if gaugeno==21418: t1 = 1400; t2 = t1+2.*3600;
        if gaugeno==21419: t1 = 4000; t2 = t1+2.*3600;
        for model in models:
            t0_min, rms_shift, rms_noshift, rms_nosim \
                  = solve_and_plot(gaugeno,model,t1,t2)
            pylab.savefig('figs/%s%s.%s' % (model,gaugeno,gext))
            t0_dict[model,gaugeno] = t0_min
            rms_shift_dict[model,gaugeno] = rms_shift
            rms_noshift_dict[model,gaugeno] = rms_noshift
            #print t1, t2, rms_0
            rms_nosim_dict[gaugeno] = rms_nosim
            print "%s %s:  %7.2f  %9.5f  %9.5f" \
               % (model.rjust(9), gaugeno, t0_min, rms_shift, rms_noshift)


    print " "
    print "t0 Values:"
    format = "%s " + len(gaugenos)*"%9i  "
    print format % tuple([" ".rjust(9)] + gaugenos)
    format = "%s " + len(gaugenos)*"%9.2f  "
    for model in models:
        print format \
           % tuple([model.rjust(9)] + [t0_dict[model,gaugeno] for gaugeno in gaugenos])
    print " "
    print "RMS Values (shifted by t0):"
    format = "%s " + len(gaugenos)*"%9i  "
    print format % tuple([" ".rjust(9)] + gaugenos)
    format = "%s " + len(gaugenos)*"%9.4f  "
    print format % tuple(["FlatWater"] + [rms_nosim_dict[gaugeno] for gaugeno in gaugenos])
    for model in models:
        print format \
           % tuple([model.rjust(9)] + [rms_shift_dict[model,gaugeno] for gaugeno in gaugenos])

    print " "
    print "RMS Values (unshifted):"
    format = "%s " + len(gaugenos)*"%9i  "
    print format % tuple([" ".rjust(9)] + gaugenos)
    format = "%s " + len(gaugenos)*"%9.4f  "
    print format % tuple(["FlatWater"] + [rms_nosim_dict[gaugeno] for gaugeno in gaugenos])

    for model in models:
        print format \
           % tuple([model.rjust(9)] + [rms_noshift_dict[model,gaugeno] for gaugeno in gaugenos])
    return t0_dict, rms_shift_dict, rms_noshift_dict, rms_nosim_dict

def plot_rms_models(gaugeno,t0_dict,rms_shift_dict,rms_noshift_dict,rms_nosim_dict):
    pylab.figure(4)
    pylab.clf()
    rms_shift = [rms_shift_dict[model,gaugeno] for model in models]
    rms_noshift = [rms_noshift_dict[model,gaugeno] for model in models]
    #import pdb; pdb.set_trace()
    rms_nosim = rms_nosim_dict[gaugeno]
    pylab.plot(rms_noshift/rms_nosim,'bo',markersize=8)
    pylab.plot(rms_shift/rms_nosim,'r^',markersize=8)
    if gaugeno==21418:
        rmslim = (0,1.2)
    else:
        rmslim = (0,1.2)
    for i in np.arange(len(models)):
        pylab.plot([i,i],rmslim,'k')
    #r = rms_nosim_dict[gaugeno]
    r = 1.  # after normalizing
    pylab.plot([0,9],[r,r],'g',linewidth=2)

    #pylab.legend(['unshifted','shifted by t0'],'lower right')
    #pylab.xlim(-1,15)
    pylab.xlim(-1,10)

    modelnums = [str(i) for i in range(1,11)]
    modelnums[-1] = '9'
    modelnums[7] = '8a'
    modelnums[8] = '8b'
    pylab.xticks(range(10),modelnums)
    pylab.ylim(rmslim)
    pylab.title('Relative RMS error at %s' % gaugeno)

def plot_rms_models_all_guages(t0_dict,rms_shift_dict,rms_noshift_dict,rms_nosim_dict):
    pylab.figure(4)
    pylab.clf()
    markers = {}
    markers[21401] = 'bo'
    markers[21413] = 'r^'
    markers[21418] = 'g<'
    markers[21419] = 'k>'
    for gaugeno in gaugenos:
        rms_shift = [rms_shift_dict[model,gaugeno] for model in models]
        rms_noshift = [rms_noshift_dict[model,gaugeno] for model in models]
        #import pdb; pdb.set_trace()
        rms_nosim = rms_nosim_dict[gaugeno]
        #pylab.plot(rms_noshift/rms_nosim,'bo',markersize=8)
        pylab.plot(rms_shift/rms_nosim,markers[gaugeno],markersize=8)
    pylab.legend([str(gaugeno) for gaugeno in gaugenos])

    rmslim = (0,1.2)
    for i in np.arange(len(models)):
        pylab.plot([i,i],rmslim,'k')
    #r = rms_nosim_dict[gaugeno]
    r = 1.  # after normalizing
    pylab.plot([0,9],[r,r],'g',linewidth=2)

    #pylab.legend(['unshifted','shifted by t0'],'lower right')
    #pylab.xlim(-1,15)
    pylab.xlim(-1,10)

    modelnums = [str(i) for i in range(1,11)]
    modelnums[-1] = '9'
    modelnums[7] = '8a'
    modelnums[8] = '8b'
    pylab.xticks(range(10),modelnums)
    pylab.ylim(rmslim)
    pylab.title('Relative RMS error of shifted data')

def plot_t0(t0_dict):
    pylab.figure(5)
    pylab.clf()

    markers = {}
    markers[21401] = 'bo'
    markers[21413] = 'r^'
    markers[21418] = 'g<'
    markers[21419] = 'k>'
    for i,model in enumerate(models):
        for gaugeno in gaugenos:
            t0val = t0_dict[model,gaugeno]
            pylab.plot([i],[t0val],markers[gaugeno],markersize=8)
        if i==0:
            pylab.legend([str(gaugeno) for gaugeno in gaugenos])

    for i in np.arange(len(models)):
        pylab.plot([i,i],[-200,450],'k')
    pylab.xlim(-1,10)

    modelnums = [str(i) for i in range(1,11)]
    modelnums[-1] = '9'
    modelnums[7] = '8a'
    modelnums[8] = '8b'
    pylab.xticks(range(10),modelnums)
    pylab.ylim(-200,450)
    pylab.title('Optimal time shift')

def plot_all(t0_dict,rms_shift_dict,rms_noshift_dict,rms_nosim_dict):

    for gaugeno in gaugenos:
        plot_rms_models(gaugeno,t0_dict,rms_shift_dict,rms_noshift_dict,rms_nosim_dict)
        fname = 'figs/rms%s.%s' % (gaugeno,gext)
        pylab.savefig(fname)
        print "Created ",fname

    plot_t0(t0_dict)
    fname = 'figs/Td.%s' % gext
    pylab.savefig(fname)
    print "Created ",fname

    plot_rms_models_all_guages(t0_dict,rms_shift_dict,rms_noshift_dict,rms_nosim_dict)
    fname = 'figs/rms_alldart.%s' % gext
    pylab.savefig(fname)
    print "Created ",fname


    
def plot_dart_multimodel(gaugeno,models,t0_dict,t1,t2):
    pylab.figure(3)
    pylab.clf()

    color = ['r','b','g','m']
    for i,model in enumerate(models):
        t_dart,eta_dart,eta_dart_fcn, t_sim,eta_sim,eta_sim_fcn = \
             get_data(gaugeno,model)

        tplot = np.arange(t1,t2,30)
        #pylab.plot(tplot,eta_sim_fcn(tplot),'b')
        t0_min = t0_dict[model,gaugeno]
        pylab.plot(tplot,eta_sim_fcn(tplot-t0_min),color[i],linewidth=2)

    pylab.plot(t_dart,eta_dart,'k',linewidth=2)
    legendstr = [modelnum[m] for m in models] + ['DART']
    pylab.legend(legendstr)
    pylab.title("DART %s " % gaugeno)
    pylab.xlim([t1,t2])
    pylab.draw()

def plot_all_multidart(t0_dict):
    #pylab.figure(13,(8,8))
    #pylab.clf()
    jsp = 0
    for gaugeno in gaugenos:
        if gaugeno==21401: t1 = 3000; t2 = t1+2.*3600;
        if gaugeno==21413: t1 = 4000; t2 = t1+2.*3600;
        if gaugeno==21418: t1 = 1400; t2 = t1+2.*3600;
        if gaugeno==21419: t1 = 4000; t2 = t1+2.*3600;

        jsp = jsp+1
        #pylab.subplot(4,3,jsp)
        models = ['GCMT', 'Hayes', 'UCSB3']
        plot_dart_multimodel(gaugeno,models,t0_dict,t1,t2)
        fname = 'figs/dart%s_1-3.%s' % (gaugeno,gext)
        pylab.savefig(fname)
        print "Created ",fname

        jsp = jsp+1
        #pylab.subplot(4,3,jsp)
        models = "Ammon Caltech Fujii".split()
        plot_dart_multimodel(gaugeno,models,t0_dict,t1,t2)
        fname = 'figs/dart%s_4-6.%s' % (gaugeno,gext)
        pylab.savefig(fname)
        print "Created ",fname

        jsp = jsp+1
        #pylab.subplot(4,3,jsp)
        models = "Saito Gusman GusmanAU pmelWei".split()
        plot_dart_multimodel(gaugeno,models,t0_dict,t1,t2)
        fname = 'figs/dart%s_7-9.%s' % (gaugeno,gext)
        pylab.savefig(fname)
        print "Created ",fname

    #pylab.savefig("all_darts.%s" % gext)

def make_all():
    t0_dict,rms_shift_dict,rms_noshift_dict, rms_nosim_dict = run_all()
    plot_all(t0_dict,rms_shift_dict,rms_noshift_dict,rms_nosim_dict)
    plot_all_multidart(t0_dict)

def check_results():
    t_dart,eta_dart,eta_dart_fcn, t_sim,eta_sim,eta_sim_fcn \
      = get_data(21418,'UCSB3')
    t1 = 1400; t2 = t1+2.*3600;
    t = np.arange(t1-500,t2+500,15)
    t0_min, rms_min = \
        minimize_rms(t,eta_dart_fcn,eta_sim_fcn,0.,t1,t2,pTd)
    print "UCSB3:  t0 = ",t0_min, "  RMS = ",rms_min
    pylab.figure(30)
    pylab.clf()
    pylab.plot(t_dart,eta_dart,'k')
    pylab.plot(t,eta_sim_fcn(t-t0_min),'b')
    err_UCSB3 = eta_dart_fcn(t) - eta_sim_fcn(t-t0_min)

    t_dart,eta_dart,eta_dart_fcn, t_sim,eta_sim,eta_sim_fcn \
        = get_data(21418,'GCMT')
    t0_min, rms_min = \
        minimize_rms(t,eta_dart_fcn,eta_sim_fcn,0.,t1,t2,pTd)
    print "GCMT:  t0 = ",t0_min, "  RMS = ",rms_min
    pylab.figure(30)
    pylab.plot(t,eta_sim_fcn(t-t0_min),'r')
    err_CGMT = eta_dart_fcn(t) - eta_sim_fcn(t-t0_min)

    return t, err_UCSB3, err_CGMT

