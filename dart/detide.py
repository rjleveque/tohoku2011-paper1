"""
Detide <gaugeno> using data from
    http://www.ndbc.noaa.gov/station_history.php?station=<gaugeno>
Run this from within python shell to see the plots.
"""

import datetime
import sys, os
import pylab
import numpy as np
import dart  

# old value from March 2011 runs, might have been centroid time but seems
# too long after quake start:
#t_quake = datetime.datetime(2011, 3, 11, 5, 48, 15) 
          
# initial time:
t_quake = datetime.datetime(2011, 3, 11, 5, 46, 24) 

def detide(gaugeno):

    fname = '%s.txt' % gaugeno
    if gaugeno == 21401:
        t1fit = -12.*3600.
        t2fit = 36.*3600.
        t1out = 0.
        t2out = 12.*3600.
        degree = 15
    else:
        t1fit = -12.*3600.
        t2fit = 36.*3600.
        t1out = 0.
        t2out = 12.*3600.
        degree = 15

    fname_notide = os.path.splitext(fname)[0] + '_notide.txt'
    t,t_sec,eta = dart.plotdart(fname, t_quake)

    c,t_notide,eta_notide = dart.fit_tide_poly(t_sec,eta,degree,\
                               t1fit,t2fit, t1out,t2out)


    # fix bad data value:
    if gaugeno==21418:
        print "Bad data: ",eta_notide[121:124]
        eta_notide[122] = -0.075
        print "Fix data: ",eta_notide[121:124]
        print "Bad data: ",eta_notide[152:155]
        eta_notide[153] = -0.07
        print "Fix data: ",eta_notide[152:155]

    d = np.vstack([t_notide,eta_notide]).T
    pylab.savetxt(fname_notide,d)
    print "Created file ",fname_notide

    dart.plot_post_quake(t_notide,eta_notide,gaugeno)


    pylab.figure(70)
    pylab.subplot(211)
    pylab.title('DART %s' % gaugeno)

    if 0:
        pylab.figure(75)
        pylab.plot(t_notide,eta_notide,'k')
        pylab.xlim(0,3.*3600)
        pylab.ylim(-1,2)

    return t_notide, eta_notide

def make_all():
    for gaugeno in [21401,21413,21418,21419]:
        detide(gaugeno)

if __name__=="__main__":
    make_all()

