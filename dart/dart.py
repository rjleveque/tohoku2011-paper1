from pylab import *
import numpy as np
import datetime

def readdart(fname, t1, t2, t_quake=None):
    """
    Read DART data from file fname.
    t1, t2 = upper and lower bounds on times to return t and eta.
    t_quake = time of earthquake, a datetime object, if known.

    Returns:
    t = array of datetime objects for time series,
    eta = surface elevation at each time,
    t_sec = times in seconds since t_quake (or since t1 if t_quake==None).
    """

    d = loadtxt(fname,skiprows=2);
    t = [datetime.datetime(int(dd[0]),int(dd[1]),int(dd[2]),int(dd[3]),int(dd[4]))\
        for dd in d]
    t = np.array(t)
    i = (t >= t1) & (t <= t2)
    t = t[i]
    eta = d[i,7]
    mask = eta < 9990
    eta = eta[mask]
    t = t[mask]
    if t_quake is not None:
        t0 = t_quake
    else:
        t0 = t1

    # Compute seconds past t0 for each element of t:
    dt = t - t0  # an array of datetime.timedelta objects
    t_sec = np.array([dtj.total_seconds() for dtj in dt])

    return t, t_sec, eta

def plot_datetime(t, eta, fig=None, ax=None):
    """
    Plot depth eta vs. time t where t is an array of datetime objects.
    fig is a matplotlib.figure.Figure object or an integer.
    ax is a matplotlib.axis object.
    Improves the formatting of tick labels to make them readable.
    """
    from matplotlib.dates import DateFormatter
    formatter = DateFormatter('%d %h %H:%M')
    if fig==None:
        fig = figure()
        clf()
    elif type(fig) is int:
        fig = figure(fig)
        clf()
    if ax==None:
        clf()
        ax = axes()
    ax.xaxis.set_major_formatter(formatter)
    xticks(rotation=25)

    class DepthFormatter(Formatter):
        def __call__(self,x,pos=None):
            return "%8.2f" % x
    ax.yaxis.set_major_formatter(DepthFormatter())

    fig.autofmt_xdate()
    plot(t,eta)
    return fig,ax




def plotdart(fname, t_quake):
    t1 = t_quake - datetime.timedelta(0, 12*3600)
    t2 = t_quake + datetime.timedelta(0, 36*3600)
    t, t_sec, eta = readdart(fname, t1, t2, t_quake)
    plot_datetime(t, eta)
    return t,t_sec,eta
    

def fit_tide_poly(t, eta,degree, t1fit,t2fit, t1out,t2out):
    """
    Fit a polynomial of the specified degree to data in the range t1fit <= t <= t2fit.
    Returns the coefficents c of c[0] + c[1]*t + ...
    and detided data eta_notide over the range t1out <= t <= t2out.
    """
    from numpy.linalg import lstsq
    
    # select subset of t, eta where fit is done:
    mask = ((t>=t1fit) & (t<=t2fit)) 
    tfit = t[mask]
    etafit = eta[mask]
    
    
    # select subset of t, eta for output:
    mask = ((t>=t1out) & (t<=t2out)) 
    tout = t[mask]
    etaout = eta[mask]
    
    # Scale data so matrix is well-conditioned:
    scale_factor = tfit[0]
    tfit = tfit/scale_factor
    tout = tout/scale_factor
    
    # Use Newton polynomial basis using these points:
    tpts = linspace(tfit.min(),tfit.max(),degree+1)
    
    # Form A matrix Afit for least squares fit and
    # Aout for applying fit to output data:
    Afit = ones((len(tfit),degree+1))
    Aout = ones((len(tout),degree+1))
    for j in range(1,degree+1):
        Afit[:,j] = Afit[:,j-1] * (tfit - tpts[j])
        Aout[:,j] = Aout[:,j-1] * (tout - tpts[j])
        
    # Performs least squares fit:
    c = lstsq(Afit,etafit)[0]
    
    #import pdb; pdb.set_trace()
    
    # evaluate polynomial at the output times:
    etaoutfit = dot(Aout,c)
    
    # evaluate polynomial at the fit times:
    etafit2 = dot(Afit,c)
    
    # Compute de-tided values by subtracting fit values from raw data:
    tout = tout*scale_factor
    tfit = tfit*scale_factor
    t_notide = tout
    eta_notide =  etaout - etaoutfit
    
    # plot fit and de-tided data:
    figure(70,figsize=(8,8))
    clf()
    subplot(211)
    plot(tfit,etafit,'b')
    plot(tout,etaout,'g')
    plot(tfit,etafit2,'k')
    plot(tout,etaoutfit,'r')
    legend(['raw data over [t1fit, t2fit]', 'raw data over [t1out, t2out]', \
            'fit to data over [t1fit, t2fit]','fit over [t1out, t2out]'], \
             loc=0)
    ymin = etafit.min() - 0.5*(etafit.max()-etafit.min())
    ymax = etafit.max() + 0.5*(etafit.max()-etafit.min())
    ylim([ymin,ymax])
    subplot(212)
    plot(t_notide, eta_notide,'k')
    title('de-tided data over [t1out, t2out]')
        
    return c, t_notide, eta_notide

def plot_post_quake(t,eta,gaugeno=''):
    thours = t/3600.
    figure(63)
    clf()
    plot(thours,eta)
    xlabel("Hours after quake")
    title("DART #%s" % gaugeno)
    
