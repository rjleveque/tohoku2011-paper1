
"""
Process fixed grid output to determine runup at each latitude for run along
the Japan coast.

Use frame 5 (before waves hit but AMR grids set) to determine shoreline as
cells that are wet but have a dry neighbor or vice versa.  Then use etamax
in frame 75 to determine eta in these cells so that these cells can be
plotted vs. latitude (in plotdata.py).
"""

from pylab import *
from pyclaw.plotters import plotfg
from pyclaw.plotters.colormaps import white_blue
from numpy import ma

def readraw():
    xraw,yraw,etaraw = loadtxt('rawdata.txt',unpack=True)
    return yraw,etaraw

def readh0():
    fg = plotfg.ClawPlotFGData()
    fg.outdir = '.'
    grid0,soln0 = fg.get_frame(5)
    h0 = soln0.h
    x = grid0.xcenter
    y = grid0.ycenter
    mx = grid0.mx
    my = grid0.my
    return x,y,mx,my,h0

def process_eta(dirname,x,y,mx,my,h0):
    fg = plotfg.ClawPlotFGData()
    fg.outdir = dirname
    grid1,soln1 = fg.get_frame(75)
    etamax1 = soln1.fg[:,:,6]

    mask0 = (h0 <= 0.)

    mask = zeros(h0.shape)
    for j in range(1,my-1):
        for i in range(1,mx-1):
            if (mask0[j,i] != mask0[j+1,i]) or (mask0[j,i] != mask0[j-1,i]) \
                   or (mask0[j,i] != mask0[j,i+1]) or (mask0[j,i] != mask0[j,i-1]):
                mask[j,i] = 1

    etamax1mask = etamax1 * mask

    if 0:
        figure(1)
        clf()
        pcolormesh(x,y,etamax1mask,cmap=white_blue)
        contour(x,y,mask,[0.5,0.5],colors='r')


    etaflat = reshape(etamax1mask,(mx*my,1))
    X,Y = meshgrid(x,y)
    Xflat = reshape(X,(mx*my,1))
    Yflat = reshape(Y,(mx*my,1))
    i = (etaflat>0.1)
    etaflat = etaflat[i]
    Xflat = Xflat[i]
    Yflat = Yflat[i]
    return Xflat, Yflat, etaflat

def plot_runup(yraw,etaraw,Yflat,etaflat):
    figure(2)
    clf()
    plot(yraw,etaraw,'k.')
    plot(Yflat,etaflat,'r-',linewidth=2)
    xlim(36,41)
    ylim(0,40)



def process_all():
    models = """GCMT Hayes UCSB3 Ammon Caltech Fujii Saito 
                Gusman GusmanAU PMELWei""".split()
    x,y,mx,my,h0 = readh0()
    for model in models:
        dirname = 'fgoutput/' + model
        xflat,yflat,etaflat = process_eta(dirname,x,y,mx,my,h0)
        d = vstack((xflat,yflat,etaflat)).T
        fname = '%s.txt' % model
        savetxt(fname,d)
        print "Created ",fname
