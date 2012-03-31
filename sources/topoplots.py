
import colormaps
from matplotlib.colors import Normalize 
import topotools
from numpy import ma

class TopoPlotData(object):
    def __init__(self, fname=None, topotype=3):
        self.fname = fname 
        self.x = None
        self.y = None 
        self.X = None
        self.Y = None
        self.Z = None 
        self.topotype = topotype
        self.neg_cmap = None
        self.pos_cmap = None
        self.cmap = None
        self.cmax = 100.
        self.cmin = -4000.
        self.climits = None
        self.figno = 200
        self.clf = False
        self.addcolorbar = False
        self.addcontour = False
        self.contour_levels = [0, 0]
        self.xlimits = None
        self.ylimits = None
        self.coarsen = 1
        self.imshow = True
        self.gridedges_show = True
        self.print_fname = True

    def read(self):
        read_topo(self)
        
    def plot(self):
        plot_topo(self)
        

    
def read_topo(topoplotdata):
    """
    Read in a topo or bathy file from topoplotdata.fname.
    """
    from topotools import topofile2griddata 
    if topoplotdata.fname is None:
        raise Exception("Need to specify fname attribute")
    try:
        X,Y,Z,x,y = topofile2griddata(topoplotdata.fname, topotype=topoplotdata.topotype)
    except:
        raise Exception("Problem reading data")
    topoplotdata.x = x
    topoplotdata.y = y
    topoplotdata.X = X
    topoplotdata.Y = Y
    topoplotdata.Z = Z
    print "Read arrays of shape %s from %s" % (topoplotdata.Z.shape, topoplotdata.fname)
    print "     and set x,y,z attributes of input object"
    #return topoplotdata
        
        
def plot_topo(topoplotdata):
    """
    Plot topo data as specified by topoplotdata.
    """

    import os
    import pylab

    if topoplotdata.Z is None:
        topoplodata = read_topo(topoplotdata)
        
    x = topoplotdata.x
    y = topoplotdata.y
    topo = topoplotdata.Z
    X = topoplotdata.X
    Y = topoplotdata.Y
    xllcorner = x[0]
    yllcorner = y[0]
    cellsize = x[1]-x[0]
    
    
    fname = topoplotdata.fname 
    
    topotype = topoplotdata.topotype
    if topoplotdata.climits:
        # deprecated option
        cmin = topoplotdata.climits[0]
        cmax = topoplotdata.climits[1]
    else:
        cmin = topoplotdata.cmin
        cmax = topoplotdata.cmax
    figno = topoplotdata.figno
    addcolorbar = topoplotdata.addcolorbar
    addcontour = topoplotdata.addcontour
    contour_levels = topoplotdata.contour_levels
    xlimits = topoplotdata.xlimits
    ylimits = topoplotdata.ylimits
    coarsen = topoplotdata.coarsen
    imshow = topoplotdata.imshow
    gridedges_show = topoplotdata.gridedges_show
    neg_cmap = topoplotdata.neg_cmap
    pos_cmap = topoplotdata.pos_cmap
    cmap = topoplotdata.cmap
    clf = topoplotdata.clf
    print_fname = topoplotdata.print_fname



    if neg_cmap is None:
        neg_cmap = colormaps.make_colormap({cmin:[0.3,0.2,0.1],
                                                -0.00001:[0.95,0.9,0.7],
                                                0.00001:[.5,.7,0]})
    if pos_cmap is None:
        pos_cmap = colormaps.make_colormap({ -0.00001:[0.95,0.9,0.7],
                                     0.00001:[.5,.7,0],
                                        cmax:[.2,.5,.2]})
    if cmap is None:
        cmap = colormaps.make_colormap({cmin:[0.3,0.2,0.1],
                                           -0.00001:[0.95,0.9,0.7],
                                           0.00001:[.5,.7,0],
                                           cmax:[.2,.5,.2]})
        #cmap = colormaps.make_colormap({-1:[0,0,1],0:[1,1,1],1:[1,0,0]})
        
    #-------------
    if 0:
        if abs(topotype) == 1:

            X,Y,topo = topotools.topofile2griddata(fname, topotype)
            topo = pylab.flipud(topo)
            Y = pylab.flipud(Y)
            x = X[0,:]
            y = Y[:,0]
            xllcorner = x[0]
            yllcorner = y[0]
            cellsize = x[1]-x[0]


        elif abs(topotype) == 3:

            file = open(fname, 'r')
            lines = file.readlines()
            ncols = int(lines[0].split()[0])
            nrows = int(lines[1].split()[0])
            xllcorner = float(lines[2].split()[0])
            yllcorner = float(lines[3].split()[0])
            cellsize = float(lines[4].split()[0])
            NODATA_value = int(lines[5].split()[0])

            print "Loading file ",fname
            print "   nrows = %i, ncols = %i" % (nrows,ncols)
            topo = pylab.loadtxt(fname,skiprows=6,dtype=float)
            print "   Done loading"

            if 0:
                topo = []
                for i in range(nrows):
                    topo.append(pylab.array(lines[6+i],))
                print '+++ topo = ',topo
                topo = pylab.array(topo)

            topo = pylab.flipud(topo)

            x = pylab.linspace(xllcorner, xllcorner+ncols*cellsize, ncols)
            y = pylab.linspace(yllcorner, yllcorner+nrows*cellsize, nrows)
            print "Shape of x, y, topo: ", x.shape, y.shape, topo.shape

        else:
            raise Exception("*** Only topotypes 1 and 3 supported so far")
        
        if topotype < 0:
            topo = -topo
        
    
    #------------------

    ncols,nrows = topo.shape
    #import pdb; pdb.set_trace()
    
    if coarsen > 1:
        topo = topo[slice(0,nrows,coarsen), slice(0,ncols,coarsen)]
        X = X[slice(0,nrows,coarsen), slice(0,ncols,coarsen)]
        Y = Y[slice(0,nrows,coarsen), slice(0,ncols,coarsen)]
        x = x[slice(0,ncols,coarsen)]
        y = y[slice(0,nrows,coarsen)]
        print "Shapes after coarsening: ", x.shape, y.shape, topo.shape


    if figno:
        pylab.figure(figno)
    if clf:
        pylab.clf()

    if topoplotdata.imshow:
            color_norm = Normalize(cmin,cmax,clip=True)
            xylimits = (x[0],x[-1],y[0],y[-1])
            #pylab.imshow(pylab.flipud(topo.T), extent=xylimits, \
            pylab.imshow(topo, extent=xylimits, \
                    cmap=cmap, interpolation='nearest', \
                    norm=color_norm)
            #pylab.clim([cmin,cmax])
            if addcolorbar:
                pylab.colorbar(shrink=0.5)
    else:
        neg_topo = ma.masked_where(topo>0., topo)
        all_masked = (ma.count(neg_topo) == 0)
        if not all_masked:
            pylab.pcolormesh(X,Y,neg_topo,cmap=neg_cmap)
            pylab.clim([cmin,0])
            if addcolorbar:
                pylab.colorbar(shrink=0.5)

        pos_topo = ma.masked_where(topo<0., topo)
        all_masked = (ma.count(pos_topo) == 0)
        if not all_masked:
            pylab.pcolormesh(X,Y,pos_topo,cmap=pos_cmap)
            pylab.clim([0,cmax])
        if addcolorbar:
            pylab.colorbar(shrink=0.5)

    pylab.axis('scaled')


    if addcontour:
        pylab.contour(X,Y,topo,levels=contour_levels,colors='k')

    patchedges_show = True
    if patchedges_show:
        pylab.plot([x[0],x[-1]],[y[0],y[0]],'k')
        pylab.plot([x[0],x[-1]],[y[-1],y[-1]],'k')
        pylab.plot([x[0],x[0]],[y[0],y[-1]],'k')
        pylab.plot([x[-1],x[-1]],[y[0],y[-1]],'k')

    if print_fname:
        fname2 = os.path.splitext(fname)[0]
        pylab.text(xllcorner+cellsize, yllcorner+cellsize, fname2, color='m')

    topodata = object()
    topodata.x = x
    topodata.y = y
    topodata.topo = topo

    return topodata

def test():
    tpd = TopoPlotData('CrescentCity.asc', topotype=-3)
    tpd.read()
    tpd.coarsen = 1
    tpd.addcolorbar = True
    tpd.cmin = -30
    tpd.cmax = 20
    tpd.clf = True
    tpd.plot()
    
    # pcolor version:
    tpd.imshow = False
    tpd.figno = 201
    tpd.plot()
    return tpd
