"""
Module for computing seafloor deformation using Okada model.

Okada model is a mapping from several fault parameters
to a surface deformation.
See Okada 1985, or Okada 1992, Bull. Seism. Soc. Am.

some routines adapted from fortran routines written by
Xiaoming Wang.

"""
#from numpy import *
from pylab import *
import os
import string
import fixdata

poisson = 0.25  # Poisson ratio

#=================================================================================
def builddeffile (okadaparamfile,faultparamfile,outfile):

    faultparams=getokadaparams(okadaparamfile)
    faultparams.update(getfaultparams(faultparamfile))

    fid=open(outfile,'w')

    X=linspace(faultparams['xlower'],faultparams['xupper'],faultparams['mx'])
    Y=linspace(faultparams['ylower'],faultparams['yupper'],faultparams['my'])

    dZ=okadamap(faultparams,X,Y)
    ind=fixdata.findbadindices(dZ)
    if ind:
        dZ=fixdata.fillbaddata(dZ,ind)

    dZ = filtermask(dZ,faultparams)
    #pdb.set_trace()
    for jj in xrange(faultparams['my']):
        j=-1-jj
        for i in xrange(faultparams['mx']) :
            fid.write('%012.6e %012.6e %012.6e \n' % (X[i],Y[j],dZ[j,i]))

    fid.close()
    return

#=================================================================================
def builddynamicdeffile (okadaparamfile,faultparamfile,outfile,t0=0.0, tend=1.0, nt = 2):

    faultparams=getokadaparams(okadaparamfile)
    faultparams.update(getfaultparams(faultparamfile))

    fid=open(outfile,'w')

    X=linspace(faultparams['xlower'],faultparams['xupper'],faultparams['mx'])
    Y=linspace(faultparams['ylower'],faultparams['yupper'],faultparams['my'])

    T=linspace(t0,tend,nt)

    dZ=okadamap(faultparams,X,Y)
    ind=fixdata.findbadindices(dZ)
    if ind:
        dZ=fixdata.fillbaddata(dZ,ind)

    dZ = filtermask(dZ,faultparams)
    #pdb.set_trace()
    for it in T:
        alpha=(it-t0)/(tend-t0)
        for jj in xrange(faultparams['my']):
            j=-1-jj
            for i in xrange(faultparams['mx']) :
                fid.write('%012.6e %012.6e %012.6e %012.6e \n' % (it,X[i],Y[j],alpha*dZ[j,i]))

    fid.close()
    return

#=================================================================================
def getokadaparams (infile):

    """
    obtain parameters necessary for okada map from config file: infile

        file format:
        parameter names and values should appear on the same single line seperated by a space
    """

    keylist=["Focal_Depth","Fault_Length","Fault_Width","Dislocation","Strike_Direction", \
             "Dip_Angle","Slip_Angle","Fault_Latitude","Fault_Longitude"]

    okadaparams={}
    fid=open(infile,'r')
    keyleft=len(keylist)
    while keyleft> 0 :
        line=string.split(fid.readline())
        if line:
            if line[0] in keylist:
                okadaparams[line[0]]=float(line[1])
                keyleft=keyleft-1
            if line[1] in keylist:
                okadaparams[line[1]]=float(line[0])
                keyleft=keyleft-1

    for key in keylist :
        if not key in okadaparams:
            print('ERROR: parameters for okada fault not fully specified in %s' % (infile))
            exit

    fid.close()
    return okadaparams
    #end getokadaparams=============================================================

#===================================================================================
def getfaultparams (infile):

    """
    obtain params from a file that specify a fault grid from infile
    params are xlower,ylower,dx,dy,mx,my, OR
    xlower,ylower,xupper,yupper,mx,my

    file format:
        parameter names and values should appear on the same single line seperated by a space
    """

    keylist=["xlower","ylower","xupper","yupper","dx","dy","mx","my"]

    faultgridparams={}
    fid=open(infile,'r')
    keyleft=len(keylist)-2
    while keyleft> 0 :
        line=string.split(fid.readline())
        if line:
            if line[0] in keylist:
                faultgridparams[line[0]]=float(line[1])
                keyleft=keyleft-1
            if line[1] in keylist:
                faultgridparams[line[1]]=float(line[0])
                keyleft=keyleft-1

    faultgridparams['mx'] = int(faultgridparams['mx'])
    faultgridparams['my'] = int(faultgridparams['my'])

    if faultgridparams.has_key('dx')& faultgridparams.has_key('dy'):
        faultgridparams['xupper'] = faultgridparams['xlower'] + faultgridparams['dx']*(faultgridparams['mx']-1)
        faultgridparams['yupper'] = faultgridparams['ylower'] + faultgridparams['dy']*(faultgridparams['my']-1)
    elif faultgridparams.has_key('xupper')&faultgridparams.has_key('yupper'):
        faultgridparams['dx'] = (faultgridparams['xupper']-faultgridparams['xlower'])/(faultgridparams['mx']-1)
        faultgridparams['dy'] = (faultgridparams['yupper']-faultgridparams['ylower'])/(faultgridparams['my']-1)
    else:
        print('ERROR: parameters for fault grid not fully specified in %s' % (infile))
        exit

    for key in keylist :
        if not key in faultgridparams:
            print('ERROR: parameters for fault grid not fully specified in %s' % (infile))
            exit

    fid.close()
    return faultgridparams
    #end getfaultparams===========================================================


#=============================================================================
def  okadamap(okadaparams,X,Y):

    """
    create displacement matrix dZ for a surface displacement
    over gridded region defined by X,Y, vectors of length nx,ny
    given okadaparams
    """

    rad = pi/180.       # conversion factor from degrees to radians
    rr = 6.378e6        # radius of earth
    lat2meter = rr*rad  # conversion factor from degrees latitude to meters

    hh =  okadaparams["Focal_Depth"]
    L  =  okadaparams["Fault_Length"]
    w  =  okadaparams["Fault_Width"]
    d  =  okadaparams["Dislocation"]
    th =  okadaparams["Strike_Direction"]
    dl =  okadaparams["Dip_Angle"]
    rd =  okadaparams["Slip_Angle"]
    y0 =  okadaparams["Fault_Latitude"]
    x0 =  okadaparams["Fault_Longitude"]
    location =  okadaparams.get("LatLong_Location", "top center")

    ang_dip = rad*dl
    ang_slip = rad*rd
    ang_strike = rad*th
    halfL = 0.5*L

    plot_plane = True
    print_xy = False

    if plot_plane:
        figure(2)
        #clf()
        plot([x0],[y0],'bo')
    if print_xy:
        print "x0,y0: ",x0,y0


    if location == "top center":

        # Convert focal depth used for Okada's model
        # from top of fault plane to bottom:
        hh = hh + w*sin(ang_dip)

        # Convert fault origin from top of fault plane to bottom:
        del_x = w*cos(ang_dip)*cos(ang_strike) / (lat2meter*cos(y0*rad))
        del_y = w*cos(ang_dip)*sin(ang_strike) / lat2meter

    elif location == "centroid":

        # Convert focal depth used for Okada's model
        # from middle of fault plane to bottom:
        hh = hh + 0.5*w*sin(ang_dip)

        # Convert fault origin from middle of fault plane to bottom:
        del_x = 0.5*w*cos(ang_dip)*cos(ang_strike) / (lat2meter*cos(y0*rad))
        del_y = 0.5*w*cos(ang_dip)*sin(ang_strike) / lat2meter

    else:
        raise ValueError("Unrecognized LatLong_Location: %s" % location)

    x0 = x0 + del_x
    y0 = y0 - del_y

    if print_xy:
        print "del_x, del_y: ",del_x,del_y
        print "new x0,y0: ",x0,y0
    if plot_plane:
        plot([x0],[y0],'ro')
        #legend(['Top','Bottom'],'upper left')
        halfLc = halfL*sin(ang_strike) / (lat2meter*cos(y0*rad))
        halfLs = halfL*cos(ang_strike) / lat2meter
        plot([x0+halfLc, x0-del_x+halfLc, x0-del_x-halfLc,x0-halfLc,x0+halfLc], \
              [y0+halfLs, y0+del_y+halfLs, y0+del_y-halfLs, y0-halfLs, y0+halfLs], 'r-')
        axis([X[0],X[-1],Y[0],Y[-1]])


    x,y = meshgrid(X,Y)

    # Convert distance from (x,y) to (x0,y0) from degrees to meters:
    xx = lat2meter*cos(rad*y)*(x-x0)   
    yy = lat2meter*(y-y0)


    # Convert to distance along strike (x1) and dip (x2):
    x1 = xx*sin(ang_strike) + yy*cos(ang_strike) 
    x2 = xx*cos(ang_strike) - yy*sin(ang_strike) 

    # In Okada's paper, x2 is distance up the fault plane, not down dip:
    x2 = -x2

    if 0:
        figure(3)
        clf()
        plot([xx[0,0],xx[0,-1],xx[-1,-1],xx[-1,0],xx[0,0]], \
             [yy[0,0],yy[0,-1],yy[-1,-1],yy[-1,0],yy[0,0]], 'k-')
        
        plot([x1[0,0],x1[0,-1],x1[-1,-1],x1[-1,0],x1[0,0]], \
             [x2[0,0],x2[0,-1],x2[-1,-1],x2[-1,0],x2[0,0]], 'b-')
        
    p = x2*cos(ang_dip) + hh*sin(ang_dip)
    q = x2*sin(ang_dip) - hh*cos(ang_dip)

    f1=strike_slip (x1+halfL,p,  ang_dip,q)
    f2=strike_slip (x1+halfL,p-w,ang_dip,q)
    f3=strike_slip (x1-halfL,p,  ang_dip,q)
    f4=strike_slip (x1-halfL,p-w,ang_dip,q)

    g1=dip_slip (x1+halfL,p,  ang_dip,q)
    g2=dip_slip (x1+halfL,p-w,ang_dip,q)
    g3=dip_slip (x1-halfL,p,  ang_dip,q)
    g4=dip_slip (x1-halfL,p-w,ang_dip,q)

    # Displacement in direction of strike and dip:
    ds = d*cos(ang_slip)
    dd = d*sin(ang_slip)

    us = (f1-f2-f3+f4)*ds
    ud = (g1-g2-g3+g4)*dd

    dZ = (us+ud)

    if 0:
        contour(x,y,dZ,linspace(-8,8,17),colors='k')

    return dZ
    #=========================================================================


#==============================================================================
def strike_slip (y1,y2,ang_dip,q):
    """
    !.....Used for Okada's model
    !.. ..Methods from Yoshimitsu Okada (1985)
    !-----------------------------------------------------------------------
    """
    sn = sin(ang_dip)
    cs = cos(ang_dip)
    d_bar = y2*sn - q*cs
    r = sqrt(y1**2 + y2**2 + q**2)
    xx = sqrt(y1**2 + q**2)
    a4 = 2.0*poisson/cs*(log(r+d_bar) - sn*log(r+y2))
    f = -(d_bar*q/r/(r+y2) + q*sn/(r+y2) + a4*sn)/(2.0*3.14159)

    return f
    #==========================================================================


#==============================================================================
def dip_slip (y1,y2,ang_dip,q):
    """
    !.....Based on Okada's paper (1985)
    !.....Added by Xiaoming Wang
    !-----------------------------------------------------------------------
    """
    sn = sin(ang_dip)
    cs = cos(ang_dip)

    d_bar = y2*sn - q*cs;
    r = sqrt(y1**2 + y2**2 + q**2)
    xx = sqrt(y1**2 + q**2)
    a5 = 4.*poisson/cs*arctan((y2*(xx+q*cs)+xx*(r+xx)*sn)/y1/(r+xx)/cs)
    f = -(d_bar*q/r/(r+y1) + sn*arctan(y1*y2/q/r) - a5*sn*cs)/(2.0*3.14159)

    return f
    #===========================================================================


#================================================================================
def filtermask (dZ,faultparams):
    """
    borrowed from code written by Xiaoming Wang and Tom Logan at ARSC

    !.....Filter the deformation using a circular mask centered
    !.....at the epicenter using a calculated radius
    !.....Removes small numerical artifacts away from the epicenter
    """
    filterindices=[]

    osixty = 0.016666666667
    rad = 0.01745329252
    rr = 6.378e6

    xo = faultparams['xlower']
    yo = faultparams['ylower']
    nx = faultparams['mx']
    ny = faultparams['my']
    spacing = faultparams['dx']

    x0 = faultparams['Fault_Longitude']
    y0 = faultparams['Fault_Latitude']
    l =  faultparams['Fault_Length']
    w =  faultparams['Fault_Width']
    dl = faultparams['Dip_Angle']


    ang_dip = rad*dl # convert degree to radian

    #!-- fault origin in pixels -----------
    ypix = (y0-yo)/spacing
    xpix = (x0-xo)/spacing

    #!-- conversion from meters to pixels ---
    tmpd=spacing*rad
    xdist = tmpd*rr

    #!-- size of the fault in pixels --------
    npix_x = l/xdist
    npix_y = w/xdist

    #!-- set the range (radius) of the filter circle --------
    #!----- for small dip angles, use the length and width --
    #!----- for larger dip angles, use only the length ------

    if dl<30.0:
        drange = 1.5 * cos(ang_dip)*sqrt(npix_x*npix_x+npix_y*npix_y)
    else:
        drange = 1.2 * npix_x

    print("Filtering deformation using a circle of radius %s" % (drange))

    #!-- Create the filtering mask ----------
    for i in xrange(nx):
        for j in xrange(ny) :
            dist = sqrt((i+1-xpix)**2+(j+1-ypix)**2)
            if dist > drange :
                filterindices.append((j,i))

    #!-- apply the filter to the actual deformation ------
    dZ = filterdata(dZ,filterindices,radius=2)

    return dZ


#============================================================================

def filterdata (Z,filterinds,radius=1):
    """
    filter data in array z, at indice tuples in list filterinds
    by averaging surrounding data, ball with radius=radius in inf-norm
    acts as a low-band pass filter and removes oscillatory data
    """

    m=shape(Z)[0]
    n=shape(Z)[1]

    for ind in filterinds :
        i=ind[0]
        j=ind[1]
        r=radius

        irange=range(max(0,i-r),min(i+r+1,m))
        jrange=range(max(0,j-r),min(j+r+1,n))
        summands=0
        sum=0.
        for ii in irange:
            for jj in jrange:
                ballind=(ii,jj)
                sum = sum + Z[ballind[0],ballind[1]]
                summands=summands+1
        if summands >0 : 
            Z[ind[0],ind[1]] = sum/summands

    return Z


def make_dz_chile2010():
    dtopo_cfg = 'testdata/usgs100227-new.cfg'

    #dtopo_fname = 'testdata/usgs100227.tt1'
    #print "Using Okada model to create %s " % dtopo_fname

    faultparams=getokadaparams(dtopo_cfg)
    faultparams.update(getfaultparams(dtopo_cfg))
    faultparams["LatLong_Location"] = "top center"
    print "Fault parameters: \n", faultparams

    X=linspace(faultparams['xlower'],faultparams['xupper'],faultparams['mx'])
    Y=linspace(faultparams['ylower'],faultparams['yupper'],faultparams['my'])

    dZ=okadamap(faultparams,X,Y)
    return X,Y,dZ






