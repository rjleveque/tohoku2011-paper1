
from pylab import *
from okada2 import okadamap
import colormaps
import html_movie


import re,sys

Rearth = 6367.5e3
class FaultModel(object):
    pass

def read_subfault_model(fname):
    lines = open(fname,'r').readlines()
    regexp_dx = re.compile(r"Dx=[ ]*(?P<dx>[^k]*)")
    regexp_dy = re.compile(r"Dy=[ ]*(?P<dy>[^k]*)")
    regexp_nx = re.compile(r"nx[^=]*=[ ]*(?P<nx>[^D]*)")
    regexp_ny = re.compile(r"ny[^=]*=[ ]*(?P<ny>[^D]*)")
    
    i_dxdy = -1
    for i,line in enumerate(lines):
        result_dx = regexp_dx.search(line)
        result_dy = regexp_dy.search(line)
        result_nx = regexp_nx.search(line)
        result_ny = regexp_ny.search(line)
        if result_dx and result_dy:
            i_dxdy = i
            dx = float(result_dx.group('dx')) * 1.e3
            dy = float(result_dy.group('dy')) * 1.e3
            nx = int(result_nx.group('nx'))
            ny = int(result_ny.group('ny'))
            print "Found dx = %s, dy = %s, nx = %s, ny = %s in line %s" \
                  % (dx, dy, nx, ny, i_dxdy)
            
    if i_dxdy == -1:
        print "*** Did not find a line containing both Dx and Dy"
        raise Exception()
    
    fault_bdry = []
    for i in range(i_dxdy+3, i_dxdy+8):    
        fault_bdry.append([float(s) for s in lines[i].split()])
    fault_bdry = array(fault_bdry)
    #print "+++ fault_bdry = ", fault_bdry

        
    
    data = loadtxt(fname,skiprows=9+i_dxdy)
    arrayshape = (ny,nx)
    fm = FaultModel()
    fm.arrayshape = arrayshape
    fm.subfault_length = dx
    fm.subfault_width = dy
    
    fm.latitude = reshape(data[:,0],arrayshape)
    fm.longitude = reshape(data[:,1],arrayshape)
    fm.depth = reshape(data[:,2],arrayshape) * 1.e3 # convert km to meters
    fm.slip = reshape(data[:,3],arrayshape) * 0.01  # convert cm to meters
    fm.rake = reshape(data[:,4],arrayshape)
    fm.strike = reshape(data[:,5],arrayshape)
    fm.dip = reshape(data[:,6],arrayshape)

    # format used by UCSB:
    fm.rupture_initial_time = reshape(data[:,7],arrayshape)
    fm.rise_time_starting = reshape(data[:,8],arrayshape)
    fm.rise_time_ending = reshape(data[:,9],arrayshape)

    fm.fault_bdry = fault_bdry
    
    plot_plane = True

    if plot_plane:
        # uncomment to plot source region and depth:
        figure(2)
        clf()
        plot(fault_bdry[:,0], fault_bdry[:,1], 'k')
        plot(fm.longitude, fm.latitude, 'go')
        #pcolor(fm.longitude, fm.latitude, fm.depth)
        pcolor(fm.longitude, fm.latitude, fm.slip)
        colorbar()
        #plot_topo(topofile='etopo1min139E147E34N41N.asc',figno=100)
        title("Slip and subfault long-lat points")

    if 0:
        # Not working...
        dx0 = dx
        dy0 = dy
        for i in range(ny):
            for j in range(nx):
                x0 = fm.longitude[i,j]
                y0 = fm.latitude[i,j]
                dx = dx0 / (Rearth * cos(y0) * pi/180.)
                dy = dy0 / (Rearth * pi/180.)
                xsf = array([x0-dx/2, x0+dx/2, x0+dx/2, x0-dx/2])
                ysf = array([y0-dy/2, y0-dy/2, y0+dy/2, y0+dy/2])
                fill(xsf, ysf, 'b')
                print xsf, ysf, "\n"
    
    return fm
    
def make_okada_dz(fm, faultparams):
    mx = faultparams['mx']
    my = faultparams['my']
    
    X=linspace(faultparams['xlower'],faultparams['xupper'],mx)
    Y=linspace(faultparams['ylower'],faultparams['yupper'],my)
    dZ = zeros((my,mx))
    
    okadaparams = {}
    
    print "Making Okada dZ for each of %s subfaults" \
            % str(fm.arrayshape[0]*fm.arrayshape[1])

    for i in range(fm.arrayshape[0]):
        for j in range(fm.arrayshape[1]): 
            sys.stdout.write("%s.." % str(j+i*fm.arrayshape[1]))
            sys.stdout.flush()
            okadaparams["Focal_Depth"] = fm.depth[i,j]
            okadaparams["Fault_Length"] = fm.subfault_length
            okadaparams["Fault_Width"] = fm.subfault_width
            okadaparams["Dislocation"] = fm.slip[i,j]
            okadaparams["Strike_Direction"] = fm.strike[i,j]
            okadaparams["Dip_Angle"] = fm.dip[i,j]
            okadaparams["Slip_Angle"] = fm.rake[i,j]
            #okadaparams["Epicenter_Latitude"] = fm.latitude[i,j]
            #okadaparams["Epicenter_Longitude"] = fm.longitude[i,j]
            okadaparams["Fault_Latitude"] = fm.latitude[i,j]
            okadaparams["Fault_Longitude"] = fm.longitude[i,j]
            okadaparams["LatLong_Location"] = "centroid"     # correct
            #okadaparams["LatLong_Location"] = "top center"  # used before
            
            dZ = dZ + okadamap(okadaparams, X, Y)
            
    sys.stdout.write("\nDone\n")
    return X,Y,dZ
    
def make_okada_dz_witht(fm, faultparams, times, fname):
    mx = faultparams['mx']
    my = faultparams['my']
    
    X=linspace(faultparams['xlower'],faultparams['xupper'],mx)
    Y=linspace(faultparams['ylower'],faultparams['yupper'],my)
    mt = len(times)
    dZ = zeros((mt,my,mx))
    first_time = -ones(fm.arrayshape)

    
    okadaparams = {}
    
    dZsubfault = []

    print "Making Okada dZ for each of %s subfaults" \
            % str(fm.arrayshape[0]*fm.arrayshape[1])

    for i in range(fm.arrayshape[0]):
        dZsubi = []
        for j in range(fm.arrayshape[1]): 
            sys.stdout.write("%s.." % str(j+i*fm.arrayshape[1]))
            sys.stdout.flush()
            okadaparams["Focal_Depth"] = fm.depth[i,j]
            okadaparams["Fault_Length"] = fm.subfault_length
            okadaparams["Fault_Width"] = fm.subfault_width
            okadaparams["Dislocation"] = fm.slip[i,j]
            okadaparams["Strike_Direction"] = fm.strike[i,j]
            okadaparams["Dip_Angle"] = fm.dip[i,j]
            okadaparams["Slip_Angle"] = fm.rake[i,j]
            okadaparams["Fault_Latitude"] = fm.latitude[i,j]
            okadaparams["Fault_Longitude"] = fm.longitude[i,j]
            okadaparams["LatLong_Location"] = "centroid"
            dZij = okadamap(okadaparams, X, Y)
            dZsubi.append(dZij)
        dZsubfault.append(dZsubi)
    sys.stdout.write("\nDone\n")
            
    dZsubfault = array(dZsubfault)
    dZ = zeros((mx,my))
    tprev = times.min() - 1.
    fid = open(fname, 'w')
    pngfiles = []

    if 0. not in times:
        times = hstack(([0.], times))

    print "Making time-dependent dZ for each of %s times" \
            % str(len(times))

    for frameno,t in enumerate(times):
        ruptured = []
        # add to dZ any slip during time interval tprev to t:
        for i in range(fm.arrayshape[0]):
            for j in range(fm.arrayshape[1]): 
                if (fm.rupture_initial_time[i,j] > tprev) and \
                   (fm.rupture_initial_time[i,j] <= t):
                     dZ = dZ + dZsubfault[i,j]
                     print "Adding subfault (%s,%s) at time %s...  rupture_initial_time = %s" \
                        % (i,j,t,fm.rupture_initial_time[i,j])
                     ruptured.append((i,j))

        # write out dZ at this time:
        for jj in range(len(Y)):
            j=-1-jj
            for i in range(len(X)) :
                fid.write('%012.6e %012.6e %012.6e %012.6e \n' \
                        % (t,X[i],Y[j],dZ[j,i]))
        print "Made dZ at time %s" % t
        plot_dZ = True
        if plot_dZ:
            figure(5)
            clf()
            plot(fm.fault_bdry[:,0], fm.fault_bdry[:,1], 'k')
            cmap = colormaps.blue_white_red
            pcolor(X,Y,dZ,cmap=cmap)
            for (i,j) in ruptured:
                # plot centroids of the fault segments that ruptured 
                # during this time increment:
                plot([fm.longitude[i,j]],[fm.latitude[i,j]],'go')
            clim(-10,10)
            colorbar()
            contour(X,Y,dZ,linspace(-13,13,14),colors='k')
            title("dZ at time t = %10.2f" % t)
            draw()
            fname = "dZframe%s.png" % str(frameno).zfill(4)
            savefig(fname)
            pngfiles.append(fname)

        tprev = t
    fid.close()
    print "Created dtopo file ", fname
    if plot_dZ:
        html_movie.make_movie(pngfiles, "index.html")
    return X,Y,dZ

    
def write_dz(fname,X,Y,dZ,tend=100.):
    fid = open(fname, 'w')
    t0 = 0.
    nt = 2
    T = linspace(t0, tend, nt)
    for it in T:
        alpha=(it-t0)/(tend-t0)
        
        for jj in range(len(Y)):
            j=-1-jj
            for i in range(len(X)) :
                fid.write('%012.6e %012.6e %012.6e %012.6e \n' % (it,X[i],Y[j],alpha*dZ[j,i]))

    fid.close()
    print "Created ",fname
    
def write_dz_witht(fname,X,Y,dZ,times):
    fid = open(fname, 'w')
    for it in range(len(times)):
        for jj in range(len(Y)):
            j=-1-jj
            for i in range(len(X)) :
                fid.write('%012.6e %012.6e %012.6e %012.6e \n' \
                     % (times[it],X[i],Y[j],dZ[it,j,i]))

    fid.close()
    print "Created ",fname
    

def plot_dz(X,Y,dZ):
    figure(200)
    clf()
    bwr = colormaps.blue_white_red
    pcolor(X,Y,dZ,cmap=bwr)
    clim(-3,3)
    colorbar()
    
    
def make_dz_honshu_ucsb3():
    fname_subfaults = 'UCSB3.txt'
    fname_dtopo = 'UCSB3.tt1'
    faultparams = {}
    faultparams['mx'] = 301
    faultparams['my'] = 301
    faultparams['xlower'] = 140
    faultparams['xupper'] = 146
    faultparams['ylower'] = 35
    faultparams['yupper'] = 41
    fm = read_subfault_model(fname_subfaults)
    X,Y,dZ = make_okada_dz(fm, faultparams)
    write_dz(fname_dtopo, X,Y,dZ,tend=1.)
    return X,Y,dZ
    
def make_dz_honshu_ucsb3_t100():
    fname_subfaults = 'UCSB3.txt'
    fname_dtopo = 'UCSB3_t100.tt1'
    faultparams = {}
    faultparams['mx'] = 301
    faultparams['my'] = 301
    faultparams['xlower'] = 140
    faultparams['xupper'] = 146
    faultparams['ylower'] = 35
    faultparams['yupper'] = 41
    fm = read_subfault_model(fname_subfaults)
    X,Y,dZ = make_okada_dz(fm, faultparams)
    write_dz(fname_dtopo, X,Y,dZ,tend=100.)
    return X,Y,dZ
    
def make_dz_honshu_ucsb3_t100_truncate():
    # This is on truncated domain that was used originally but
    # too small for total source.  To check against past results.
    fname_subfaults = 'UCSB3.txt'
    fname_dtopo = 'UCSB3_t100_truncate.tt1'
    faultparams = {}
    faultparams['mx'] = 201
    faultparams['my'] = 301
    faultparams['xlower'] = 140
    faultparams['xupper'] = 144
    faultparams['ylower'] = 35
    faultparams['yupper'] = 41
    fm = read_subfault_model(fname_subfaults)
    X,Y,dZ = make_okada_dz(fm, faultparams)
    write_dz(fname_dtopo, X,Y,dZ,tend=100.)
    return X,Y,dZ
    
    
def make_dz_witht_honshu_ucsb3():
    fname_subfaults = 'UCSB3.txt'
    fname_dtopo = 'UCSB3_dynamic.tt1'
    faultparams = {}
    faultparams['mx'] = 301
    faultparams['my'] = 301
    faultparams['xlower'] = 140
    faultparams['xupper'] = 146
    faultparams['ylower'] = 35
    faultparams['yupper'] = 41
    fm = read_subfault_model(fname_subfaults)
    times = linspace(0,170,18)
    X,Y,dZ = make_okada_dz_witht(fm, faultparams, times, fname_dtopo)
    #write_dz(fname_dtopo, X,Y,dZ)
    return X,Y,dZ

    
if __name__=="__main__":
    print "No main program"
