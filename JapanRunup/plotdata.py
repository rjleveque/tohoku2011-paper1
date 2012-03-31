"""
Run this script to create plots of runup from:
  rawdata.txt   - observations
  NAME.txt - from source NAME
The latter can be created from fixed grid output using 
  runup.process_all()
"""

from pylab import *
import runup

xraw,yraw,etaraw = loadtxt('rawdata.txt',unpack=True)

figure(1,(9,8))
clf()
subplot(4,3,1)
plot(yraw,etaraw,'k.',markersize=2)
ylabel("Elevation (m)")
xlim(36,41)
ylim(0,40)
yticks([0,10,20,30,40])
text(36.5,32,"Data",fontsize=15)

models = """GCMT Hayes UCSB3 Ammon Caltech Fujii Saito Gusman GusmanAU PMELWei
          """.split()
modelnums = [str(i) for i in range(1,11)]
modelnums[-1] = '9'
modelnums[7] = '8a'
modelnums[8] = '8b'


for j,model in enumerate(models):
    subplot(4,3,j+2)
    plot(yraw,etaraw,'k.',markersize=2)

    xflat,yflat,etaflat = loadtxt('%s.txt' % model,unpack=True)
    plot(yflat,etaflat,color=[1.0,0.6,0.6],linewidth=2)
    xlim(36,41)
    ylim(0,40)
    yticks([0,10,20,30,40])

    text(36.5,32,modelnums[j],fontsize=15)
    if j in (2,5,8):
        ylabel("Elevation (m)")
    if j in (7,8,9):
        xlabel("Latitude (N)")


fname = "all_japan_runup.png"
savefig(fname)
print "Created ",fname


# Separate figures:

figure(2)
clf()
plot(yraw,etaraw,'k.',markersize=2)
ylabel("Elevation (m)")
xlabel("Latitude (N)")
xlim(36,41)
ylim(0,40)
yticks([0,10,20,30,40])
text(36.5,32,"Data",fontsize=15)
fname = "japan_runup_observations.png"
savefig(fname)
print "Created ",fname


for j,model in enumerate(models):
    clf()
    plot(yraw,etaraw,'k.',markersize=2)

    dirname = 'fgoutput/' + model
    xflat,yflat,etaflat = loadtxt('%s.txt' % model,unpack=True)
    plot(yflat,etaflat,color=[1.0,0.6,0.6],linewidth=2)

    xlim(36,41)
    ylim(0,40)
    yticks([0,10,20,30,40])

    text(36.5,32,modelnums[j],fontsize=15)
    ylabel("Elevation (m)")
    xlabel("Latitude (N)")

    fname = "japan_runup_source_%s.png" % model
    savefig(fname)
    print "Created ",fname

