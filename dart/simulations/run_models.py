
import os
from pyclaw import runclaw
from pyclaw.plotters.data import ClawPlotData
import setrun
import setplot
import numpy as np

dtopodir = '/Users/rjl/git/tohoku2011a/sources/'

gaugenos = [21401, 21413, 21418, 21419]
models = """GCMT Hayes UCSB3 Ammon Caltech Fujii Saito Gusman GusmanAU
            pmelWei""".split()

models = ["GCMT"]

for model in models:
    print "Model: ",model

    rundata = setrun.setrun()
    rundata.geodata.dtopofiles = [[1,4,4,dtopodir+'%s.txydz' % model]]
    rundata.write()  # create *.data files

    outdir = '_output_%s' % model
    runclaw.runclaw('xgeoclaw',outdir)

    pd = ClawPlotData()
    pd = setplot.setplot(pd)
    pd.outdir = outdir
    for gaugeno in gaugenos:
        g = pd.getgauge(gaugeno)
        d = np.vstack([g.t,g.q[:,3]]).T
        fname = '../simulation_results/%s_%s.txt' % (model,gaugeno)
        np.savetxt(fname, d)
        print "Created ",fname


