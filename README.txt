
Simulation codes and data for the paper...

Comparison of earthquake source models for the 2011 Tohoku-oki event using
tsunami simulations and near field observations

by Breanyn T MacInnes, Aditya Riadi Gusman, Randall J LeVeque, Yuichiro
Tanioka 

Submitted to Bulletin of the Seismological Society of America, 31 March 2012

See also the electronic supplements at:
  http://faculty.washington.edu/rjl/macinnes-esupp/


topo directory contains topo files (bathymetry/topography)

sources directory contains dtopo files (change to bathy due to quake)

dart directory contains scripts for running code and processing DART buoy data.

    GaugeErrors2.py is module containing all analysis codes.  
        >>> import GaugeErrors2 as GE
        >>> GE.make_all()

    dart/simulation directory contains code to run the dart simulations.
       run_all.py  to run all simulations (with each source)
      
       These runs were performed with GeoClaw from
        https://github.com/clawpack/clawpack-4.x
       Revision: 9317a769dfab7f5b78fa1a6b037750353f029976


inundation directory

       These runs were performed with GeoClaw from
        https://github.com/clawpack/clawpack-4.x
       Revision: 6a0021ceb7b1200220ddffb57c8d1a538281343f


JapanRunup directory contains runup results along the coast of Japan
   plotdata.py plots this data and creates datapoints.png

