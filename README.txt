
Simulation codes and data for the paper...

Comparison of earthquake source models for the 2011 Tohoku-oki event using
tsunami simulations and near field observations

by Breanyn T MacInnes, Aditya Riadi Gusman, Randall J LeVeque, Yuichiro
Tanioka 

Submitted to Bulletin of the Seismological Society of America, 31 March 2012
Published in 103(2013), pp. 1256-1274. See
* http://www.bssaonline.org/content/103/2B/1256.abstract?stoc

See also the electronic supplements at:
* http://www.bssaonline.org/content/103/2B/1256/suppl/DC1


topo directory contains topo files (bathymetry/topography)

sources directory contains dtopo files (change to bathy due to quake)
for the sources we have permission to redistribute:
  
    GDMT.txydz (source 1) from www.globalcmt.org/
    Hayes.txydz (source 2) (as in Hayes, 2011) is available at
      earthquake.usgs.gov/earthquakes/eqinthenews/2011/usc0001xgp/results/static_out
    UCSB3.txydz (source 3) (as in Shao et al., 2011) is available at
      www.geol.ucsb.edu/faculty/ji/big_earthquakes/2011/03/0311_v3/Honshu.html
    Ammon.txydz (source 4) (as in Ammon et al., 2011) is available at
      eqseis.geosc.psu.edu/~cammon/Japan2011EQ/
    Caltech.txydz (source 5) is available at
      http://www.tectonics.caltech.edu/slip_history/2011_taiheiyo-oki/
    Fujii.txydz (source 6) can be found in Fujii et al., 2011.
    Gusman.txydz (source 8a) can be found in (Gusman et al., in press) 
    GusmanAU.txydz (source 8b) can be found in (Gusman et al., in press) 
￼￼￼￼￼

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

   These runs were done with bathymetry that cannot be redistributed.
