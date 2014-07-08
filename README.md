upgradeAnalysis
===============
photon.h loads the ntuple and variables used in the analysis.
photon.cxx loops over all events and produces vectors/histograms of relevant quantities
params.h contains global variables (like cut values) also contains a few macros (dr calculation and effective sigma calculation)

sigmaEff.C isn't supposed to be there...
ntuple_PU50_age0.root is the hgg signal sample for testing our analysis code and plotting software
testout.root saves all histograms created by the analysis (not formatted)
