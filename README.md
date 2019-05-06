MM 21/06/2017
MCAL front end simulator

-----------------------------

tgf_dead_time_2_1.C

SW ready for production.

int nrdm = 2,                                  number of simulations (tested up to 1000, 10000 in progress)
char *outfilename="./sim.root”,      name of the output root file
int idoffset=0                                 offset to add to the simulation number *
float minFluence=0.1,                    minimum fluence to simulate (cm^-1)
float maxFluence=2.,                     maximum fluence to simulate (cm^-1) 
float minT50=0.02,                        minimum duration T50 (ms) to simulate  <— NB: it’s milliseconds!
float maxT50=0.05,                       maximum duration T50 (ms) to simulate  <— NB: it’s milliseconds!
int verbose = 0,                             verbose level, if !0 print stuff for debug
int sourcetype=0,           type of source. 0: Gaussian shaped TGF (default); 
                                            1: unform distribution with constant rate, to compare with real bkg measurements; 
                                            2: fake patterns to compare with PB pspice simulations
float rate= 20000.          constant rate (Hz) for sourcetype==1
float aeff= 150.            effective area (cm2) for the chosen geometry



------------------------------

mytgfdt_2_2-A.C

Compiled version by MG, build using custom Makefile.txt


