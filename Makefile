# Last revision : 11-jul-2017 
#
# Modified for tgfdt version 2.1 
#(tgf pileup-dead-time simulation)
#
#  From a simple Makefile, for hits2hist 
#  by   M.Galli oct-2005
# 
# ------------------------------------------
#
#  defining options and general symbols :
#
#     directory of  the agile simulation system 
#     (in some versione from the environment) 
#
#MAINDIR = /home/galli/AGILE-SIMULATION
#MAINDIR = /home/mgalli/Agile_MCAL/localsoft  # for Nos 
MAINDIR = /home/marisaldi/localsoft

#   ----------- Dirs where include files and libs are placed 

LIBDIR = $(MAINDIR)/lib 
INCDIR = $(MAINDIR)/include
BINDIR = $(MAINDIR)/bin

#    ----------  libraries and includes from local soft 

#INCDIR =  
#LIBS = -L$(LIBDIR) -lhitlist  -lparametri 
#LIBS = -lz -L$(LIBDIR) -lhitlist -lgzstream 

# ----- cern root, includes and libs (root 5.34, in Deb 8, jessie) 

#LIBSROOT = -L/usr/lib/x86_64-linux-gnu -lRIO -lCore -lCint -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic
LIBSROOT = -L/prod_iasfbo/root/root_v5.34.24/lib -lRIO  -lCore -lCint -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic

INCSROOT = $(ROOTSYS)/include

# --- root 6.06, Deb8, jessie, Nos, compilato dai sources 

# da: root-config --libs
#-L/soft/Root/root60602/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic

#LIBSROOT = -L/soft/Root/root60602/lib  -lCore -lRIO -lCling -lTree -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic 
#INCSROOT = /soft/Root/root60602/include


#            --------  the fits libraries :
#    now from the environment 
#FITSDIR = /usr/local/fitsio

#MLIBFITSDIR += $(FITSDIR)/lib 
#MINCFITSDIR += $(FITSDIR)/include

#MFITSLIB = -L$(LIBFITSDIR) -lcfitsio 
#MFITSINC = -I$(INCFITSDIR) 

# FITS = -I$(INCFITSDIR)  $(FITSLIB)


#   --- general definitions 

GC = g++ 
# per root 6:  secondo: root-config --cflags
# FFLAGS = 
FFLAGS = -pthread -std=c++11 -Wno-deprecated-declarations -m64 

LDFLAGS = -m64
DEBUG = -g -O0 
NODEBUG = -O3 

#  ------------------------ project files: 

#  the header files 

HFILE = 

#  the cpp files  

#CPPFILE = mytgfdt_2_1B-comp-stats.cpp 
#CPPFILE = mytgfdt_2_1-comp.cpp 
#CPPFILE = mytgfdt_2_1-compth.cpp 

# the binary file 

#EXEFILE = tgfdt21Bstats
#EXEFILE = tgfdt21A
#EXEFILE = tgfdt21Ath

# caso con timers, t50 e fluence fissi, args to main 
#EXEFILE = tgfdt21Bmg1
EXEFILE = tgfdt22A
CPPFILE = tgf_dead_time_2_2-A.C



# --------------------------- the rules 

#                                assegnazione condizionale 
debug  :  FFLAGS += $(DEBUG)
tgfdt  :  FFLAGS += $(NODEBUG) 

# the  target is the default 

def : tgfdt ; @echo " =======>  All  done ...  "

#
#
#         ----------------------------------------------------
#                              dependencies:
#         ----------------------------------------------------
#

#.PHONY : debug
#debug: reset hits2root 

#$(EXEFILE) : $(CPPFILE)  

tgfdt  : $(EXEFILE)
#$(EXEFILE) : $(CPPFILE) MCVFileInfo.cpp MCVFileInfoDict.cpp
$(EXEFILE) : $(CPPFILE) 
	$(GC) $(DEBUG) $(FFLAGS) $(LDFLAGS) -I$(INCDIR) -I$(INCSROOT) -o  $@  $^  $(LIBS) $(LIBSROOT)

# generazione funzioni per I/O per classi custom, con rootcint 
# occorre il file header con le classi custom (qui MCVFileInfo.h )
# da cui rootcint crea i dictionary per l'I/O : *InfoDict.cpp  InfoDict.h 

#MVCFileInfoDict.cpp: MCVFileInfo.h ; rootcint -v MCVFileInfoDict.cpp -c MCVFileInfo.h

#  this target resets dependencies and force make to re-make all

.PHONY : reset                  # dico che reset non e' un file ... 
reset : ; touch   *.cpp *.h 

#  this target deletes object files 

.PHONY : clear
clear : ; rm *.o   
	
# install :

install:   ; cp $(EXEFILE) $(BINDIR) 

