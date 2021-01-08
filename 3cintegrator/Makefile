

ifndef MACH
	MACH := $(shell uname)
endif

ifeq ($(MACH), Linux)
#==== Linux system
# choose compiler
	include Linux_g++_make.inc
#	include Linux_icpc_make.inc
#if unusal loction specify below (here is for nersc)
#	FFTWLIB   = -L$(FFTW_DIR)
#	FFTWINC   = -I$(FFTW_INC)
else
#==== Darwin system
# for clang needs OpenMP support tesed with clang 3.9
#	include Darwin_clang_make.inc
	include Darwin_g++_omp_make.inc
	FFTWDIR  = /opt/local
	FFTWLIB   = -L$(FFTWDIR)/lib
	FFTWINC   = -I$(FFTWDIR)/include
endif

# ==== Profiling ====
PROFIL=
ifdef PROFILING
	PROFIL = -DPROFILING
endif


FFTWLIBN  = $(FFTWLIB) -lfftw3 -lfftw3_threads
#FFTWLIBN = -L$(FFTWLIB) -lfftw3 -lfftw3_omp



OBJ = ./Objs/
EXE = ./bin/
INC = ./inc/
INCEXT = $(INC)3CInt/
SRC = ./src/
LIB = ./lib/


BOOSTDIR = $(INC)
BOOSTINC = -I$(BOOSTDIR)


CPPFLAGS  += $(PROFIL) $(BOOSTINC) $(FFTWINC)
LDFLAGS   += $(FFTWLIBN) -lm


#  Define our target list
.PHONY: default
default: makedir lib 3cint

.PHONY: all 
all : makedir lib 3cint

.PHONY: tidy
tidy : 
	find . -name "*~" | xargs -I {} rm {}

.PHONY: clean
clean :
	rm -rf ./$(OBJ)/* ./$(EXE)/* ./$(LIB)/*

.PHONY: makedir
makedir :
	@mkdir -p $(OBJ)
	@mkdir -p $(EXE)
	@mkdir -p $(LIB)


#C++ common Objects
CXXOBJ = $(OBJ)walltimer.o \
	$(OBJ)walltime_c.o


CXXSHOBJ = walltimer.o \
	walltime_c.o


#C++ common Headers
CXXHDR =  $(INCEXT)walltimer.h \
	$(INCEXT)walltime_c.h \
	$(INCEXT)angpow_numbers.h \
	$(INCEXT)angpow_exceptions.h \
	$(INCEXT)angpow_func.h  \
	$(INCEXT)angpow_fft.h \
	$(INCEXT)3cint_chefunc.h \
	$(INCEXT)3cint_chealgo.h




#C++ rule for compiling
$(OBJ)%.o: $(SRC)%.cc $(CXXHDR)
	echo "compile... $@"
	$(CXXCOMPILE) -I$(INC) -c $< -o $@

######################
.PHONY: sharelib
sharelib : $(CXXOBJ)
	echo "make shared library..."
	cd $(OBJ); \
	$(CMDSHLCXX) -o ../$(LIB)lib3cint.$(SLEXT) $(CXXSHOBJ) $(LDFLAGS)

.PHONY: lib
lib : $(LIB)lib3cint.a

$(LIB)lib3cint.a : $(CXXOBJ)
	$(AR) $(ARFLAGS) $@ $^

###############
.PHONY:  3cint
3cint: $(EXE)3cint
	echo '---- 3cint made'

$(OBJ)3cint.o: 3cint.cc $(CXXHDR)
	echo "compile... $@"
	$(CXXCOMPILE)  -I$(INC) -c $< -o $@ 

$(EXE)3cint :  $(OBJ)3cint.o $(LIB)lib3cint.a
	echo "Link..."
	$(CXXLINK) -o $@ $(OBJ)3cint.o -L$(LIB) -l3cint  $(LDFLAGS)








