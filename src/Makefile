#################################################
#           set parameters and options          #
#################################################
MOD=../mod/
LOBJ=../obj/
BIN=../bin/
VPATH = $(BIN):$(LOBJ)

# set the compiler 
F95 = gfortran

# fftw3 locations
FFTW3_api = /usr/lib/x86_64-linux-gnu/
FFTW3 = /usr/include/

# netcdf.mod location
NETCDF = 

# set the flags
#FFLAGS =   -I../mod  -J../mod -O3  -I$(FFTW3) -I$(NETCDF)
FFLAGS =   -I../mod  -J../mod -O3  -I$(FFTW3) -I$(FFTW3_api) 

# libraries to link to
#LIBS =  -lfftw3 -lblas -llapack -lnetcdff
LIBS =  -lfftw3 -lblas -llapack 

# compile the testing programs? 
TEST=true

# set some parameters
MV=mv -f
.PHONY = clean

#####################################################
#             list of excutable programs            #
#####################################################

# main programs
SRCS = love_numbers.f90          \
       sea_level_fingerprint.f90 

# add in testing programs if desired
SRCS += test_special_functions.f90    \
        test_interp.f90               \
        test_quadrature.f90           \
        test_fftw3.f90                \
        test_spherical_harmonics.f90  \
        test_wigner.f90               \
        test_spherical_model.f90      \
        test.f90                      \
        test_mesh.f90                 \
        test_matrix.f90               


###########################################################
#         list of modules (in order of dependency)        #
###########################################################
MOD_SRCS = module_constants.f90            \
           module_physical_constants.f90   \
	   module_error.f90                \
           module_util.f90                 \
	   module_fftw3.f90                \
           module_interp.f90               \
	   module_special_functions.f90    \
	   module_quadrature.f90           \
           module_spherical_harmonics.f90  \
           module_spherical_model.f90      \
	   module_PREM.f90                 \
	   module_DECK.f90                 \
           module_mesh.f90                 \
	   module_matrix.f90               \
           module_force.f90                \
           module_linear_solver.f90        \
           module_sea_level.f90


###########################################################
#       These sections should not need to be edited       #
###########################################################


# set some variables automatically
EXCT     = $(SRCS:%.f90=%)
OBJS     = $(SRCS:%.f90=$(LOBJ)%.o)
MOD_OBJS = $(MOD_SRCS:%.f90=$(LOBJ)%.o)
BINEXCT  = $(EXCT:%=$(BIN)%)


all: $(EXCT) 

$(EXCT): %:  $(MOD_OBJS)  %.o
	$(F95) $(FFLAGS) $^ $(LIBS) -o $(BIN)$@

$(OBJS): $(LOBJ)%.o: %.f90 
	$(F95) $(FFLAGS) $^  -c
	@$(MV) *.o $(LOBJ)

$(MOD_OBJS): $(MOD_SRCS)	
	$(F95)   $(FFLAGS) -c $(MOD_SRCS) 
	@$(MV) *.o $(LOBJ)

clean:
	@rm $(LOBJ)*.o
	@rm $(MOD)*.mod
	@rm $(BINEXCT)
