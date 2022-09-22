#################################################
#             define directory paths            #
#################################################
MOD=../mod/
LOBJ=../obj/
BIN=../bin/
VPATH = $(BIN):$(LOBJ)

# set the compiler 
F95 = gfortran


# fftw3 f90 location
FFTW3 = /home/david/lib/fftw-3.3.3/api
NETCDF = /usr/include/

# set the flags
FFLAGS =   -I../mod  -J../mod -O3  -I$(FFTW3) -I$(NETCDF)


# libraries to link to
LIBS =  -lfftw3 -lblas -llapack -lnetcdff



# set some parameters
MV=mv -f
.PHONY = clean

#####################################################
#             list of excutable programs            #
#####################################################
SRCS = test_special_functions.f90    \
       test_quadrature.f90           \
       test_fftw3.f90                \
       test_spherical_harmonics.f90  \
       test_wigner.f90               \
       test_spherical_model.f90      \
       test.f90                      \
       test_mesh.f90                 \
       test_matrix.f90               \
       load_love_numbers.f90         \
       tide_love_numbers.f90         \
       netcdf_test.f90


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
           module_linear_solver.f90


###########################################################
# These sections should not need to be edited unless a    #
#program needs non-standard files to be linked            #
###########################################################


# set some variables automatically
EXCT     = $(SRCS:%.f90=%)
OBJS     = $(SRCS:%.f90=$(LOBJ)%.o)
MOD_OBJS = $(MOD_SRCS:%.f90=$(LOBJ)%.o)
BINEXCT     = $(EXCT:%=$(BIN)%)


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
