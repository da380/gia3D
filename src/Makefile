#################################################
#             define directory paths            #
#################################################
MOD=../mod/
LOBJ=../obj/
BIN=/home/david/codes/bin/
VPATH = $(BIN):$(LOBJ)

# set the compiler and flags
F95 = gfortran
FFLAGS =   -I../mod  -J../mod  -O5 -fbounds-check
LIBS =  -l lapack

# set some parameters
MV=mv -f
.PHONY = clean

#####################################################
#             list of excutable programs            #
#####################################################
SRCS = test.f90       


###########################################################
#         list of modules (in order of dependency)        #
###########################################################
MOD_SRCS = module_constants.f90            \
           module_physical_constants.f90   \
	   module_error.f90                \
           module_util.f90                 \
           module_interp.f90               \
	   module_quadrature.f90           \
           module_spherical_harmonics.f90  \
           module_spherical_model.f90


###########################################################
# These sections should not need to be edited unless a    #
#program needs non-standard files to be linked            #
###########################################################


# set some variables automatically
EXCT     = $(SRCS:%.f90=%)
OBJS     = $(SRCS:%.f90=$(LOBJ)%.o)
MOD_OBJS = $(MOD_SRCS:%.f90=$(LOBJ)%.o)


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

