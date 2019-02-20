# Set variables FC the Fortran compiler of your choice. Only tested in gfortran
FC = gfortran

# Set the compilation flags
FFLAGS  = -O2 -ffixed-line-length-0 -march=nocona


# Set variable POPDYN with the absolute path of your code and binaries.
# The value shown expands to the path of current working directory.
POPDYN = $(CURDIR)

#BIN      = $(POPDYN)/bin 
#SRC      = $(POPDYN)/src

# If there are several versions of the code, determine here which
# files should be compiled without extension. Also, provide a name to the output file
FILEOUT      = BeyondMeanField_gfortran
FILE1        = 1Beyond-MeanField_2-3
FILE2        = 2Parameters
FILE3        = 3Integration_Adaptive
FILESOURCE   = $(FILE1).f $(FILE2).f $(FILE3).f

# Stop your modifications here
# Build the code
$(FILEOUT): $(FILESOURCE)
	$(FC) $(FFLAGS) -o $(FILEOUT) $(FILESOURCE)
$(FILE1).o: $(FILE1).f
	$(FC) -c $(FFLAGS) $(FILE1).f
$(FILE2).o: $(FILE2).f
	$(FC) -c $(FFLAGS) $(FILE2).f
$(FILE3).o: $(FILE3).f
	$(FC) -c $(FFLAGS) $(FILE3).f
#Cleaning everything
#clean:
#	rm $(FILESOURCE)
# End of the makefile
