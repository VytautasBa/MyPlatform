# Compiles and links project MyPlatform
PROJ = myplatform
OBJS = std_types.o std_lapack.o numer_matrix.o numer_interp.o numer_fft_row.o numer_fft.o resources.o spectroscopy.o prep_rules.o HEOM.o MAIN.o	# list of object files
FLIBS = -llapack -lblas # list of libraries
#LIBDIR = /usr/local/lib

CC = gfortran		# name of the compiler

# execution of the makefile

$(PROJ).exe : $(OBJS)
	$(CC) -o $(PROJ) $(OBJS) $(FLIBS)

std_types.o : std_types.F90
	$(CC) -c $(.TARGET) std_types.F90

std_lapack.o : std_lapack.F90
	$(CC) -c $(.TARGET) std_lapack.F90

numer_matrix.o : numer_matrix.F90
	$(CC) -c $(.TARGET) numer_matrix.F90
	
numer_interp.o : numer_interp.F90
	$(CC) -c $(.TARGET) numer_interp.F90

numer_fft_row.o : numer_fft_row.F90
	$(CC) -c $(.TARGET) numer_fft_row.F90

numer_fft.o : numer_fft.F90
	$(CC) -c $(.TARGET) numer_fft.F90

resources.o : resources.F90
	$(CC) -c $(.TARGET) resources.F90

spectroscopy.o : spectroscopy.F90
	$(CC) -c $(.TARGET) spectroscopy.F90
	
prep_rules.o : prep_rules.F90
	$(CC) -c $(.TARGET) prep_rules.F90

HEOM.o : HEOM.F90
	$(CC) -c $(.TARGET) HEOM.F90

MAIN.o : MAIN.F90
	$(CC) -c $(.TARGET) MAIN.F90

# cleaning up with 'make clean'
clean :
	rm -f *.mod *.o

# runing myplatform
run :
	./$(PROJ).exe

# cleans up stuff from previous runs	
fresh :
	rm -f ./out/*.dat