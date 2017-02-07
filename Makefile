### Compilers & flags
F90=mpifort

FFTWLIBS=~/bin/FFTW/lib/libfftw3.a
VECLIBSMACOSX=
LAPACKLIB=-L/opt/local/lib/lapack-3.5.0 -llapack -lblas
BLASLIB=/opt/local/lib/lapack-3.5.0/librefblas.a
PNETCDFLIBS=

LIBS    = 


EXE = exec
F90SRC = main.f90 constants.f90 modMPI.f90 MatrixVector.f90 modSpecies.f90 modMesh.f90 modAssign.f90 modRecord.f90 modPM1D.f90 random.f90 ArMCC.f90 init.f90 modBC.f90 modTarget.f90 modSource.f90 timeStep.f90 modAdj.f90 modQoI.f90  testmodule.f90 modFSens.f90
F90OBJ = main.o constants.o modMPI.f90 MatrixVector.o modSpecies.o modMesh.o modAssign.o modRecord.o modPM1D.o random.o ArMCC.o init.o modBC.o modTarget.o modSource.o timeStep.o modAdj.o modQoI.f90 testmodule.o modFSens.o

### Targets
all: $(EXE)
run: $(EXE) 
	./$(EXE)

# Link object files to executables
$(EXE): $(F90OBJ)
	$(F90) -o $(EXE) $(F90OBJ) $(LIBS)

# All .o files depend on the corresponding .f90 file
%.o: %.f90
	$(F90) -c $<

# Dependencies
random.o : constants.o
MatrixVector.o : constants.o
modMPI.o : constants.o
modSpecies.o : constants.o
modMesh.o : MatrixVector.o
modAssign.o : modSpecies.o modMesh.o
modPM1D.o : modSpecies.o modMesh.o modAssign.o
modRecord.o : modPM1D.o
modAdj.o : modPM1D.o
modQoI.o : modAdj.o
ArMCC.o : modPM1D.o modRecord.o random.o
modBC.o : modPM1D.o random.o
modFSens.o : modBC.o
init.o : modFSens.o
modTarget.o : modAdj.o random.o
modSource.o : modPM1D.o random.o
timeStep.o : modTarget.o modSource.o modFSens.o modRecord.o ArMCC.o modAdj.o modQoI.o
testmodule.o : init.o timeStep.o modMPI.o
main.o : testmodule.o

clean:
	rm *.o *.mod $(EXE)

.PHONY: all run clean


