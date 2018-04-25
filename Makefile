### Compilers & flags
F90=mpifort

TOPDIR = $(dir $(firstword $(MAKEFILE_LIST)))
SRCDIR = src
OBJDIR = obj
MODDIR = mod

FFTWLIBS=/g/g92/chung28/Programs/FFTW/lib/libfftw3.a
VECLIBSMACOSX=
LAPACKLIB=-L/opt/local/lib/lapack-3.5.0 -llapack -lblas
BLASLIB=/opt/local/lib/lapack-3.5.0/librefblas.a
PNETCDFLIBS=

LIBS    = $(FFTWLIBS)

EXE = landau_Jtheta3
F90SRC = main.f90 \
		modInputHelper.f90 \
		constants.f90 \
		modMPI.f90 \
		MatrixVector.f90 \
		modSpecies.f90 \
		modMesh.f90 \
		modAssign.f90 \
		modRecord.f90 \
		modPM1D.f90 \
		modAdj.f90 \
		random.f90 \
		ArMCC.f90 \
		init.f90 \
		modBC.f90 \
		modQoI.f90 \
		modTarget.f90 \
		modSource.f90 \
		timeStep.f90 \
		timeStepAdj.f90 \
		timeStepFSens.f90 \
		testmodule.f90 \
		modFSens.f90 \
		PlasmaProblems.f90 \
		AdjointProblems.f90 \
		MCCProblems.f90 \
		FSensProblems.f90
F90OBJ = $(F90SRC:%.f90=$(OBJDIR)/%.o)

### Targets
all: dir $(EXE)
run: $(EXE) 
	./$(EXE)
dir: $(OBJDIR) $(MODDIR)

# Sub-directories
$(OBJDIR):
	mkdir $(TOPDIR)$(OBJDIR)
$(MODDIR):
	mkdir $(TOPDIR)$(MODDIR)

# Link object files to executables
$(EXE): $(F90OBJ)
	@echo $(dir $(firstword $(MAKEFILE_LIST)))
	$(F90) -o $(EXE) $(F90OBJ) $(LIBS)

# All .o files depend on the corresponding .f90 file
$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	$(F90) -J$(MODDIR) -c -o $@ $<

# Dependencies
$(OBJDIR)/random.o : $(OBJDIR)/constants.o
$(OBJDIR)/MatrixVector.o : $(OBJDIR)/constants.o
$(OBJDIR)/modMPI.o : $(OBJDIR)/constants.o
$(OBJDIR)/modInputHelper.o : $(OBJDIR)/modMPI.o
$(OBJDIR)/modSpecies.o : $(OBJDIR)/constants.o
$(OBJDIR)/modMesh.o : $(OBJDIR)/MatrixVector.o
$(OBJDIR)/modAssign.o : $(OBJDIR)/modSpecies.o \
						$(OBJDIR)/modMesh.o
$(OBJDIR)/modBC.o : $(OBJDIR)/modSpecies.o \
					$(OBJDIR)/modMesh.o \
					$(OBJDIR)/random.o
$(OBJDIR)/ArMCC.o : $(OBJDIR)/modSpecies.o \
					$(OBJDIR)/random.o
$(OBJDIR)/modPM1D.o : $(OBJDIR)/modBC.o \
						$(OBJDIR)/modAssign.o \
						$(OBJDIR)/ArMCC.o
$(OBJDIR)/modRecord.o : $(OBJDIR)/modPM1D.o
$(OBJDIR)/modAdj.o : $(OBJDIR)/modPM1D.o
$(OBJDIR)/modQoI.o : $(OBJDIR)/modAdj.o
$(OBJDIR)/modFSens.o : $(OBJDIR)/modBC.o \
						$(OBJDIR)/modRecord.o
$(OBJDIR)/init.o : $(OBJDIR)/modFSens.o
$(OBJDIR)/modTarget.o : $(OBJDIR)/modAdj.o \
						$(OBJDIR)/random.o
$(OBJDIR)/modSource.o : $(OBJDIR)/modPM1D.o \
						$(OBJDIR)/random.o
$(OBJDIR)/timeStep.o : $(OBJDIR)/modTarget.o \
						$(OBJDIR)/modSource.o \
						$(OBJDIR)/modFSens.o \
						$(OBJDIR)/modRecord.o \
						$(OBJDIR)/ArMCC.o \
						$(OBJDIR)/modAdj.o \
						$(OBJDIR)/modQoI.o
$(OBJDIR)/timeStepAdj.o : $(OBJDIR)/timeStep.o \
							$(OBJDIR)/modTarget.o
$(OBJDIR)/timeStepFSens.o : $(OBJDIR)/timeStep.o \
							$(OBJDIR)/modFSens.o
$(OBJDIR)/testmodule.o : $(OBJDIR)/init.o \
							$(OBJDIR)/timeStep.o \
							$(OBJDIR)/timeStepAdj.o \
							$(OBJDIR)/timeStepFSens.o \
							$(OBJDIR)/modMPI.o
$(OBJDIR)/PlasmaProblems.o : $(OBJDIR)/init.o \
								$(OBJDIR)/modSource.o \
								$(OBJDIR)/modTarget.o \
								$(OBJDIR)/timeStep.o \
								$(OBJDIR)/modMPI.o \
								$(OBJDIR)/modInputHelper.o
$(OBJDIR)/MCCProblems.o : $(OBJDIR)/init.o \
								$(OBJDIR)/modTarget.o \
								$(OBJDIR)/modMPI.o \
								$(OBJDIR)/timeStep.o
$(OBJDIR)/AdjointProblems.o : $(OBJDIR)/init.o \
								$(OBJDIR)/timeStepAdj.o \
								$(OBJDIR)/modMPI.o \
								$(OBJDIR)/modInputHelper.o
$(OBJDIR)/FSensProblems.o : $(OBJDIR)/init.o \
								$(OBJDIR)/timeStepAdj.o \
								$(OBJDIR)/timeStepFSens.o \
								$(OBJDIR)/modMPI.o \
								$(OBJDIR)/modInputHelper.o
$(OBJDIR)/main.o : $(OBJDIR)/testmodule.o \
					$(OBJDIR)/PlasmaProblems.o \
					$(OBJDIR)/AdjointProblems.o \
					$(OBJDIR)/MCCProblems.o \
					$(OBJDIR)/FSensProblems.o \
					$(OBJDIR)/modMPI.o \
					$(OBJDIR)/modInputHelper.o

clean:
	rm -rf $(OBJDIR) $(MODDIR) $(EXE)

.PHONY: dir all run clean


