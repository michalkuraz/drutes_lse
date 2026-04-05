# Makefile


FC     = gfortran
FFLAGS = -fimplicit-none  -fcoarray=single -fbounds-check -fbacktrace -g -g3 -fdefault-real-8 -O0 -finit-real=nan -Wsurprising
FFLAGS += -Werror=line-truncation

OBJS = typy.o  globals.o core_tools.o debug_tools.o printtools.o  tools.o \
       initvals.o hydrofnc.o readtools.o main.o 

TARGET = nour_model

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJS)

typy.o: src/typy.f90
	$(FC) $(FFLAGS) -c src/typy.f90
	
globals.o: src/globals.f90 typy.o
	$(FC) $(FFLAGS) -c src/globals.f90

core_tools.o: src/core_tools.f90 typy.o globals.o
	$(FC) $(FFLAGS) -c src/core_tools.f90

debug_tools.o: src/debug_tools.f90 typy.o globals.o core_tools.o
	$(FC) $(FFLAGS) -c src/debug_tools.f90

readtools.o: src/readtools.f90 typy.o globals.o core_tools.o debug_tools.o
	$(FC) $(FFLAGS) -c src/readtools.f90

printtools.o: src/printtools.f90 typy.o globals.o core_tools.o debug_tools.o
	$(FC) $(FFLAGS) -c src/printtools.f90


tools.o: src/tools.f90 typy.o globals.o geom_tools.o
	$(FC) $(FFLAGS) -c src/tools.f90

initvals.o: src/initvals.f90 typy.o globals.o tools.o
	$(FC) $(FFLAGS) -c src/initvals.f90

hydrofnc.o: src/hydrofnc.f90 typy.o globals.o tools.o
	$(FC) $(FFLAGS) -c src/hydrofnc.f90

main.o: src/main.f90 typy.o globals.o tools.o initvals.o hydrofnc.o geom_tools.o
	$(FC) $(FFLAGS) -c src/main.f90

clean:
	rm -f *.o *.mod $(TARGET)
	# on pure Windows/PowerShell, you can use instead:
	# del /Q *.o *.mod $(TARGET).exe 2>nul || true

.PHONY: all clean run

run: $(TARGET)
	@echo "Running $(TARGET) with ARGS='$(ARGS)'..."
	@./$(TARGET) $(ARGS)
