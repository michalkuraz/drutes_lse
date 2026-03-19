# Makefile


FC     = gfortran
FFLAGS = -fimplicit-none  -fcoarray=single -fbounds-check -fbacktrace -g -g3 -fdefault-real-8 -O0 -finit-real=nan -Wsurprising
FFLAGS += -Werror=line-truncation

OBJS = typy.o globals.o core_tools.o debug_tools.o printtools.o geom_tools.o tools.o \
       initvals.o hydrofnc.o readtools.o main.o

TARGET = nour_model

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJS)

typy.o: typy.f90
	$(FC) $(FFLAGS) -c typy.f90

globals.o: globals.f90 typy.o
	$(FC) $(FFLAGS) -c globals.f90

core_tools.o: core_tools.f90 typy.o globals.o
	$(FC) $(FFLAGS) -c core_tools.f90

debug_tools.o: debug_tools.f90 typy.o globals.o core_tools.o
	$(FC) $(FFLAGS) -c debug_tools.f90

readtools.o: readtools.f90 typy.o globals.o core_tools.o debug_tools.o
	$(FC) $(FFLAGS) -c readtools.f90

printtools.o: printtools.f90 typy.o globals.o core_tools.o debug_tools.o
	$(FC) $(FFLAGS) -c printtools.f90

geom_tools.o: geom_tools.f90 typy.o globals.o printtools.o core_tools.o debug_tools.o
	$(FC) $(FFLAGS) -c geom_tools.f90

tools.o: tools.f90 typy.o globals.o geom_tools.o
	$(FC) $(FFLAGS) -c tools.f90

initvals.o: initvals.f90 typy.o globals.o tools.o
	$(FC) $(FFLAGS) -c initvals.f90

hydrofnc.o: hydrofnc.f90 typy.o globals.o tools.o
	$(FC) $(FFLAGS) -c hydrofnc.f90

main.o: main.f90 typy.o globals.o tools.o initvals.o hydrofnc.o geom_tools.o
	$(FC) $(FFLAGS) -c main.f90

clean:
	rm -f *.o *.mod $(TARGET)
	# on pure Windows/PowerShell, you can use instead:
	# del /Q *.o *.mod $(TARGET).exe 2>nul || true

.PHONY: all clean run

run: $(TARGET)
	@echo "Running $(TARGET) with ARGS='$(ARGS)'..."
	@./$(TARGET) $(ARGS)
