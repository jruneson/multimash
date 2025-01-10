MODULE=mashf90
# If you add any further potential module files, add them to the next line (must appear later than pes.o)
OBJ=types.o maths.o pes.o linvib.o frexc.o tully.o mash.o

COMPILER=gfortran
FCOMPILER=gnu95
LIBS=-llapack -lblas

# Choose between debug (tells which line the code breaks on) and fast (parallel and optimized)
debug: FC=$(COMPILER) -g -fPIC -Wall -Wextra -Wconversion -pedantic --check=all
debug: F2PY=f2py --fcompiler=$(FCOMPILER)
debug: main

fast: FLAGS=-O2 -fopenmp -fPIC
fast: FC=$(COMPILER) $(FLAGS)
fast: F2PY=f2py --fcompiler=$(FCOMPILER) --f90flags='$(FLAGS)' -lgomp --quiet
fast: main

main: $(MODULE).so

$(MODULE).so: $(OBJ) f2py.f90
	$(F2PY) -c -m $(MODULE) $^ $(LIBS)

%.o: %.f90
	$(FC) -c $^

%.o: %.f
	$(FC) -c $^

clean:
	rm -f -r *.o *.so *.mod
