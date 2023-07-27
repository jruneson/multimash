MODULE=mashf90
# If you add any further potential module files, add them to the next line (must appear later than pes.o)
OBJ=types.o pes.o linvib.o bilinvib.o frexc.o tully.o mash.o
LIBS=-llapack

debug: FC=gfortran -g -fPIC -Wall -Wextra -Wconversion -pedantic --check=all
debug: F2PY=f2py --fcompiler=gnu95
debug: main

fast: FLAGS=-O2 -fopenmp -fPIC
fast: FC=gfortran $(FLAGS)
fast: F2PY=f2py --fcompiler=gnu95 --f90flags='$(FLAGS)' -lgomp --quiet
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
