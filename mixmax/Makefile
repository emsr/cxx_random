ifndef CC
CC=gcc
endif

ifdef _N
NFLAG=-D_N=${_N}
endif

ifndef PREFIX
PREFIX=/usr/local
endif

ifeq ($(CC), icc)
CFLAGS=-std=c99 -O3 -funroll-loops  -Wall
else
CFLAGS=-std=c99 -O3 -funroll-loops -Wall -I ${PREFIX}/include/
endif

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
DYNLIBFLAGS=-dynamiclib -Wl,-dead_strip_dylibs -Wl,-headerpad_max_install_names,-undefined,dynamic_lookup,-compatibility_version,1.0,-current_version,1.0,-install_name,/usr/local/lib/libmixmax.dylib
dynlib=libmixmax.dylib
else
CFLAGS+= -fPIC -lc
DYNLIBFLAGS=-shared -fPIC -Wl,-dead_strip_dylibs
dynlib=libmixmax.so
endif


objects = mixmax.o
mainobj = driver_main.o
skipobj = driver_skip_charpoly.o
verobj = driver_verification.o

.c.o:
	${CC} ${CFLAGS}  ${NFLAG} -c $<

mixmax: $(objects) $(mainobj) mixmax.c driver_main.c mixmax.h
	${CC} ${CFLAGS} -o mixmax $(objects) $(mainobj)

# multithreaded demo
multi: mixmax.c mixmax.h driver_threads.c 
	${CC} ${CFLAGS}  ${NFLAG}  -o multi mixmax.c driver_threads.c -lpthread

all: mixmax verify testU01 lib check

verify: $(objects) $(verobj) mixmax.c mixmax.h driver_verification.c 
	${CC} ${CFLAGS}  ${NFLAG} -o verify $(objects) $(verobj)
	
check:  clean mixmax verify
	@ echo ""
	@ echo [Test 1] - Testing correct floating point output
	@ if [ -e N240_12000.diff ]; then rm N240_12000.diff; fi
	@ echo "12000" | ./mixmax 2> /dev/null > N240_12000_local.out
	@ echo "Outputting the first 12000 numbers into N240_12000_local.out"
	@ diff -w N240_12000.out N240_12000_local.out > N240_12000.diff
	@ echo "The diff file N240_12000.diff should be empty."
	@ /bin/ls -l N240_12000.diff
	@ if [ -s N240_12000.diff ]; then  echo "<<Files differ!>> Test failed."; else echo "Test 1 passed"; fi
	@ echo ""
	@ echo "[Test 2] - Writing and reading the state from file"
	${CC} ${CFLAGS}  ${NFLAG}  -o testIO driver_iotest.c mixmax.c
	@ ./testIO
	@ echo ""
	@ echo "[Test 3] - Checking the explicit MIXMAX matrix"
	@ if [ -e mathematica.txt.gz ]; then gunzip mathematica.txt.gz; fi
	@ echo 100 | ./verify 2> /dev/null > mathematica.local.txt
	@ diff -w mathematica.N240.txt mathematica.local.txt > mathematica.txt.diff
	@ echo "The diff file mathematica.txt.diff should be empty."
	@ /bin/ls -l mathematica.txt.diff
	@ if [ -s mathematica.txt.diff ]; then  echo "<<Files differ!>> Test failed."; else echo "Test 3 passed"; fi
	@ echo ""
	@ echo "[Test 4] - Optional"
	@ echo "You may now run the contents of mathematica.local.txt in Mathematica (yourself) to get an independent check of correctness."

gsl: mixmax.c mixmax.h driver_gsl.c
	${CC} ${CFLAGS} driver_gsl.c mixmax.c -lgsl -L${PREFIX}/lib -DHOOKUP_GSL=1  -o gsl
	@ echo "Run it with:\n GSL_RNG_SEED=123 GSL_RNG_TYPE=MIXMAX ./gsl"

testU01: mixmax.c mixmax.h driver_testU01.c 
	${CC} ${CFLAGS}  ${NFLAG} -o testU01 driver_testU01.c -ltestu01  mixmax.c

lib: mixmax.o mixmax.h
	@ echo "Dynamic library is ${dynlib}"
	ar rcs libmixmax.a mixmax.o
	gcc ${DYNLIBFLAGS} -o ${dynlib} mixmax.o

install: lib
	cp mixmax.h ${PREFIX}/include/
	cp libmixmax.a ${PREFIX}/lib/
	cp ${dynlib} ${PREFIX}/lib/

asm: mixmax.c mixmax.h
	${CC} ${CFLAGS}  ${NFLAG} -S mixmax.c

src: 
	zip -9 src.zip mixmax.c Makefile mixmax.c mixmax.h driver_testU01.c \
	driver_main.c mixmax_skip_N256.c mixmax_skip_N256.oldS.c mixmax_skip_N240.c mixmax_skip_N17.c mixmax_skip_N8.c\
	driver_verification.c driver_threads.c driver_gsl.c driver_iotest.c README.pdf N256_12000.out N240_12000.out mathematica.txt mathematica.N240.txt mixmax.cpp mixmax.hpp example.cpp
	
clean:
	- rm $(objects) $(mainobj) $(skipobj) $(verobj) verify N256_12000_local.out N256_12000.diff mixmax gsl \
	testU01 libmixmax.* a.out

plusplus:
	g++ -std=c++11 mixmax.cpp example.cpp
