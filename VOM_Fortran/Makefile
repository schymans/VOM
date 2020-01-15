#Makefile to compile the VOM model
#Usage (linux systems):
#Adjust the compiler to your own compiler and set the paths
#Open terminal; go to the directory and type "make"
#To remove old compiled code, type "make clean" 
#Written by: R.C. Nijzink june 2014
#References: 
#Adjustments:  

#The name of the exetuble (the final program)
PROGRAM  = model.x


FC      := /usr/bin/gfortran #compiler name and directory
#FC      := /usr/bin/x86_64-w64-mingw32-gfortran #64bit


SRC_PATH = ./src		#path of source files
OBJ_PATH = ./tmp		#path .o and .mod files



SRC_PATH    := $(abspath $(SRC_PATH:~%=${HOME}%))
OBJ_PATH   := $(abspath $(OBJ_PATH:~%=${HOME}%))

export OMP_NUM_THREADS=2

# compiler flags and options
FCFLAGS = 
FLFLAGS =

NC_INCLUDE=/usr/include
NC_LIB := /usr/lib/

#source objects
SRC = VOM-code/modules.f90 VOM-code/readdata.f90 VOM-code/writedata.f90 VOM-code/coreprog.f90 VOM-code/sample.f90 VOM-code/sce.f90 VOM-code/transpmodel.f90 VOM-code/watbal.f90 

TESTS = VOM_test/test_tmp/model.x

#mod files
MSRC = *.mod 
    

all: $(PROGRAM) 
		
$(PROGRAM): $(SRC) 
	$(FC) -fopenmp -o $@ $^ -static-libgfortran -I${NC_INCLUDE} -L${NC_LIB} -lnetcdff; rm $(MSRC)    #;mv $(MSRC) $(OBJ_PATH)
#	$(FC) -o $@ $^ -I${NFDIR}/include -L${NFDIR}/lib  -static ;mv $(MSRC) $(OBJ_PATH)
%.o: %.f
	$(FC) -o $@ $<




.PHONY: clean

clean:
	rm -f *.o *.mod *.MOD *.exe $(OBJ_PATH)/*.o $(OBJ_PATH)/*.mod $(OBJ_PATH)*.MOD $(OBJ_PATH)*.exe $(PROGRAM)   

check:
	@echo "Running checks..."; \
    mkdir VOM_test/test_tmp ; \
    mkdir VOM_test/test_tmp/input ; \
    mkdir VOM_test/test_tmp/output ; \
    cd VOM_test/test_tmp/ ; \
    cp ../ref_test/dailyweather.prn input/dailyweather.prn ; \
	cp ../ref_test/vom_namelist vom_namelist;\
	cp ../ref_test/check_outcomes.py check_outcomes.py;\
	cp ../ref_test/pars.txt input/pars.txt; \
	cp ../../model.x model.x; \
    ./model.x > /dev/null; \
    echo "Testing if results are different ..."; \
	python check_outcomes.py || exit 1; cd ../../ ; \
	rm -r VOM_test/test_tmp/; \
 










