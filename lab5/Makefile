# File: Makefile
# Author: T. Dufaud & G. Dos Santos
#
# Makefile for the code developed in TP5.

################################################################################
# Variables for this Makefile
################################################################################
# -- Compiler
CC=gcc
# -- Compiler Option
OPTC=-march=native -mtune=native -O3 -fomit-frame-pointer -fPIC -mavx -DAdd_ -DF77_INTEGER=int -DStringSunStyle
# -- Directories
TPDIR=.
TPDIRSRC=$(TPDIR)/src
# -- Librairies
LIBS=-llapacke -lblas -lm
# -- Include directories
INCLBLASLAPACK= -I /usr/include
INCL= -I $(TPDIR)/include $(INCLBLASLAPACK) 

################################################################################
# Makefile
################################################################################
OBJENV=deps/tp_env.o
OBJTP2ITER=deps/lib_poisson1D.o deps/tp2_poisson1D_iter.o
OBJTP2DIRECT=deps/lib_poisson1D.o deps/tp2_poisson1D_direct.o
OBJTP2LU=deps/lib_poisson1D.o deps/lib_lu_tridiag.o deps/tp2_lu_tridiag.o

all: bin/tp_testenv bin/tp2poisson1D_iter bin/tp2poisson1D_direct bin/tp2_lu_tridiag

testenv: bin/tp_testenv

tp2poisson1D_iter: bin/tp2poisson1D_iter

tp2poisson1D_direct: bin/tp2poisson1D_direct

tp2_lu_tridiag: bin/tp2_lu_tridiag

deps/tp_env.o: $(TPDIRSRC)/tp_env.c
	@mkdir -p deps/
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/tp_env.c -o $@ 

deps/lib_poisson1D.o: $(TPDIRSRC)/lib_poisson1D.c
	@mkdir -p deps/
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/lib_poisson1D.c -o $@
	
deps/lib_lu_tridiag.o: $(TPDIRSRC)/lib_lu_tridiag.c
	@mkdir -p deps/
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/lib_lu_tridiag.c -o $@

deps/tp2_poisson1D_iter.o: $(TPDIRSRC)/tp2_poisson1D_iter.c
	@mkdir -p deps/
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/tp2_poisson1D_iter.c -o $@

deps/tp2_poisson1D_direct.o: $(TPDIRSRC)/tp2_poisson1D_direct.c
	@mkdir -p deps/
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/tp2_poisson1D_direct.c -o $@

deps/tp2_lu_tridiag.o: $(TPDIRSRC)/tp2_lu_tridiag.c
	@mkdir -p deps/
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/tp2_lu_tridiag.c -o $@

bin/tp_testenv: $(OBJENV)
	$(CC) -o bin/tp_testenv $(OPTC) $(OBJENV) $(LIBS)

bin/tp2poisson1D_iter: $(OBJTP2ITER)
	$(CC) -o bin/tp2poisson1D_iter $(OPTC) $(OBJTP2ITER) $(LIBS)

bin/tp2poisson1D_direct: $(OBJTP2DIRECT)
	$(CC) -o bin/tp2poisson1D_direct $(OPTC) $(OBJTP2DIRECT) $(LIBS)

bin/tp2_lu_tridiag: $(OBJTP2LU)
	$(CC) -o bin/tp2_lu_tridiag $(OPTC) $(OBJTP2LU) $(LIBS)

run_testenv: bin/tp_testenv
	@$^

run_tp2poisson1D_iter: bin/tp2poisson1D_iter
	@$^

run_tp2poisson1D_direct: bin/tp2poisson1D_direct
	@$^

run_tp2_lu_tridiag: bin/tp2_lu_tridiag
	@$^

clean:
	rm -rf deps/ bin/* data/* *.dat
