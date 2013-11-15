DEFSDIR=../

BINDIR = ../../bin/
SRCDIR = ./src
OBJDIR = ./obj
VPATH = $(SRCDIR):$(OBJDIR)\
        ../select_precision/src \
        ../derivatives/src \
        ../misc/src \
        ../burgers-2d/src \
	../springsystem2d_force/fortran

LIB = conmin
LIBDIR= $(DEFSDIR)../lib

#define source files and output files ---------------
EXE =  burgers2d-te-conmin.exe
OBJ =  select_precision.o \
       misc_func.o \
       burgers2d_functions.o \
       derivative_transform.o \
       fd_derivatives.o \
       functional_integral_te_squared.o \
       derivatives.o \
       functional.o \
       burgers2d_conmin.o \
       springsystem2d_force.o
       
#include standard definitions for compilation -------

include $(DEFSDIR)std_defs.mk

#define dependencies ---------------------------------
misc_func.o : select_precision.o

functional_integral_te_squared.o : select_precision.o \
                                   misc_func.o

burgers2d_functions.o:select_precision.o \
                      fd_derivative_calc.o \
                      derivative_transform.o\
                      misc_func.o \
                      fd_derivatives.o \
                      derivatives.o

derivatives.o : select_precision.o \
                fd_derivatives.o

derivative_transform.o:select_precision.o

fd_derivative_calc.o:select_precision.o misc_func.o

fd_derivatives.o : select_precision.o


functional.o : select_precision.o \
               burgers2d_functions.o \
               functional_integral_te_squared.o

springsystem2d_force.o : select_precision.o \
                         misc_func.o
              
burgers2d_conmin.o : select_precision.o \
                     burgers2d_functions.o \
                     functional.o \
		     springsystem2d_force.o		     


$(EXE) : $(OBJ)
