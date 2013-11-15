DEFSDIR=../

BINDIR = ../../bin/
SRCDIR = ./src
OBJDIR = ./obj
VPATH = $(SRCDIR):$(OBJDIR)\
        ../select_precision/src \
        ../derivatives/src \
        ../misc/src

#define source files and output files ---------------
EXE =  burgers1d-te.exe
OBJ =  select_precision.o \
       misc_func.o \
       burgers1d_functions.o \
       derivative_transform.o \
       fd_derivatives.o \
       fd_derivative_calc.o \
       burgers1d-te.o \
       burgers1d_testroutines.o \
       functional_integral_te_squared.o\
       functional_integral_te_4th.o \
       functional_integral_te_1st.o \
       derivatives.o\
       functional.o \
       calc_grid.o
       
#include standard definitions for compilation -------

include $(DEFSDIR)/std_defs.mk

#define dependencies ---------------------------------
misc_func.o:select_precision.o

functional_integral_te_squared.o : select_precision.o \
                                   misc_func.o

functional_integral_te_4th.o : select_precision.o \
                               misc_func.o

functional_integral_te_1st.o : select_precision.o \
                                  misc_func.o

burgers1d_functions.o:select_precision.o \
                      fd_derivative_calc.o \
                      derivative_transform.o\
                      derivatives.o

derivatives.o : select_precision.o \
                fd_derivatives.o

burgers1d_testroutines.o : select_precision.o

derivative_transform.o:select_precision.o

fd_derivative_calc.o:select_precision.o misc_func.o

burgers1d-te.o : select_precision.o \
                 burgers1d_functions.o \
                 derivative_transform.o \
                 fd_derivatives.o \
                 fd_derivative_calc.o \
                 burgers1d_testroutines.o \
                 misc_func.o \
                 functional_integral_te_squared.o \
                 functional_integral_te_4th.o \
                 functional_integral_te_1st.o \
                 calc_grid.o \
                 functional.o

functional.o : select_precision.o \
               burgers1d_functions.o \
               derivatives.o \
               functional_integral_te_squared.o \
               fd_derivative_calc.o \
               fd_derivatives.o


calc_grid.o : select_precision.o \
              fd_derivative_calc.o \
              burgers1d_functions.o

$(EXE) : $(OBJ)
