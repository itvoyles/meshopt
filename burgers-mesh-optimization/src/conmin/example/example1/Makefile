DEFSDIR=../../../

BINDIR = $(DEFSDIR)/../bin/
SRCDIR = ./src
OBJDIR = ./obj
LIB = libconmin.so
VPATH = $(SRCDIR):$(OBJDIR)\
        $(DEFSDIR)/select_precision/src \
        $(DEFSDIR)/../lib/

#define source files and output files ---------------
EXE =  conmin_example1.exe
OBJ =  select_precision.o \
       conmin_ex1.o
       
#include standard definitions for compilation -------

include $(DEFSDIR)/std_defs.mk

#define dependencies ---------------------------------
conmin_ex1.f95 : 

$(EXE) : $(OBJ)
