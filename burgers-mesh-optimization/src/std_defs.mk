FC = gfortran
MV = mv
FFLAGS =


.PHONY : build
build:$(EXE)

$(EXE):$(OBJ)
	@echo
	@echo 'Building' $(EXE)
	$(FC) $(FFLAGS) -o $(EXE) -J$(OBJDIR) -I$(OBJDIR) $(addprefix $(OBJDIR)/,$(OBJ)) -L$(LIBDIR) -l$(LIB)

.PHONY:install

install:$(EXE)
	@echo
	@echo 'Installing' $(EXE) 'to' $(BINDIR)
	@$(MV) $(EXE) $(BINDIR)


.PHONY : compile
compile : $(OBJ)

$(OBJ) : %.o : %.f95
	@echo 'Compiling' $<
	@echo '      -->' $(OBJDIR)/$@
	@$(FC) $(FFLAGS) -c $< -o $(OBJDIR)/$@ -I$(OBJDIR) -J$(OBJDIR)



.PHONY : print
print :
	@echo
	@echo $VPATH
	@echo 'Objects:' $(OBJ)


.PHONY : clean
clean :
	@echo
	@echo 'Cleaning up...' 
	@rm -fv $(OBJDIR)/*.o
	@rm -fv $(OBJDIR)/*.mod
