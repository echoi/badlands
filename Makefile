# ============================================================================
# Name        : Makefile
# Author      : tristan salles
# Copyright (C) 2015 
# Description : Makefile for BADLANDS
# ============================================================================
TOP=$(shell pwd)
CONFFILE= $(TOP)/config/Makefile.inc

include $(CONFFILE)

DIRMODS= GlobalUtils GeoMesh OceanModel EarthModel SpmModel CouplerModel

SOURCES = badlands_App.f90
OBJS=$(SOURCES:.f90=.o)

.PHONY : all dist plugin dust clobber

all: dist

dist: 
	@echo
	@echo "*************************************************"
	@echo "BADLANDS Author Tristan Salles "
	@echo "*************************************************"
	@echo
	@mkdir -p $(BUILDDIR)	
	@mkdir -p $(OBJDIR)
	@mkdir -p $(MODDIR)
	@mkdir -p $(LIBDIR)
	@mkdir -p bin
	for i in $(DIRMODS) ; do   \
    	  ( cd $$i ; make dist) ;       \
	done
	@echo "*************************************************"
	@echo	
	@echo "Build BADLANDS binary."
	@echo	
	@echo "*************************************************"
	@$(if $(wildcard badlands_App.o),rm -f badlands_App.o,)	
	make $(EXEC)

$(EXEC) :	$(OBJS)
	$(ESMF_F90LINKER) $(ESMF_F90LINKOPTS) $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS) $(FFLAGS) $(FOXFLAGS) $(H5FLAGS) $(ZOLTANFLAGS) -o $@ $^ \
	$(LDFLAGS) -lBADLANDS $(H5LDFLAGS) $(H5LIBS) $(ZOLTANLDFLAGS) $(ZOLTANLIBS) $(METISLDFLAGS) $(METISLIBS) $(LDFOXFLAGS) $(ESMF_F90ESMFLINKLIBS)
	@echo "*************************************************"
	@echo	
	@echo "BADLANDS updated in ./bin/."
	@echo	
	@echo "*************************************************"

%.o : %.f90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREENOCPP) $(FFLAGS) $(FOXFLAGS) $(H5FLAGS) $(ZOLTANFLAGS) $< -o $@ 
	$(AR) $(LIBDIR)/libBADLANDS.a $(OBJDIR)/*.o
	
dust :
	for i in $(DIRMODS) ; do   \
    	( cd $$i ; make dust) ;       \
	done
	$(foreach module,$(MODULES),cd $(module); cd - ; )
	rm -fv *~ *.bak *.o *.mod *.original

clobber : dust
	for i in $(DIRMODS) ; do   \
    	( cd $$i ; make clobber) ;   \
	done	
	rm -rfv $(BUILDDIR)
