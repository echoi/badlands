# ============================================================================
# Name        : Makefile
# Author      : tristan salles
# Copyright (C) 2015 
#
# Description : Makefile for GeoMesh
# ============================================================================
UP=..
TOP=$(shell pwd)/$(UP)
CONFFILE= $(TOP)/config/Makefile.inc

include $(CONFFILE)

FTREE = ./TreeSearch
FTRIA = ./Triangulator

TARGET = GeoMesh
TARGETLIB = $(LIBDIR)/lib$(TARGET).a

FSOURCES = $(FTREE)/kdtree.f90 $(FTREE)/kdtree2.f90 \
	ReadTopology.f90 ReadGeomorphology.f90 \
	ReadForces.f90 VisualiseGeo.f90 Geometry.f90

CSOURCES = $(FTRIA)/triangle.c $(FTRIA)/TriangleCall.c

FOBJS=$(FSOURCES:.f90=.o)
COBJS=$(CSOURCES:.c=.o)

.PHONY: dist module dust clobber 

dist : module

module :	
	@echo "*************************************************"
	@echo	
	@echo "Build GeoMesh modules and library."
	@echo	
	@echo "*************************************************"
	make $(TARGETLIB)
	
$(TARGETLIB): $(COBJS) $(FOBJS)
	cp *.o $(OBJDIR)
	cp $(FTREE)/*.o $(OBJDIR)
	cp $(FTRIA)/*.o $(OBJDIR)
	for f in *.mod ; do   \
    	if [ -f $$f ]; then \
    		mv $$f $(MODDIR);       \
		fi; 		\
	done
	for f in $(FTREE)/*.mod ; do   \
    	if [ -f $$f ]; then \
    		mv $$f $(MODDIR);       \
		fi; 		\
	done

%.o : %.f90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREENOCPP) $(FFLAGS) $(FOXFLAGS) $(H5FLAGS) $(ZOLTANFLAGS) $< -o $@ 
	
%.o : %.c
	$(BADLANDS_C) -c $(ESMF_COMPILEOPTS) $(ESMF_COMPILEPATHSLOCAL) $(ESMF_COMPILEPATHS) $(ESMF_COMPILECPPFLAGS) $(CFLAGS) $< -o $@ 

dust:
	rm -fv *.bak *~ *.o *.mod $(FTRIA)/*.o *.original $(FTRIA)/*.original $(FTREE)/*.o $(FTREE)/*.original
	
clobber: dust
