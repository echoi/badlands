CFGFILE = $(TOP)/config/Makefile.cfg.docker_linux
include $(CFGFILE)

EXEC=bin/$(EXE)
SHAREDLIB=lib/libbadlands.so

BUILDDIR=$(TOP)/build
OBJDIR=$(BUILDDIR)/objects
MODDIR=$(BUILDDIR)/modules
LIBDIR=$(BUILDDIR)/lib

FFLAGS= $(FCFLAGS) -I$(MODDIR)
	
LDFLAGS= -L$(LIBDIR)/

