#   Makefile for CLASTOOL code.
#  
# Notice:
#
# At CEBAF you must execute the following before this will compile and link.
#
# SETUP YOUR CLAS ENVIRONMENT (source /group/clas/builds/../packages/cms/jlab.csh )
# SETUP YOUR ROOT ENVIRONMENT (use root       )
#
#
# REWRITE:
#
#	Maurik.Holtrop@unh.edu	May-12-2004
# 
#
# NOTES:
#
# SRC_CLASS: contains the source files (*.cc) of the CLASSES. A dictionary is made.
# SRC_FILES: contains the source files (*.cc) for the FUNCTIONS and PROGRAMS. No dictionary.
# INC_FILES: contains include files. This is only used to detect if we need to remake the Makefile_depends.
#
# The depends file will only be used if you keep the last include: 'include Makefile_depends'
# This scheme will still try to generate one (in Makefile_bottom).
#
#
#
SRC_CLASS = $(wildcard T*.cc)

include ../Makefile_top

#
CREATED_SHLIB=libClasBanks.$(DllSuf)
CREATED_LIB=libClasBanks.a

OTHERLOADLIBS = -L$(ROOTLIBDIR)  -lMatrix
#
# Specify the default
#
all: shlib
#
include ../Makefile_bottom
include Makefile_depends
