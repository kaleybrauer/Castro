PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE

DIM        = 2

COMP	   = gnu

USE_MPI    = TRUE
USE_OMP    = FALSE

USE_MHD    = FALSE

USE_FORT_MICROPHYSICS := FALSE
BL_NO_FORT := TRUE

# define the location of the CASTRO top directory
CASTRO_HOME  := ../../..

# This sets the EOS directory in $(MICROPHYSICS_HOME)/EOS
EOS_DIR     := gamma_law

# This sets the network directory in $(MICROPHYSICS_HOME)/Networks
NETWORK_DIR := general_null
NETWORK_INPUTS = simple.net

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_HOME)/Exec/Make.Castro
