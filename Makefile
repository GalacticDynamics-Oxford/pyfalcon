# makefile for the basic version of Walter Dehnen's falcON gravity code
# and two interfaces for it: pure Python extension module and a community module for the AMUSE framework

all: pyfalcon amuse

# Python extension module
pyfalcon:
	python setup.py install --user

# if AMUSE is installed, compile the plugin and put it into the AMUSE community code folder
ifeq ($(origin AMUSE_DIR), undefined)
  AMUSE_DIR := $(shell amusifier --get-amuse-dir)
# AMUSE_DIR is not necessarily the root folder of the python package itself, but a folder containing build.py
endif

ifdef AMUSE_DIR
-include $(AMUSE_DIR)/config.mk

# AMUSE_ROOT will contain the root folder of the python package (e.g. [...]/site-packages/amuse/ )
ifeq ($(AMUSE_ROOT),)
AMUSE_ROOT = $(shell python3 -c "import amuse; print(amuse.__path__[0])")
endif
ifeq ($(AMUSE_ROOT),)
# "python3" didn't work - try "python"
AMUSE_ROOT = $(shell python -c "import amuse; print(amuse.__path__[0])")
endif
ifeq ($(AMUSE_ROOT),)
# still didn't work...
$(warning Cannot determine AMUSE root folder: compiled module will remain in the current folder)
AMUSE_WORKER_DIR = $(shell pwd)
else
# compiled module will be placed in this folder
AMUSE_WORKER_DIR = $(AMUSE_ROOT)/community/falcon
endif
MPICXX ?= mpicxx
CODE_GENERATOR  ?= $(AMUSE_DIR)/build.py
AMUSE_WORKER     = $(AMUSE_WORKER_DIR)/falcon_worker
AMUSE_INTERFACE  = $(AMUSE_WORKER_DIR)/interface.py
AMUSE_WORKER_INIT= $(AMUSE_WORKER_DIR)/__init__.py

CODELIB = src/libfalcON.a

amuse: $(AMUSE_WORKER)

$(CODELIB):
	cd src && make

$(AMUSE_WORKER): interface.py amuse_interface.cpp $(CODELIB)
	-mkdir -p $(AMUSE_WORKER_DIR)
	touch $(AMUSE_WORKER_DIR)/__init__.py
	-cp interface.py $(AMUSE_INTERFACE)
	-cp amuse_example.py $(AMUSE_WORKER_DIR)
	$(CODE_GENERATOR) --type=H $(AMUSE_INTERFACE) FalconInterface -o "$(AMUSE_WORKER_DIR)/worker_code.h"
	$(CODE_GENERATOR) --type=c $(AMUSE_INTERFACE) FalconInterface -o "$(AMUSE_WORKER_DIR)/worker_code.cpp"
	$(MPICXX) -o "$@" "$(AMUSE_WORKER_DIR)/worker_code.cpp" amuse_interface.cpp -Isrc/inc -Isrc/inc/utils -DfalcON_SINGLE -O2 $(CODELIB) $(MUSE_LD_FLAGS)

else
amuse:
	echo "AMUSE is not installed"
endif
