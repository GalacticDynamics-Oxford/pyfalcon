# makefile for the basic version of Walter Dehnen's falcON gravity code
# and two interfaces for it: pure Python extension module and a community module for the AMUSE framework

CODELIB = src/libfalcON.a

all: $(CODELIB) amuse pyfalcon

$(CODELIB):
	cd src && make

# Python extension module
pyfalcon: $(CODELIB)
	python setup.py install --user


# if AMUSE is installed, compile the plugin and put it into the AMUSE community code folder
ifdef AMUSE_DIR

# standard amuse configuration include (do we need it?):
#-include $(AMUSE_DIR)/config.mk
AMUSE_WORKER_DIR = $(AMUSE_DIR)/src/amuse/community/falcon
AMUSE_WORKER     = $(AMUSE_WORKER_DIR)/falcon_worker
AMUSE_INTERFACE  = $(AMUSE_WORKER_DIR)/interface.py
AMUSE_WORKER_INIT= $(AMUSE_WORKER_DIR)/__init__.py
MPICXX   ?= mpicxx

amuse: $(AMUSE_WORKER)

$(AMUSE_WORKER_INIT):
	@mkdir -p $(AMUSE_WORKER_DIR)
	echo>>$(AMUSE_WORKER_DIR)/__init__.py

$(AMUSE_WORKER): amuse_interface.py amuse_interface.cpp $(CODELIB) $(AMUSE_WORKER_INIT)
	cp amuse_interface.py $(AMUSE_INTERFACE)
	cp amuse_example.py   $(AMUSE_WORKER_DIR)
	$(AMUSE_DIR)/build.py --type=H $(AMUSE_INTERFACE) FalconInterface -o "$(AMUSE_WORKER_DIR)/worker_code.h"
	$(AMUSE_DIR)/build.py --type=c $(AMUSE_INTERFACE) FalconInterface -o "$(AMUSE_WORKER_DIR)/worker_code.cpp"
	$(MPICXX) -o "$@" "$(AMUSE_WORKER_DIR)/worker_code.cpp" amuse_interface.cpp -Isrc/inc -Isrc/inc/utils -DfalcON_SINGLE -O2 $(CODELIB) $(MUSE_LD_FLAGS)

else
amuse:
	echo "AMUSE is not installed"
endif
