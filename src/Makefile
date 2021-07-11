INCDIR   = inc
SRCDIR   = src
OBJDIR   = obj
LIBNAME  = libfalcON.a
CXXFLAGS = -I$(INCDIR) -I$(INCDIR)/utils -O2 -fPIC -mfpmath=sse -march=native -DfalcON_SINGLE -Wall

SOURCES  = \
    exception.cc \
    numerics.cc \
    basic.cc \
    body.cc \
    gravity.cc \
    kernel.cc \
    tree.cc

OBJECTS  = $(patsubst %.cc,$(OBJDIR)/%.o,$(SOURCES))

$(LIBNAME):  $(OBJECTS)
	$(AR) rc $(LIBNAME) $(OBJECTS)

$(OBJDIR)/%.o:  $(SRCDIR)/%.cc Makefile
	@mkdir -p $(OBJDIR)
	$(CXX) -c $(CXXFLAGS) -o "$@" "$<"

