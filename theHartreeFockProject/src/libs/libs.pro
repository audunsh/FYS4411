include(../../defaults.pri)
TEMPLATE = lib
TARGET=myapp

SOURCES= myclass.cpp \
         basis.cpp \
         integrator.cpp \
         lib.cpp \
         hfsolve.cpp \
         boysfunction.cpp \
         primitive.cpp \
         returnhermitecoeffs.cpp \
    kineticenergy.cpp

HEADERS= myclass.h \
         basis.h \
         integrator.h \
         primitive.h \
         lib.h \
         hfsolve.h \
         boysfunction.h \
         returnhermitecoeffs.h \
    kineticenergy.h

LIBS += -larmadillo -lblas -llapack
