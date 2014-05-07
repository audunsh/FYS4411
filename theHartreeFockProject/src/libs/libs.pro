include(../../defaults.pri)
TEMPLATE = lib
TARGET=myapp

SOURCES= \
         basis.cpp \
         integrator.cpp \
         lib.cpp \
         hfsolve.cpp \
         boysfunction.cpp \
         primitive.cpp \
    contracted.cpp

HEADERS= \
         basis.h \
         integrator.h \
         primitive.h \
         lib.h \
         hfsolve.h \
         boysfunction.h \
    contracted.h

LIBS += -larmadillo -lblas -llapack
