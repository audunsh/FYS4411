TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -larmadillo -lblas -llapack

SOURCES += main.cpp \
hfsolve.cpp \
lib.cpp \

HEADERS += \
hfsolve.h \
testingHFSolve.h \
lib.h \

release {
    QMAKE_CXXFLAGS_RELEASE -= -O2
    QMAKE_CXXFLAGS_RELEASE += -O3
}

COMMON_CXXFLAGS = -std=c++0x
QMAKE_CXXFLAGS += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_RELEASE += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_DEBUG += $$COMMON_CXXFLAGS

#INCLUDEPATH += /home/goranbs/goran/CompPhys/programs/cppLibrary/
