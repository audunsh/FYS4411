TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    HartreeFock.cpp \
    ../../programs/cppLibrary/lib.cpp

OTHER_FILES += \
    ../../programs/cppLibrary/lib.i

HEADERS += \
    ../../programs/cppLibrary/lib.h

