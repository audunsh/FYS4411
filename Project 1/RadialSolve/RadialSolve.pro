TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += /usr/local/Cellar/armadillo/4.000.0/include

SOURCES += main.cpp
SOURCES += ../lib.cpp
LIBS += -llapack -armadillo
