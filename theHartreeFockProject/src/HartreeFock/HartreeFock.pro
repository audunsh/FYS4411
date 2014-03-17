include(../../defaults.pri)
TEMPLATE = app
SOURCES = main.cpp
#LIBS += -L../src/libs/ -lmyapp
LIBS += -L$$TOP_OUT_PWD/src/libs -lmyapp
