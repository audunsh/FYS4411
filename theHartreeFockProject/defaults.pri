# Directories
INCLUDEPATH += $$TOP_PWD/src/libs
SRC_DIR = $$TOP_PWD

DEFINES += ARMA_USE_CXX11

release {
    QMAKE_CXXFLAGS_RELEASE -= -O2
    QMAKE_CXXFLAGS_RELEASE += -O3
    DEFINES += ARMA_NO_DEBUG       # faster code, but need to be sure that indexing is correct.
    #DEFINES += ARMA_EXTRA_DEBUG
    #DEFINES += ARMA_PRINT_ERRORS
}

debug {
    #DEFINES += ARMA_NO_DEBUG
    DEFINES += ARMA_EXTRA_DEBUG
    DEFINES += ARMA_PRINT_ERRORS
}

COMMON_CXXFLAGS = -std=c++0x
QMAKE_CXXFLAGS += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_RELEASE += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_DEBUG += $$COMMON_CXXFLAGS
