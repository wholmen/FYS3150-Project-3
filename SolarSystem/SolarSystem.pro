TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
LIBS += -larmadillo -lblas -llapack

SOURCES += main.cpp \
    constants.cpp \
    solver.cpp \
    planet.cpp

HEADERS += \
    constants.h \
    solver.h \
    planet.h

