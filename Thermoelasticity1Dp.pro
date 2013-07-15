TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt


QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS *= -fopenmp

LIBS += -lquadmath


SOURCES += main.cpp \
    tdma.cpp \
    imc_dfm.cpp \
    functor1d.cpp \
    diffusion1dp.cpp \
    boundary.cpp \
    thermoelasticity1dp.cpp

HEADERS += \
    tdma.h \
    imc_dfm.h \
    functor1d.h \
    diffusion1dp.h \
    boundary.h \
    thermoelasticity1dp.h

