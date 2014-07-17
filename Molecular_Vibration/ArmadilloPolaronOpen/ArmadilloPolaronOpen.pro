TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
INCLUDEPATH += D:\DataNumerics\ArmadilloPolaronOpen\armaInclude

LIBS += -LD:\DataNumerics\ArmadilloPolaronOpen\armaInclude -lblas
LIBS += -LD:\DataNumerics\ArmadilloPolaronOpen\armaInclude -llapack
LIBS += -LD:\DataNumerics\ArmadilloPolaronOpen\armaInclude -llibf2c
SOURCES += main.cpp

QMAKE_CXXFLAGS += -openmp
