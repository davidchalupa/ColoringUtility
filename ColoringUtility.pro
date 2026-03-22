#-------------------------------------------------
#
# Project created by QtCreator 2014-12-11T17:47:47
#
#-------------------------------------------------

QT       += core

TARGET = ColoringUtility
TEMPLATE = app
CONFIG += console

SOURCES += main.cpp\
    algorithm_greedyclique.cpp \
    algorithm_igcol.cpp \
    graphs.cpp \
    common.cpp \
    cli.cpp \
    statistics.cpp \
    random_generator.cpp \
    tabu_base.cpp \
    tabucol.cpp \
    algorithm.cpp \
    algorithm_brelaz.cpp \
    h2col.cpp

HEADERS  += graphs.h \
    algorithm_greedyclique.h \
    algorithm_igcol.h \
    common.h \
    cli.h \
    statistics.h \
    random_generator.h \
    tabu_base.h \
    tabucol.h \
    algorithm.h \
    algorithm_brelaz.h \
    h2col.h

FORMS    += mainwindow.ui

QMAKE_CFLAGS_RELEASE -= -O2
QMAKE_CFLAGS -= -O1
QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS -= -O1
QMAKE_CFLAGS  *= -O3
QMAKE_LFLAGS  *= -O3
QMAKE_CXXFLAGS *= -O3
