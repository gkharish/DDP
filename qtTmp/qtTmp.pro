TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    costfunction.cpp \
    dynamicmodel.cpp \
    ilqrsolver.cpp \
    romeosimpleactuator.cpp
SOURCES +=

include(deployment.pri)
qtcAddDeployment()

DISTFILES += \
    qtTmp.pro.user

HEADERS += \
    costfunction.h \
    dynamicmodel.h \
    ilqrsolver.h \
    romeosimpleactuator.h

INCLUDEPATH += /usr/include/eigen3
