include (../../shared.pri)

HEADERS += \
    $$VCGDIR/vcg/complex/algorithms/clean.h\
    mesh_segmentation.h\
    Matrix.h

SOURCES += \
    mesh_segmentation.cpp\
    Matrix.cpp

TARGET = filter_mesh_segmentation


