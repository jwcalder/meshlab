# Copyright 2019 Collabora, Ltd.
# SPDX-License-Identifier: BSL-1.0

### Generated file! Edit the templates in src/templates then re-run ./make-cmake.py

set(SOURCES mesh_segmentation.cpp Matrix.cpp)

set(HEADERS mesh_segmentation.h Matrix.h)

#set(RESOURCES meshlab.qrc)

add_library(filter_mesh_segmentation MODULE ${SOURCES} ${HEADERS} ${RESOURCES})

target_include_directories(filter_mesh_segmentation PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(filter_mesh_segmentation PUBLIC common)

set_property(TARGET filter_mesh_segmentation PROPERTY FOLDER Plugins)

set_property(TARGET filter_mesh_segmentation PROPERTY RUNTIME_OUTPUT_DIRECTORY
                                           ${MESHLAB_PLUGIN_OUTPUT_DIR})

set_property(TARGET filter_mesh_segmentation PROPERTY LIBRARY_OUTPUT_DIRECTORY
                                           ${MESHLAB_PLUGIN_OUTPUT_DIR})

install(TARGETS filter_mesh_segmentation DESTINATION ${MESHLAB_PLUGIN_INSTALL_DIR}
                COMPONENT Plugins)
