cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

project(jspec2-python)

find_package(pybind11 REQUIRED)
find_package(GSL REQUIRED)

set(INCLUDE_DIR
    "${CMAKE_SOURCE_DIR}/include"
    )

set(SOURCE_FILES
    src/arbitrary_electron_beam.cc
    src/electron_beam.cc
    src/ecooling.cc
    src/force.cc
    src/functions.cc
    src/ibs.cc
    src/ion_beam.cc
    src/luminosity.cc
    src/other_effects.cc
    src/particle_model.cc
    src/ring.cc
    src/rms_dynamic.cc
    src/simulator.cc
    src/turn_by_turn.cc
)

set(HEADER_FILES
        ${INCLUDE_DIR}/jspec2/arbitrary_electron_beam.h
        ${INCLUDE_DIR}/jspec2/electron_beam.h
        ${INCLUDE_DIR}/jspec2/constants.h
        ${INCLUDE_DIR}/jspec2/cooler.h
        ${INCLUDE_DIR}/jspec2/datasink.h
        ${INCLUDE_DIR}/jspec2/ecooling.h
        ${INCLUDE_DIR}/jspec2/force.h
        ${INCLUDE_DIR}/jspec2/functions.h
        ${INCLUDE_DIR}/jspec2/ibs.h
        ${INCLUDE_DIR}/jspec2/ion_beam.h
        ${INCLUDE_DIR}/jspec2/luminosity.h
        ${INCLUDE_DIR}/jspec2/other_effects.h
        ${INCLUDE_DIR}/jspec2/particle_model.h
        ${INCLUDE_DIR}/jspec2/ring.h
        ${INCLUDE_DIR}/jspec2/rms_dynamic.h
        ${INCLUDE_DIR}/jspec2/simulator.h
        ${INCLUDE_DIR}/jspec2/turn_by_turn.h
        ${INCLUDE_DIR}/jspec2/twiss.h
        ${INCLUDE_DIR}/jspec2/rate.h
)

set(PYTHON_FILES
        src/electron_beam_pb.cc
        src/cooler_pb.cc
	src/datasink_pb.cc
        src/ecooling_pb.cc
        src/force_pb.cc
        src/ibs_pb.cc
        src/ion_beam_pb.cc
        src/jspec2_pb.cc
        src/luminosity_pb.cc
        src/ring_pb.cc
        src/simulators_pb.cc
)

pybind11_add_module(jspec
	${SOURCE_FILES}
	${HEADER_FILES}
	${PYTHON_FILES}
)
target_include_directories(
    jspec PUBLIC
    ${INCLUDE_DIR}
    ${GSL_INCLUDE_DIRS}
)

target_link_libraries(jspec PUBLIC ${GSL_LIBRARIES})

