cmake_minimum_required(VERSION 3.11)
project(soft_robot)

set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_BUILD_TYPE Release)
# set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_FLAGS "-O3 -DNDEBUG -Wno-deprecated-declarations")
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../")

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

set(THREADS_PREFER_PTHREAD_FLAG ON)

# Per https://eigen.tuxfamily.org/dox/TopicUsingIntelMKL.html eigen should be used
# with MKL LAPACK only with lp64 set. The default may be ILP64, which we don't want.
set(MKL_INTERFACE_FULL intel_lp64)

add_definitions(-DEIGEN_USE_MKL_ALL)

find_package(SymEngine REQUIRED)
find_package(GLUT REQUIRED)
find_package(OpenGL REQUIRED)
find_package(MKL CONFIG REQUIRED)
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)
find_package(Threads REQUIRED)
find_package(Eigen3 3.4 REQUIRED)
find_package(OpenMP)

if (${HAVE_SYMENGINE_LLVM})
    message(STATUS "Testing SymEngine LLVM & SBML support - found")
else ()
    message(WARNING "SymEngine library is missing LLVM and/or SBML support")
    message(WARNING "${SYMENGINE_LLVM_ERROR_LOG}")
endif ()

include_directories(${SYMENGINE_INCLUDE_DIRS}
        ${EIGEN_INCLUDE_DIRS}
        ${OPENGL_INCLUDE_DIRS}
        ${GLUT_INCLUDE_DIRS}
        MKL::MKL)

link_directories(${MKL_LIBRARY_DIR})


# Run with cmake --fresh -DCREATE_DOCS=ON .. only when documentation needs to be updated.

# If docs change and you want to preview the updates in index.html, you'll need to rm -rf * and
# remake everything in build. Don't make clean && make. Kind of annoying but this should only be
# done for previewing. All updates to the actual docs will be done more fluidly.

# If wanting to rebuild source code, will have to cmake --fresh .. again.

# Probably a good idea to split the documentation cmake process into its own directory later.
option(CREATE_DOCS OFF)

if(CREATE_DOCS)
    add_subdirectory("docs")
else()
    add_executable(simDER
        src/main.cpp
        
        src/timeStepper.cpp
        src/world.cpp

        src/setInput.cpp

        src/elasticPlate.cpp

        src/mechanics/base_force.cpp

        src/mechanics/forces/inertialForce.cpp
        src/mechanics/forces/dampingForce.cpp        
        src/mechanics/forces/elasticBendingForce.cpp
        src/mechanics/forces/elasticStretchingForce.cpp
        src/mechanics/forces/elasticTwistingForce.cpp
        src/mechanics/forces/externalGravityForce.cpp

        src/mechanics/constraints/elasticTwistingBound.cpp
        src/mechanics/constraints/elasticBendingBound.cpp 
        src/mechanics/constraints/symbolicEquations.cpp
        src/mechanics/constraints/constrainedForce.cpp
    )

    target_include_directories(simDER PRIVATE
            src
    )

    target_link_libraries(simDER PRIVATE
            ${OPENGL_LIBRARIES}
            ${GLUT_LIBRARIES}
            MKL::MKL
            Eigen3::Eigen
            ${LAPACK_LIBRARIES}
            ${BLAS_LIBRARIES}

            lapack

            # gfortran
            Threads::Threads
            m
            ${SYMENGINE_LIBRARIES}
    )

    if (OpenMP_CXX_FOUND)
        target_link_libraries(simDER PUBLIC OpenMP::OpenMP_CXX)
    endif ()
endif()

