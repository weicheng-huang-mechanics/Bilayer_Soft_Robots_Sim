# A Virtual Playground for the Bilayer Soft Robotics

### [Paper]() | [Video]()

## Overview

This study introduces a novel simulation framework based on the Discrete Elastic Rod (DER) model to accurately capture the dynamic behavior of bilayer soft robots, particularly in contact interactions. By leveraging discrete differential geometry, the approach enables efficient modeling of complex deformations, facilitating the design and control of advanced soft robotic systems.

<div align="center">
  <img src="assets/videos/demo.gif" alt="Bilayer robot">
</div>

## Prerequisites

- [Ubuntu 18.04 or above](https://ubuntu.com/tutorials/install-ubuntu-desktop#1-overview)
- C++ dependencies

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/DezhongT/Bilayer_Soft_Robots_Sim.git
   cd Bilayer_Soft_Robots_Sim
   ```
   
2. Install C++ dependencies

- **Note**: Some of these packages are installed to the system library for convenience. You may want to install locally to e.g., `~/.local` to avoid conflicts with system libraries. Add the `cmake` flag: `-D CMAKE_INSTALL_PREFIX=~/.local`. Then `sudo` is not required to install. You'll need to ensure subsequent builds know where to find the build libraries.

- [Eigen 3.4.0](http://eigen.tuxfamily.org/index.php?title=Main_Page)
  - Eigen is used for various linear algebra operations.
  - The project is built with Eigen version 3.4.0 which can be downloaded [here](https://gitlab.com/libeigen/eigen/-/releases/3.4.0). After downloading the source code, install through cmake as follows.
    ```bash
    cd eigen-3.4.0 && mkdir build && cd build
    cmake ..
    sudo make install
    ```

- [SymEngine](https://github.com/symengine/symengine)
  - SymEngine is used for symbolic differentiation and function generation.
  - Before installing SymEngine, LLVM is required which can be installed most easily via a package manager:
    - **Ubuntu**: `sudo apt-get install llvm`
  - Afterwards, install SymEngine from source using the following commands:
    ```bash
    git clone https://github.com/symengine/symengine
    cd symengine && mkdir build && cd build
    cmake -D WITH_LLVM=on -D BUILD_BENCHMARKS=off -D BUILD_TESTS=off ..
    make -j4
    sudo make install
    ```

- [Intel oneAPI Math Kernel Library (oneMKL)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html?operatingsystem=linux&distributions=webdownload&options=online)
  - Necessary for access to Pardiso, which is used as a sparse matrix solver.
  - Intel MKL is also used as the BLAS / LAPACK backend for Eigen.
  - **Ubuntu**: Follow the below steps.
    ```bash
    cd /tmp
    wget https://registrationcenter-download.intel.com/akdlm/irc_nas/18483/l_onemkl_p_2022.0.2.136.sh

    # This runs an installer, simply follow the instructions.
    sudo sh ./l_onemkl_p_2022.0.2.136.sh
    ```
  - Add one of the following to your .bashrc so that cmake can find the MKL library. Change the directory accordingly if your MKL version is different. 
   Note that older versions require setting `MKLROOT` while newer versions require `MKL_DIR`.
   You can find out which one from the cmake error message.
    ```bash
    export MKLROOT=/opt/intel/oneapi/mkl/2022.0.2   # for older versions
    export MKL_DIR=/opt/intel/oneapi/mkl/2024.2     # for newer versions
    ```

- [OpenGL / GLUT](https://www.opengl.org/)
  - OpenGL / GLUT is used for rendering the knot through a simple graphic.
  - Simply install through apt package manager:
    - **Ubuntu**: `sudo apt-get install libglu1-mesa-dev freeglut3-dev mesa-common-dev`

- [Lapack](https://www.netlib.org/lapack/) (*included in MKL*)

3. Configure the simulation engine
   ```bash
   mkdir build && cd build
   cmake ..
   make -j4
   cd ..
   ```

4. To simulate the bilayer robot with customized setting parameters, run
   ```bash
   ./simDER ./option.txt
   ```
   The parameters are specified in the ```option.txt``` with specifications as follows (SI units):
   - ```render (0 or 1) ```- Flag indicating whether OpenGL visualization should be rendered.
   - ```saveData (0 or 1)``` - Flag indicating whether positions should be recorded.
   - ```YoungM``` - Young's modulus.
   - ```totalTime``` - Total simulation time.
   - ```deltaTime``` - Time step size.
   - ```rodRadius``` - Cross-sectional radius of the beam.
   - ```density``` - Material density.
   - ```stol``` - A small number used in solving the linear system.
   - ```forceTol``` - Force tolerance.
   - ```maxIter``` - Maximum iteration.
   - ```viscosity``` - Viscosity.
   - ```scaleRendering``` - Dimension scale of rendering.
   - ```gVector``` - Gravitational vector.
   - ```Possion``` - Possion ratio.

   The simulation data of robot configuration will be saved to `datafiles/` directory, with each column in the file corresponding to `x`, `y`, `z`.

### Citation
If our work has helped your research, please cite the following paper.
```
@article{tong2024inverse,
  title={Inverse Design of Snap-Actuated Jumping Robots Powered by Mechanics-Aided Machine Learning},
  author={Tong, Dezhong and Hao, Zhuonan and Liu, Mingchao and Huang, Weicheng},
  journal={arXiv preprint arXiv:2408.10470},
  year={2024}
}

```
