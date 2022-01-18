# A software for converting CAD models to hybrid grid

The program is used to convert CAD models(STEP and IGES files) to meshes and render using OpenGL.

## External Libraries

* [OpenGL](https://www.opengl.org/resources/)            Render
* [Qt](https://www.qt.io/)                               GUI
* [Eigen](http://eigen.tuxfamily.org/)                   Linear Algebra
* [OpenMesh](https://www.openmesh.org/)                  Mesh 
* [Open CADCADE](https://dev.opencascade.org/release)    Polyon Mesh Data Structure
* [tbb](https://github.com/oneapi-src/oneTBB)            Parallel Computation
* [triangle](http://www.cs.cmu.edu/~quake/triangle.html) Planar Triangle Mesh Generation

## Usage

```
git clone https://github.com/yanyiss/HybridGridProcessor
cd HybridGridProcessor
```

Edit lines 36-45 of CmakeLists.txt to set the values of library path
```
mkdir build && cd build
cmake -A x64 ..
```

Open **HybridGridProcessor.sln**, select **HybridGridProcessor** as launch project, and run.


## Supported File Formats

Mesh Files Format: .obj .off .ply .stl
CAD Files Format: .stp .STP .STEP .igs .IGS .IGES