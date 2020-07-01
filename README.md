Code for parametrization-based vectorization.

# IGSV -- Integer-grid sketch vectorization

## Dependencies
- Boost
- CGAL
- CoMISo
- libigl
- OpenCV
- Triangle
- (optional) OpenMP

MacOS: if cmake fails to find OpenMP, try executing the following before running cmake:

```bash
export LDFLAGS="-L/usr/local/opt/llvm/lib"
```

# OLD STUFF

## Build

```bash
cd PolyVectorization
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j12
```

First build is gonna take a while since `libigl` is gonna be built as a static lib.


## Run
The project has two main parts -- executables: parametrization and extraction.

The parametrization executable is used to load a line drawing, compute a frame field, and parametrize it.

```bash
# in build dir
./Param_gui -d cloud
```

The computed parametrization can be exported (`e` key) into the `serialized` folder in the root dir. The folder `serialized` is created automatically if it does not exist.

Extraction uses the computed parametrization to extract the final vectorization.

```bash
# in build dir, do either
./Extraction_gui ../serialized/cloud__2018-12-17__13h28m41s
# or simply
./Extraction_gui cloud__2018-12-17__13h28m41s
```

## Parameters

Key : **parameter** (*default value*) : what does it control?

### Sketch
- **binarization threshold** (*165*) : amount of binarization of the input.
- **stroke width** (*autocomputed*) : i. length of the cross for counting black pixels; ii. generation of the scalar field
- **dilation size** (*1*) : width of the cross for counting black pixels
- **tangent dir. threshold** (*0.45*) : threshold for the ratio of black pixel counts when deciding on the tangent direction

### Mesh
- **max area** (*10.0*) : maximal area of triangles around black pixels
- **sigma** (*0.002*) : variance of Gaussian weights from the distance transform
- **sigma inv** (*0.1*) : variance of Gaussian weights from the *inverse* distance transform

### Frame field
- **weight : smoothness** (*50.0*)
- **weight : alignment** (*1.0*)
- **weight : regularization** (*0.02*)

### Parametrization
- **min scale** (*10.0*) : lower bound of the absolute local scale.
- **max scale** (*40.0*) : upper bound of the absolute local scale.
- **weight : Poisson** (*1.0*)
- **weight : transition fn** (*10.0*)
- **weight : fitting** (*0.1*)
- **weight : alignment** (*10.0*)
- **CG iters** (*10*) : number of CG iterations to be performed at each optimization round
