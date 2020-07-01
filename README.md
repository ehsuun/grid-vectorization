# Grid Vectorization

This repository contains the C++ implementation of the paper [Integer-Grid Sketch Simplification and Vectorization](https://repo-sam.inria.fr/d3/grid-vectorization/) (SGP 2020).

The code was developed on MacOS and tested on Fedora and Windows.

## Quick start
```
git clone --recursive https://gitlab.inria.fr/D3/grid-vectorization.git
cd grid-vectorization
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -GNinja
ninja
./IGSV_bin -i ../data/fish.png
```

## Usage

*Available soon.*

For now, just execute `./IGSV_bin -h` to show help.

## Dependencies (included)
- [alglib](https://www.alglib.net/)
- [CGAL](https://www.cgal.org/)
- [CoMISo](https://www.graphics.rwth-aachen.de/software/comiso/)
- [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
- [libigl](https://libigl.github.io/)
- [simple-svg](https://code.google.com/archive/p/simple-svg/)
- [Triangle](https://www.cs.cmu.edu/~quake/triangle.html)

We also provide an optional gui, which is built on top of [Dear ImGui](https://github.com/ocornut/imgui), [GLFW](https://www.glfw.org/) and [Glad](https://glad.dav1d.de/).

## License
GNU GPLv3
