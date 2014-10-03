Polymer
=======


## Overview
This project, which concerns dynamics of exciton and polaron motion in conjugated polymer (currently one dimensional), is temporarily built on **C++** with **Qt** IDE. In order to implement efficient matrix manipulation, **LAPACK**, **BLAS**, etc with the interface of **Amardillo** is included.

## Windows Version [Click here to view project page](https://github.com/scott-chong-wong/Polymer)
* [**Must!** *C++ Compiler of VS2010 Express*](http://www.visualstudio.com/en-us/downloads#d-2010-express)
* [**Not necessary!** *Qt IDE*](http://qt-project.org/downloads)
* [**Already contained in Git Project!** *Armadillo with LAPACK, BLAS*](http://arma.sourceforge.net/)

## Linux Armadillo Page [Click here to view project page](https://github.com/scott-chong-wong/CppPrimer/tree/master/Armadillo_00)
# Armadillo for Ubuntu

## install package

* `sudo apt-get install cmake`

* `sudo apt-get install libopenblas-dev`

* `sudo apt-get install liblapack-dev`

* `sudo apt-get install libarpack-dev`

## install
It is recommended to read README.txt before installing.

* `cmake .`

* `make`

* `sudo make install`

## Invokation from C++ Source

* `#include <armadillo>`
* `using namespace arma;`
* `g++ -g -I/usr/include -c main.cpp`
* `g++ -g -L/usr/lib -o executable main.o -O1 -larmadillo`

