# Stochastic simulations

This directory is the code needed to run the stochastic simulations, and the data from it.

Running the stochastic simulations requires this [simple C++ library for the Gillespie algorithm](https://github.com/smrfeld/gillespie-cpp), which you must install separately.

<img src="figures_readme/trajs.png" alt="drawing" width="600"/>

## Load data from DVC

To just load the data from DVC, run `dvc pull`. Otherwise you must run the stochastic simulations from scratch - this can take a while!

## Running stochastic simulations

**Running the stochastic simulations requires this [simple C++ library for the Gillespie algorithm](https://github.com/smrfeld/gillespie-cpp), which you must install separately.**

The source files are in the `src` directory. They are built with `CMake` from a dedicated build directory:
```
mkdir build
cd build
cmake ..
```
or, you can use your favorite generator like `cmake .. -GXcode`.

The binaries are in the bin directory:
```
cd bin
./run_gillespie
```

The default output directories are the `data_gillespie` and `data_tau_leaping` folders. You can see how to plot the trajectories using the Mathematica notebooks in the [mathematica](mathematica) directory.