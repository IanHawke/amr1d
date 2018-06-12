# Adaptive Mesh Refinement

Berger-Oliger Adaptive Mesh Refinement is a nested grid technique for solving hyperbolic PDEs. This is an old toy code implementing AMR for one dimensional conservation laws.

## Usage

```
make amr
./amr <parameter_file>
```

Parameter files are Fortran namelists; have a look at the examples, and `parameters.f90` to see possibilities.

The only dependency should be a Fortran compiler: you will need to set the flags in the Makefile.

## Outputs

The output files are flat text and originally designed to be used with `gnuplot`. You should be able to plot them in standard tools, but it will be annoying.

## Notes

* This code is a decade old or so and chunks may not work.
* AMR needs an error estimate in order to place new grids. This code uses self-shadowing (see the work of Choptuik, Pretorius and co).
