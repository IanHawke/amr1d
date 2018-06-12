# Adaptive Mesh Refinement

Berger-Oliger Adaptive Mesh Refinement is a nested grid technique for solving hyperbolic PDEs. This is an old toy code implementing AMR for one dimensional conservation laws.

## Usage

```
make amr
./amr
```

Right now the code is reading its parameters from `input.par`, which seems to be hard-coded.

Parameter files are Fortran namelists; have a look at the examples, and `parameters.f90` to see possibilities.

The only dependency should be a Fortran compiler: you will need to set the flags in the Makefile.

## Outputs

The output files are flat text and originally designed to be used with `gnuplot`. You should be able to plot them in standard tools, but it will be annoying. With Python you could try something like (assuming standard `numpy` and `matplotlib` imports)

```
from glob import glob
files = glob('amr_l*dat')
grids = [np.genfromtext(f) for f in files]
for g in grids:
    plt.plot(g[:,0], g[:,1], 'x')
plt.show()
```

## Notes

* This code is a decade old or so and chunks may not work.
* AMR needs an error estimate in order to place new grids. This code uses self-shadowing (see the work of Choptuik, Pretorius and co).
* This was written in part as an exercise in Fortran OOP programming, in part to test some ideas ahead of building a different code, and in part to get some highly accurate 1d results. It has not been updated for a long time, and will not be updated in future.
