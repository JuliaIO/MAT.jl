# Support for reading and writing MATLAB files in Julia.

## Requirements

Until julia_hdf5 and MAT.jl are available as packages:

* MAT.jl must be cloned to a directory named "MAT" that resides somewhere in the load path.
* The "src" directory of [julia_hdf5](https://github.com/timholy/julia_hdf5/) must be in the load path.

## Quick start

To load the module:

```julia
include("MAT.jl")
using MAT
```

To read all variables from a MAT file:

```julia
file = matopen("matfile.mat")
read(file)
close(file)
```

To read a single variable from a MAT file:

```julia
file = matopen("matfile.mat")
read(file, "varname")
close(file)
```

To write a variable to a MAT file:

```julia
file = matopen("matfile.mat", "w")
write(file, "varname", variable)
close(file)
```

## Caveats

* All files are written in MATLAB v7.3 format.
* MATLAB v4 files are not currently supported.

## Credits

The MAT_HDF5 module, which provides read/write support for MATLAB v7.3 files, was written primarily by [Tim Holy](https://github.com/timholy/). The MAT_v5 module, which provides read support for MATLAB v5/v6/v7 files, was written primarily by [Simon Kornblith](https://github.com/simonster/).