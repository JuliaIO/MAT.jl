# Support for reading and writing MATLAB files in Julia.

[![Build Status](https://travis-ci.org/JuliaIO/MAT.jl.svg?branch=master)](https://travis-ci.org/JuliaIO/MAT.jl)

This library can read MATLAB `.mat` files, both in the older v5/v6/v7 format, as well as the newer v7.3 format.

## Installation

This is installed using the standard tools of the [package manager](https://julialang.github.io/Pkg.jl/v1/getting-started/):

```julia
pkg> add MAT
```
where you get the `pkg>` prompt by hitting `]` as the first character of the line. (Exit `pkg` mode by hitting backspace or Ctrl-C as the first character of the line.)

See also the requirements for the [HDF5](https://github.com/timholy/HDF5.jl/) module, used for "v7.3" files and for writing \*.mat files.

## Usage

To load the module:

```julia
using MAT
```

To read a single variable from a MAT file (compressed files are detected and handled automatically):

```julia
file = matopen("matfile.mat")
read(file, "varname") # note that this does NOT introduce a variable ``varname`` into scope
close(file)
```

To write a variable to a MAT file:

```julia
file = matopen("matfile.mat", "w")
write(file, "varname", variable)
close(file)
```

To read all variables from a MAT file as a Dict:

```julia
vars = matread("matfile.mat")
```

To write a Dict to a MAT file, using its keys as variable names.
The `compress` argument is optional, and compression is off by default:

```julia
matwrite("matfile.mat", Dict(
	"myvar1" => 0,
	"myvar2" => 1
); compress = true)
```

To get a list of variable names in a MAT file:

```julia
file = matopen("matfile.mat")
varnames = names(file)
close(file)
```

To check for the presence of a variable name in a MAT file:

```julia
file = matopen("matfile.mat")
if exists(file, "variable")
    # something
end
close(file)
```

## Caveats

* All files are written in MATLAB v7.3 format.
* MATLAB v4 files are not currently supported.

## Credits

The MAT_HDF5 module, which provides read/write support for MATLAB v7.3 files, was written primarily by [Tim Holy](https://github.com/timholy/). The MAT_v5 module, which provides read support for MATLAB v5/v6/v7 files, was written primarily by [Simon Kornblith](https://github.com/simonster/).
