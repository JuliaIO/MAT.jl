# Support for reading and writing MATLAB files in Julia.

## Installation

Within Julia, use the package manager:
```julia
Pkg.init()     # if you've never installed a package before
Pkg.add("MAT")
```

See also the requirements for the [HDF5](https://github.com/timholy/HDF5.jl/) module, used for "v7.3" files and for writing \*.mat files.

## Usage

To load the module:

```julia
using MAT
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

To read all variables from a MAT file as a Dict:

```julia
vars = matread("matfile.mat")
```

To write a Dict to a MAT file, using its keys as variable names:

```julia
matwrite("matfile.mat", {
	"myvar1" => 0
	"myvar2" => 1
})
```

## Caveats

* All files are written in MATLAB v7.3 format.
* MATLAB v4 files are not currently supported.

## Credits

The MAT_HDF5 module, which provides read/write support for MATLAB v7.3 files, was written primarily by [Tim Holy](https://github.com/timholy/). The MAT_v5 module, which provides read support for MATLAB v5/v6/v7 files, was written primarily by [Simon Kornblith](https://github.com/simonster/).
