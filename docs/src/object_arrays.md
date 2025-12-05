# Objects and struct arrays

To better handle special cases we have these types since MAT 0.11:
* [`MatlabStructArray`](@ref MatlabStructArray)
* [`MatlabClassObject`](@ref MatlabClassObject)

## Struct arrays vs Cell arrays

Cell arrays are written for `Array{Any}` or any other unsupported element type:

```julia
sarr = Any[
    Dict("x"=>1.0, "y"=>2.0),
    Dict("x"=>3.0, "y"=>4.0)
]
matwrite("matfile.mat", Dict("cell" => sarr))

```

Inside MATLAB you will find:

```matlab
>> load('matfile.mat')
>> cell

cell =

  2×1 cell array

    {1×1 struct}
    {1×1 struct}
```

Read and write behavior for struct arrays is different. For struct arrays we use the `MatlabStructArray` type. You can also write with MAT.jl using Dict arrays `AbstractArray{<:AbstractDict}` if all the Dicts have equal keys, which will automatically convert internally to `MatlabStructArray`.

```julia
sarr = Dict{String, Any}[
    Dict("x"=>1.0, "y"=>2.0),
    Dict("x"=>3.0, "y"=>4.0)
]
matwrite("matfile.mat", Dict("s" => sarr))
# which is the same as:
matwrite("matfile.mat", Dict("s" => MatlabStructArray(sarr)))
# which is the same as:
matwrite("matfile.mat", Dict("s" => MatlabStructArray(["x", "y"], [[1.0, 3.0], [2.0, 4.0]])))
```

Now you'll find the following inside MATLAB:

```matlab
>> load('matfile.mat')
>> s

s =

[2x1 struct, 576 bytes]
x: 1
y: 2
```

Note that when you read the file again, you'll find the `MatlabStructArray`, which you can convert back to the Dict array with `Array`:

```julia
julia> sarr = matread("matfile.mat")["struct_array"]
MatlabStructArray{1} with 2 columns:
 "x": Any[1.0, 3.0]
 "y": Any[2.0, 4.0]

julia> sarr["x"]
2-element Vector{Any}:
 1.0
 3.0

julia> Array(sarr)
2-element Vector{Dict{String, Any}}:
 Dict("x" => 1.0, "y" => 2.0)
 Dict("x" => 3.0, "y" => 4.0)

```

Note that before v0.11 MAT.jl will read struct arrays as a Dict with concatenated arrays in the fields/keys, which is equal to `Dict(sarr)`.

## Object Arrays

You can write an old class object with the `MatlabClassObject` and arrays of objects with `MatlabStructArray` by providing the class name. These are also the types you obtain when you read files.

> Please note that the order of the fields is important for MatlabClassObjects to be read properly in MATLAB. You may get `Warning: Fields of object 'tc' do not match the current constructor definition for class 'TestClass'. The object has been converted to a structure.`. Consider using `OrderedDict` if you have multiple fields.

Write a single class object:
```julia
d = Dict{String,Any}("foo" => 5.0)
obj = MatlabClassObject(d, "TestClass")
matwrite("matfile.mat", Dict("tc" => obj))
```

A class object array
```julia
class_array = MatlabStructArray(["foo"], [[5.0, "bar"]], "TestClass")
matwrite("matfile.mat", Dict("class_array" => class_array))
```

Also a class object array, but will be converted to `MatlabStructArray` internally:
```julia
class_array = MatlabClassObject[
    MatlabClassObject(Dict{String,Any}("foo" => 5.0), "TestClass"),
    MatlabClassObject(Dict{String,Any}("foo" => "bar"), "TestClass")
]
matwrite("matfile.mat", Dict("class_array" => class_array))
```

A cell array:
```julia
cell_array = Any[
    MatlabClassObject(Dict{String,Any}("foo" => 5.0), "TestClass"),
    MatlabClassObject(Dict{String,Any}("a" => "bar"), "AnotherClass")
]
matwrite("matfile.mat", Dict("cell_array" => cell_array))
```
