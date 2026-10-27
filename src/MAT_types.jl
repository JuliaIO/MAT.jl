# MAT_types.jl
# Internal types used by MAT.jl
#
# Copyright (C) 2012   Matthijs Cox
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

###########################################
## Reading and writing MATLAB .mat files ##
###########################################

module MAT_types

using StringEncodings: StringEncodings
import StringEncodings: Encoding
import Dates
import Dates: DateTime, Second, Millisecond
import PooledArrays: PooledArray, RefArray
using Tables: Tables
import OrderedCollections: OrderedDict
using StructUtils

export MatlabStructArray, StructArrayField, convert_struct_array
export MatlabClassObject
export MatlabOpaque, convert_opaque
export MatlabTable
export ScalarOrArray
export FunctionHandle
export MATStyle, MATReadStyle, MATWriteStyle
export MATConstructionError

include("types/utils.jl")
include("types/struct_array.jl")
# after struct_array.jl, whose type it dispatches on, and before the opaque
# converters, whose signatures use the `ScalarOrArray` it defines
include("types/structutils.jl")
include("types/function_handle.jl")
include("types/class_object.jl")
include("types/opaque.jl")
include("types/opaque_string.jl")
include("types/opaque_datetime.jl")
include("types/opaque_duration.jl")
include("types/opaque_categorical.jl")
include("types/table.jl")
include("types/char_array.jl")

end
