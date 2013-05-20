# MAT_macros.jl
# Macros for emulating MATLAB-like saving and loading functionality
#
# Copyright (C) 2013   Jon Malmaud
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

module MAT_macros

export @save, @load

macro save(filename, vars...)    
    filename=ensure_mat(filename)
    esc(quote
        let d=Dict{String,Any}()
            for var in $(vars)
                d[string(var)] = eval(var)
            end
            MAT.matwrite($(filename), d)
        end
    end)
end


macro load(filename)        
    filename=ensure_mat(filename)    
    esc(quote
        let k, v
            for (k,v) in MAT.matread($(filename))
                eval(:($(symbol(k))=$(v)))
            end
        end
    end)
end


function ensure_mat(filename)
    filename=string(filename)
    if !ismatch(r"\.mat$", filename)
        filename=string(filename,".mat")
    end
    filename
end

end
