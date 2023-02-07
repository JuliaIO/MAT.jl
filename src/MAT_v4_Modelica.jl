module MAT_v4_Modelica

import Base: read, write, close

mutable struct ModelicaV4
    ios::IOStream
    varnames::Dict{String, Int64}
    ModelicaV4(ios) = new(ios)
end


function read(mat::Matlabv4File)
    seekstart(mat.ios)


end

end #module