using MAT, Test

function check(filename, result)
    matfile = matopen(filename)
    for (k, v) in result
        @test haskey(matfile, k)
        got = read(matfile, k)
        if !isequal(got, v) || (typeof(got) != typeof(v) && (!isa(got, String) || !(isa(v, String))))
            close(matfile)
            error("""
                Data mismatch reading $k from $filename ($format)

                Got $(typeof(got)):

                $(repr(got))

                Expected $(typeof(v)):

                $(repr(v))
                """)
        end
    end
    @test union!(Set(), keys(matfile)) == union!(Set(), keys(result))
    close(matfile)

    mat = matread(filename)
    if !isequal(mat, result)
        error("""
            Data mismatch reading $filename ($format)

            Got:

            $(repr(mat))

            Expected:

            $(repr(result))
            """)
        close(matfile)
        return false
    end

    return true
end
cd(dirname(@__FILE__))
for filename in readdir("v4")
    #println("testing $filename")
    d = matread("v4/$filename")
    matwrite("v4/tmp.mat", d; version="v4")
    check("v4/tmp.mat", d)
    rm("v4/tmp.mat")
end
