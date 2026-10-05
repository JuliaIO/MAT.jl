using SnoopCompileCore: @snoop_invalidations

invalidations = @snoop_invalidations begin
    import MAT
end

using SnoopCompile: SnoopCompile, filtermod, invalidation_trees, uinvalidated
inv_owned = length(filtermod(MAT, invalidation_trees(invalidations)))
inv_total = length(uinvalidated(invalidations))
inv_deps = inv_total - inv_owned
@show inv_total, inv_deps

import PrettyTables
SnoopCompile.report_invalidations(;
    invalidations,
    process_filename = x -> last(split(x, ".julia/packages/")),
    n_rows = 0,
)
for methodinstance in sort!(string.(collect(uinvalidated(invalidations))))
    println(methodinstance)
end

using Printf
open(ENV["GITHUB_OUTPUT"], "a") do io
    println(io, @sprintf("total=%09d", inv_total))
    println(io, @sprintf("deps=%09d", inv_deps))
end
