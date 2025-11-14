execute = isempty(ARGS) || ARGS[1] == "run"

org, repo = :JuliaIO, :MAT
eval(:(using $repo))
using Documenter

# https://juliadocs.github.io/Documenter.jl/stable/man/syntax/#@example-block
ENV["GKSwstype"] = "100"
ENV["GKS_ENCODING"] = "utf-8"

isci = get(ENV, "CI", nothing) == "true"

format = Documenter.HTML(;
    prettyurls = isci,
    edit_link = "main",
    canonical = "https://$org.github.io/$repo.jl/stable/",
    assets = ["assets/custom.css"],
)

makedocs(;
    modules = [MAT],
    authors = "Contributors",
    sitename = "$repo.jl",
    format,
    pages = [
        "Home" => "index.md",
        "Object Arrays" => "object_arrays.md",
        "Methods" => "methods.md",
    ],
    warnonly = [:missing_docs,],
)

if isci
    deploydocs(;
        repo = "github.com/JuliaIO/MAT.jl",
        devbranch = "master",
        devurl = "dev",
        versions = ["stable" => "v^", "dev" => "dev"],
        forcepush = true,
        push_preview = true,
        # see https://$org.github.io/$repo.jl/previews/PR##
    )
end
