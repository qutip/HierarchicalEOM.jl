push!(LOAD_PATH, "../src/")
using Documenter, Heom

PAGES = Any[
    "Home" => "index.md"
]

makedocs(
    sitename = "Heom.jl",
    pages=PAGES
)

deploydocs(
    repo="github.com/NCKU-QFort/Heom.jl.git",
)