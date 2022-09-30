push!(LOAD_PATH, "../src/")
using Documenter, Heom

const PAGES = Any[
    "Home" => Any[
        "Introduction" => "index.md",
        "Installation" => "install.md",
        "QuickStart"   => "quick_start.md",
        "Citation"     => "cite.md"
    ],
    "Bosonic Bath" => "bosonic_bath.md",
    "Fermionic Bath" => "fermionic_bath.md",
    "HEOM Matrix" => "heom_matrix.md",
    "Auxiliary Density Operators" => "ADOs.md",
    "Time Evolution" => "evolution.md",
    "Steady State" => "steady.md",
    "Spectrum" => "spectrum.md",
    "Examples" => Any[
        "JC" => "examples/JC_model.md"
    ],
    "Library" => Any[
        "Heom Base" => "lib/heom_base.md",
        "Bath" => "lib/bath.md",
        "Bath Correlation Functions" => "lib/corr_func.md",
        "Physical Analysis Functions" => "lib/phys_analysis.md"
    ]
]

makedocs(
    sitename = "Heom.jl",
    pages=PAGES
)

deploydocs(
    repo="github.com/NCKU-QFort/Heom.jl.git",
)