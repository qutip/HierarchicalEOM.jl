push!(LOAD_PATH, "../src/")
using Documenter, Heom

const PAGES = Any[
    "Home" => Any[
        "Introduction" => "index.md",
        "Installation" => "install.md",
        "QuickStart"   => "examples/quick_start.md",
        "Cite Heom"    => "cite.md"
    ],
    "Bosonic Bath" => "bosonic_bath.md",
    "Fermionic Bath" => "fermionic_bath.md",
    "Auxiliary Density Operators" => "ADOs.md",
    "HEOMLS Matrices" => Any[
        "Introduction" => "heom_matrix/intro.md",
        "HEOMLS for Schrödinger Equation" => "heom_matrix/schrodinger_eq.md",
        "HEOMLS for Bosonic Bath" => "heom_matrix/M_Boson.md",
        "HEOMLS for Fermionic Bath" => "heom_matrix/M_Fermion.md",
        "HEOMLS for Bosonic and Fermionic Bath" => "heom_matrix/M_Boson_Fermion.md",
        "HEOMLS for Master Equation" => "heom_matrix/master_eq.md",
    ],
    "Hierarchy Dictionary" => "hierarchy_dictionary.md",
    "Time Evolution" => "time_evolution.md",
    "Stationary State" => "stationary_state.md",
    "Spectrum" => "spectrum.md",
    # "Examples" => Any[],
    "Library" => Any[
        "Heom API" => "lib/heom_api.md",
        "Bath" => "lib/bath.md",
        "Bath Correlation Functions" => "lib/corr_func.md",
        "Physical Analysis Functions" => "lib/phys_analysis.md",
        "Misc." => "lib/misc.md"
    ]
]

makedocs(
    sitename = "Documentation | Heom.jl",
    pages=PAGES
)

deploydocs(
    repo="github.com/NCKU-QFort/Heom.jl.git",
)