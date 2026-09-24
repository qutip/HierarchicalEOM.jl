const testdir = dirname(@__FILE__)

# Define the paths to the extension tests
const EXTENSION_PATH = Dict(
    "CUDA" => joinpath(testdir, "cuda"),
)
const EXTENSION_LIST = collect(keys(EXTENSION_PATH))

# the list of test groups
const GROUP_LIST = String[
    "All",
    "Main",
    "Code-Quality",
    EXTENSION_LIST...,
]

# prettier display of the list of test groups
const SHOW_GROUP_LIST = join("- " .* GROUP_LIST, "\n")

# if this file is executed directly, print the available test groups
if abspath(PROGRAM_FILE) == @__FILE__
    println("Available test GROUP:\n$(SHOW_GROUP_LIST)")
    nothing
end
