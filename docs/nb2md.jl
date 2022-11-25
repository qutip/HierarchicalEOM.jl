# Convert Jupyter Notebooks (ipynb) into Markdown (md) files

# path for directories
NBdir = "src/notebooks/"
MDdir = "src/examples/"

function nb2md(file)
    NBpath = joinpath(NBdir, file)
    run(`jupyter nbconvert --to markdown --ExecutePreprocessor.kernel_name=julia-1.8 --ExecutePreprocessor.timeout=200 --output-dir=$MDdir --template=markdown_template.tpl --execute $NBpath `)
end

print("Deleting all files in src/examples/...")
rm("src/examples/", recursive=true)
mkdir("src/examples")
println("[DONE]")

files = filter(f -> endswith(f, ".ipynb"), readdir(NBdir))
for f in files
    nb2md(f)
end