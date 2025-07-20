push!(LOAD_PATH,"../src/")
push!(LOAD_PATH, "../../CoordinateSystems.jl/src/")

using Documenter, CoordinateSystems, BlochSphereBigFloat

makedocs(
    sitename="BlochSphereBigFloat.jl",
    authors = "Lilith Emilia Höddinghaus",
    pages = [
        "index.md",
        "PulseSequenceRandomizedRotations.md"
    ],
)

# show documentation locally after generating with 
#   python3 -m http.server --bind localhost