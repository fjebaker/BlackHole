push!(LOAD_PATH,"src")

using Documenter, BlackHole

makedocs(
    modules=[BlackHole],
    clean=false,
    sitename="BlackHole Documentation",

    pages = [
        "Home" => "index.md",
        "Internal Methods" => "internals.md"
    ]
)

deploydocs(
    repo="github.com/fjebaker/BlackHole"
)