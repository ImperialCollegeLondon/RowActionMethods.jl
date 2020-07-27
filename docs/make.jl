using Documenter, RowActionMethods

makedocs(;
    sitename="RowActionMethods.jl",
    modules=[RowActionMethods],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
        "Manual" => [
        ],
        "Reference" => [
            "Public API" => "ref/public.md",
            "Internal" => "ref/internal.md",
        ],
    ],
    repo="https://github.com/ImperialCollegeLondon/RowActionMethods.jl/blob/{commit}{path}#L{line}",
    authors="EdwardStables <edward.stables1198@gmail.com>",
    assets=String[],
)

deploydocs(;
    repo="github.com/ImperialCollegeLondon/RowActionMethods.jl",
)
