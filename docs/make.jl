using ElectricMachines
using Documenter

DocMeta.setdocmeta!(ElectricMachines, :DocTestSetup, :(using ElectricMachines); recursive=true)

makedocs(;
    modules=[ElectricMachines],
    authors="Nicolás Estay Loo",
    sitename="ElectricMachines.jl",
    format=Documenter.HTML(;
        canonical="https://ing-nestay.github.io/ElectricMachines.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ing-nestay/ElectricMachines.jl",
    devbranch="master",
)
