PC = "nicol"
include(joinpath("C:\\Users", PC, "OneDrive\\PhD\\GitHub\\Networks\\set-up.jl"))

using Plots
fw = Foodweb([1 => 2, 2 => 3]) # 1 eats 2, and 2 eats 3.