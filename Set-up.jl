using Pkg
Pkg.activate(joinpath("C:\\Users", PC, "OneDrive\\PhD\\GitHub\\Networks"))
cd(joinpath("C:\\Users", PC, "OneDrive\\PhD\\GitHub\\Networks"))
meta_path = joinpath("C:\\Users", PC, "OneDrive\\PhD\\Metaweb Modelling")
simbio_path = joinpath("C:\\Users", PC, "OneDrive\\PhD\\GitHub\\simBio")
using ArchGDAL #, Shapefile, NCDatasets
using CSV, DataFrames
using NamedArrays, StaticArrays, OrderedCollections
using Rasters, RasterDataSources #, DimensionalData
using DynamicGrids, Dispersal
using Dates, Distributions, Serialization, StatsBase, JLD2
using ColorSchemes, Colors #Crayons, 
using Makie, WGLMakie, Plots # ImageMagick
# using OrdinaryDiffEq, EcologicalNetworksDynamics 
using DifferentialEquations
# using EcologicalNetworksDynamics
const DG, MK, AG, RS, Disp, DF, PL = DynamicGrids, Makie, ArchGDAL, Rasters, Dispersal, DataFrames, Plots
const COLORMAPS = [:magma, :viridis, :cividis, :inferno, :delta, :seaborn_icefire_gradient, :seaborn_rocket_gradient, :hot]