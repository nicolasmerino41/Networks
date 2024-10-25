using Pkg
Pkg.activate(joinpath("C:\\Users", PC, "OneDrive\\PhD\\GitHub\\Networks"))
cd(joinpath("C:\\Users", PC, "OneDrive\\PhD\\GitHub\\Networks"))
meta_path = joinpath("C:\\Users", PC, "OneDrive\\PhD\\Metaweb Modelling")
using ArchGDAL #, Shapefile, NCDatasets
using CSV, DataFrames
using NamedArrays, StaticArrays, OrderedCollections
using Rasters, RasterDataSources #, DimensionalData
using DynamicGrids, Dispersal
using Dates, Distributions, Serialization, StatsBase, JLD2
using ColorSchemes, Colors #Crayons, 
using Makie, WGLMakie # ImageMagick, 
# using EcologicalNetworksDynamics
const DG, MK, AG, RS, Disp, DF = DynamicGrids, Makie, ArchGDAL, Rasters, Dispersal, DataFrames
const COLORMAPS = [:magma, :viridis, :cividis, :inferno, :delta, :seaborn_icefire_gradient, :seaborn_rocket_gradient, :hot]
