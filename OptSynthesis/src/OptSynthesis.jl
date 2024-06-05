#__precompile__(false)   #uncomment if you want to compile the module
module OptSynthesis
using ModelingToolkit, MethodOfLines, OrdinaryDiffEq, DomainSets, Measures, DataInterpolations#, Plots 
using DiffEqParamEstim, Optimization, ForwardDiff, OptimizationNLopt, OptimizationOptimJL
using XLSX
using Serialization
using Plots#; plotly()
using DataFrames
using Metaheuristics


#= include loads the content of the respective files
into this module =#
include("DataHandling.jl")
include("ReactionsModels.jl") # most functions used in main.jl are in this file
include("OptOperation.jl")
include("plots.jl")

# Dict is a dictionary containing Pairs. Pairs consist of two variables: first => second
matdict = Dict( "ester" => 1, 
                "amine" => 2, 
                "TBD" => 3, 
                "product" => 4, 
                "intermediate" => 5, 
                "alcohol" => 6, 
                "4-(1-HE)-CYHEXOH" => 7, 
                "4-E-CYHEXOH" => 8,
                "sum" => 9,
                )


# maps the number of materials to variable indices of the ODE state --> we want to know the 
# concentration at the reactor outlet
# num_material => idx
# ATTENTION: This depends on the discretization step!!! (maybe should be assigned dynamically?)
outletidx = Dict( 
                    3 => [64,77,90],                         # 3 is the number of compounds
                    4 => [228,232,236,240],                  # 4 is the number of compounds
                    5 => [98,111,124,137,150],               # 5 is the number of compounds
                    6 => [115,128,141,154,167,180],          # 6 is the number of compounds
                )                 

# register functions 
@register q(t)
@register T(t)  
@register Cin(t)
@register C_ester_in(t)
@register C_amine_in(t)
@register C_TBD_in(t)
@register C_product_in(t)
@register C_intermediate_in(t)
@register C_alcohol_in(t)

export importexperiment, matdict, getreactionmodel, steadystate
export SynthesisModelSS
export conversion, selectivity, productivity, spacetimeyield
export plot_paracoords

end # module OptSynthesis
