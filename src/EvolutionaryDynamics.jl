
module EvolutionaryDynamics

include("Types.jl")
# using .Types

include("Methods.jl")
# using .Methods

include("TabularSystems.jl")
using .TabularSystems

#: TrajectoryData, VaryingEnvironmentTrajectoryData, generateTabularSystemsTrajectories, generateVaryingTabularSystemsITrajectories, IsVaryingEnvironment
# export TabularEvotype, TabularEnvironment, init_Population
# export TrajectoryData, VaryingEnvironmentTrajectoryData, generateTabularSystemsTrajectories, generateVaryingTabularSystemsITrajectories, IsVaryingEnvironment

# include("EvolvingSystems.jl")

# include("EDAnalysis.jl")

end
