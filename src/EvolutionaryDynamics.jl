
module EvolutionaryDynamics

include("Types.jl")
using .Types: FitnessReplication, NeutralReplication, WithoutReplication
using .Types: StandardMutation, RescaledMutation
using .Types: FitnessSelection, NeutralSelection, ElitismSelection

include("Methods.jl")
# using .Methods

include("TabularSystems.jl")
using .TabularSystems

# include("SATlandscapes.jl")
# using .SATlandscapes

# include("EvolvingSystems.jl")

# include("EDAnalysis.jl")

end
