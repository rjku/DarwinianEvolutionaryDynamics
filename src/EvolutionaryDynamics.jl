
module EvolutionaryDynamics

include("Types.jl")
using .Types: ReplicationType, FitnessReplication, NeutralReplication, WithoutReplication
using .Types: MutationType, StandardMutation, RescaledMutation
using .Types: SelectionType, FitnessSelection, NeutralReplication

include("Methods.jl")
# using .Methods

include("TabularSystems.jl")
using .TabularSystems

# include("EvolvingSystems.jl")

# include("EDAnalysis.jl")

end
