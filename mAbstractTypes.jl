
module mAbstractTypes

# ===================
# system types
abstract type atSystem end
abstract type atThermoSystem <: atSystem end

export atSystem, atThermoSystem

# ===================
# environment types
abstract type atEnv end
abstract type atCompEnv{T} <: atEnv end		# computational: input--ideal-output | T relates to the representation ...

export atEnv, atCompEnv

# ===================
# population dynamics types
abstract type atPopulation end
abstract type atEvotype end
abstract type atGenotype{T} end			# T relates to the representation of the genotypic variables
abstract type atPhenotype{T} end		# T relates to the representation of the phenotypic variables

abstract type at1dGty{T} <: atGenotype{T} end

export atPopulation, atEvotype, atGenotype, atPhenotype, at1dGty

end
