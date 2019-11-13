
using Revise, BenchmarkTools, PyPlot, MATLAB
using Statistics, LinearAlgebra
using mEvoTypes
import mEvoFunc, mUtils

# +
# data constants
const NEPOCHS, LOGFASMPL = 1, 1/5

# evolution and population constants
const NGEN, NPOP = Int32(10), Int32(1)
const REPFACTOR, PNTMUTFACTOR, NMUTMAX = 10.0, 0.07, Int32(15)

# metagenotypic constants
const SYSTEMSIZE, PVIABILITY, KOUT = Int32(6), 1.0, 100.0
const GMIN, GMAX = -1.5, 1.5

# environmental constants
const SELSTRENGTH, NIO = 1.0, Int32(4);
const HIGHFLOW, LOWFLOW = 1.0, 10.0^-4;

# +
# genotypic and metagenotypic types
aDisChnMGty = [ mEvoFunc.tDisChnMGty( SYSTEMSIZE, PVIABILITY, KOUT ) ]

# uniformly distributed initial genotypic variables
aDisChnGty = [ tCntGty( [aDisChnMGty[1]], GMIN .+ rand( aDisChnMGty[1].dG ) .* (GMAX - GMIN), [GMIN,GMAX] ) for i in 1:REPFACTOR*NPOP ];

# normally distributed initial genotypic variables
# aDisChnGty = [ tCntGty( [aDisChnMGty[1]], randn( aDisChnMGty[1].dG ), [GMIN,GMAX] ) for i in 1:REPFACTOR*NPOP ];

# niched array of genotypes
aDisChnGtyNiche = deepcopy(aDisChnGty[1:NPOP]);

# metagenotypes and genotypes from file
# aIsingMGty, aIsingGty, Npop = mEvoFunc.read_aIsingSigTransGty(isingDTMCprm,"test_#1")

# environmental and evotypic types: (DESIGNED/NONRANDOM) BOOLEAN AND-GATE IO
aIO = [
    [ [HIGHFLOW, HIGHFLOW, HIGHFLOW, HIGHFLOW, HIGHFLOW, HIGHFLOW], [HIGHFLOW, HIGHFLOW, LOWFLOW, LOWFLOW, LOWFLOW, LOWFLOW] ],
    [ [HIGHFLOW, HIGHFLOW, HIGHFLOW, HIGHFLOW, LOWFLOW, LOWFLOW], [LOWFLOW, LOWFLOW, LOWFLOW, LOWFLOW, HIGHFLOW, HIGHFLOW] ],
    [ [LOWFLOW, LOWFLOW, HIGHFLOW, HIGHFLOW, HIGHFLOW, HIGHFLOW], [LOWFLOW, LOWFLOW, LOWFLOW, LOWFLOW, HIGHFLOW, HIGHFLOW] ],
    [ [LOWFLOW, LOWFLOW, HIGHFLOW, HIGHFLOW, LOWFLOW, LOWFLOW], [LOWFLOW, LOWFLOW, LOWFLOW, LOWFLOW, HIGHFLOW, HIGHFLOW] ]
]
for io in aIO io[2] = io[2] ./ sum(io[2]) end

disChnEnv = tCompEnv(aIO, SELSTRENGTH)
disChnEty = tPntMutEty(REPFACTOR,PNTMUTFACTOR,[ aDisChnMGty[1].dG ],NMUTMAX);
disChnEtyNiche = tPntMutEty(10REPFACTOR,PNTMUTFACTOR,[ aDisChnMGty[1].dG ],NMUTMAX);

# print IOideal
for io in disChnEnv.IOidl
    println( replace(e -> e == 1.0 ? 1 : 0, io[1])," → ", replace( x -> x >= 1/SYSTEMSIZE ? 1 : 0, io[2]) )
end
# -

@time disChnPop = mEvoFunc.initLivingPop( NPOP,disChnEty,disChnEnv,aDisChnMGty,aDisChnGty );
@time disChnPopNiche = mEvoFunc.initLivingPop( NPOP,disChnEtyNiche,disChnEnv,aDisChnMGty,aDisChnGtyNiche );

aDisChnData = tEvoData[];
aDisChnDataNiche = tEvoData[];

for i in 1:NEPOCHS
	# push!( aDisChnData, tEvoData(NGEN, LOGFASMPL) )
	push!( aDisChnDataNiche, tEvoData(NGEN, LOGFASMPL) )
    # mEvoFunc.evolutionGKP!(disChnPop,aDisChnData[end],elite=false)
    mEvoFunc.evolutionOneNiches!(disChnPopNiche,aDisChnDataNiche[end])
end

plot( vcat([dataBtc.aveFitness for dataBtc in aDisChnDataNiche]...) ); gcf()
