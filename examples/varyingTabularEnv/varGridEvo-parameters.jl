# -*- coding: utf-8 -*-

# evolution and population constants
const RELNGENRELAX, RELNGENSAMPLE, NPOP = Int32(10^0), Int32(10^0), Int32(10^3)
const NSAMPLES = 10^3

# genotypic variables
GRIDSIZE = 12
DIMGSPACE = GRIDSIZE^2
BLOCKSIZE = GRIDSIZE÷3

# environmental constants
const REPFACTOR, REPCOEF = 0.0, 3.0
const NENV = 2
const LF, MF, HF = 1.0, 3.0, 5.0

aMutFactor = [ 0.02, 0.06, 0.2 ]
aSelStrength = [ 0.6, 1.0, 1.8 ]
aEnvTransRate = [ 0.001, 0.01,  0.1 ]

aFitnessTbl = [ [ LF for i in 1:DIMGSPACE ] for i in 1:NENV ]
for i in eachindex(aFitnessTbl[1])
    # top left block
    if (i-1)%GRIDSIZE < BLOCKSIZE && (i-1)÷GRIDSIZE < BLOCKSIZE
        aFitnessTbl[1][i] = HF
    end
end
for i in eachindex(aFitnessTbl[2])
    # bottom right block
    if (i-1)%GRIDSIZE >= GRIDSIZE - BLOCKSIZE && (i-1)÷GRIDSIZE >= GRIDSIZE - BLOCKSIZE
        aFitnessTbl[2][i] = HF
    end
end
