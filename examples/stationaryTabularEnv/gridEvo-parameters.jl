# -*- coding: utf-8 -*-
# evolution and population constants
const NGENRELAX, NGENSAMPLE, NPOP = Int32(10^4), Int32(10^4), Int32(10^3)
const NSAMPLES = 10^4


# genotypic variables
const GRIDSIZE = 12
const DIMGSPACE = GRIDSIZE^2
const BLOCKSIZE = GRIDSIZE÷3

# environmental constants
const REPFACTOR, REPCOEF = 0.0, 4.0
const LF, MF, HF = 1.0, 3.0, 5.0

aMutFactor = [ 0.02, 0.06, 0.2 ]
aSelStrength = [ 0.6, 1.0, 1.8 ]
# aMutFactor = [ 0.03 ]
# aSelStrength = [ 1.0 ]

fitnessEstimation = exp(aSelStrength[1]*LF), exp(aSelStrength[1]*HF)

fTbl = [ LF for i in 1:DIMGSPACE ]
for i in eachindex(fTbl)
    # top left block
    if (i-1)÷GRIDSIZE < BLOCKSIZE && (i-1)%GRIDSIZE < BLOCKSIZE
        fTbl[i] = HF
    end

    # bottom right line
    if ( (i-1)÷GRIDSIZE == GRIDSIZE - BLOCKSIZE && (i-1)%GRIDSIZE > GRIDSIZE - BLOCKSIZE ) ||
        ( (i-1)÷GRIDSIZE >= GRIDSIZE - BLOCKSIZE && (i-1)%GRIDSIZE == GRIDSIZE - BLOCKSIZE )
        fTbl[i] = HF
    end

    # bottom left point
    if (i-1)÷GRIDSIZE == BLOCKSIZE - 2 && (i-1)%GRIDSIZE == GRIDSIZE - BLOCKSIZE + 1
        fTbl[i] = HF
    end

    # top right point
    if (i-1)÷GRIDSIZE == GRIDSIZE - BLOCKSIZE + 1 && (i-1)%GRIDSIZE == BLOCKSIZE - 2
        fTbl[i] = HF
    end

    # top right point
#     if GRIDSIZE - BLOCKSIZE <= (i-1)÷GRIDSIZE <= GRIDSIZE - BLOCKSIZE + 1 && (i-1)%GRIDSIZE == BLOCKSIZE - 1
#         f[i] = HF
#     end

    # top right block
#     if (i-1)÷GRIDSIZE >= GRIDSIZE - BLOCKSIZE && (i-1)%GRIDSIZE < BLOCKSIZE
#         fTbl[i] = MF
#     end
end
