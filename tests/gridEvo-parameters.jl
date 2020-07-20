
# evolution and population constants
const RELNGENRELAX, RELNGENSAMPLE, NPOP = Int32(10^2), Int32(2*10^2), Int32(10^3)
const NSAMPLES = 10^3

# genotypic variables
GRIDSIZE = 8
DIMGSPACE = GRIDSIZE^2
BLOCKSIZE = GRIDSIZE÷3

# environmental constants
const REPSTRENGTH, REPFACTOR, MINREPCOEF = 0.0, 0.0, 3.0
const NENV = 3
const LF, MF, HF = 1.0, 3.0, 5.0

aMutFactor = [ 0.01, 0.03, 0.1 ]
aSelStrength = [ 0.6, 1.0, 1.8 ]
aEnvTransRate = [ 0.001, 0.01,  0.1 ]
# aλM = [ 1.0, 3.0, 10.0 ]
aλM = [ 1.0 ]

aW = Vector{Matrix{Float64}}(undef,NENV)
for (i,w) in enumerate(aEnvTransRate)
    aW[i] = [ [ 1.0 - w, w, 0.0 ] [ 0.0, 1.0 - w, w ] [ w, 0.0, 1.0 - w ] ]
end

afTb = [ [ LF for i in 1:DIMGSPACE ] for i in 1:NENV ]
for i in eachindex(afTb[1])
    # top left block
    if (i-1)÷GRIDSIZE < BLOCKSIZE && (i-1)%GRIDSIZE < BLOCKSIZE
        afTb[1][i] = HF
    end
end
for i in eachindex(afTb[2])
    # top right block
    if (i-1)÷GRIDSIZE >= GRIDSIZE - BLOCKSIZE && (i-1)%GRIDSIZE < BLOCKSIZE
        afTb[2][i] = HF
    end
end
for i in eachindex(afTb[3])
    # top bottom center block
    if GRIDSIZE÷2 - 1 <= (i-1)÷GRIDSIZE <= GRIDSIZE÷2 && (i-1)%GRIDSIZE > GRIDSIZE - BLOCKSIZE - 1
        afTb[3][i] = HF
    end
end
