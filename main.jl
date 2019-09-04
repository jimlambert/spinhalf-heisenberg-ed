include("EDutils.jl")

using LinearAlgebra

N=4

lowLyingStates=EDutils.initMagBlock(0.0, 0.0, true,  N, 0)
diagLow = eigen(lowLyingStates[1])
println(diagLow.values[1])

