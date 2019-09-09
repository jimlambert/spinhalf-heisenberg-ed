include("EDutils.jl")

using LinearAlgebra

N=10

lowLyingStates=EDutils.initMagBlock(0.0, 1.0, true,  N, 0)
diagLow = eigen(lowLyingStates[1])
println(diagLow.values[1])

