include("EDutils.jl")

using LinearAlgebra

# ------------------------------------------------------------------------------
# helper functions 
# ------------------------------------------------------------------------------

function make_magop(q::Float64, n::Int64, oparr)
  magmat=zeros(2^n, 2^n)
  for i=1:n
    magmat+=(oparr[i]*exp(im*q*i))/n
  end
  return magmat
end


# ------------------------------------------------------------------------------
N=4

qmin=0
qmax=pi
qstps=10

tmin=0
tmax=5
tstps=11

#lowLyingStates=EDutils.initMagBlock(0.0, 0.0, true,  N, 0)
#diagLow = eigen(lowLyingStates[1])
#println(diagLow.values[1])

# compute the initial ground state
initham=EDutils.inith(0.0, 0.0, true, N)
initdiag=eigen(initham)
initgs=initdiag.vectors[1]

# create hamiltonian that will generate dynamics
evolham=EDutils.inith(0.0, 1.0, true, N)

tarr=range(tmin, stop=tmax, length=tstps)
qarr=range(q)

σzops=[]
σxops=[]
σyops=[]

for i=1:N
  push!(σzops, EDutils.initσz(N, i))
  push!(σxops, EDutils.initσz(N, i))
  push!(σyops, EDutils.initσz(N, i))
end





