include("EDutils.jl")

using LinearAlgebra
#using PyPlot
#using Plots

# ------------------------------------------------------------------------------
# TYPES
# ------------------------------------------------------------------------------

ComplexOrFloat=Union{Complex{Real}, Float64}

# ------------------------------------------------------------------------------
# HELPER FUNCTIONS 
# ------------------------------------------------------------------------------

function make_magop(momentum::Float64, size::Int64, oparr)
  magmat=zeros(2^size, 2^size)
  for i=1:size
    magmat+=(oparr[i]*exp(im*momentum*i))/size
  end
  return magmat
end


#function compute_dyn(ψ::Array{ComplexOrFloat, 1}, evolutionop::Array{ComplexOrFloat, 2}, 
#                     operator::Array{ComplexOrFloat, 2}, time::Float64)
#  ψ_t=exp(-1im*time*evolutionop)*ψ
#  return adjoint(ψ)*operator*ψ
#end

function compute_dyn(ψ, evolutionop, operator, time)
  ψ_t=exp(-1im*time*evolutionop)*ψ
  return adjoint(ψ)*operator*ψ
end

# ------------------------------------------------------------------------------
# SETUP
# ------------------------------------------------------------------------------

N=2

qmin=0
qmax=pi
qstps=N

tmin=0
tmax=10
tstps=101

tarr=range(tmin, stop=tmax, length=tstps)
qarr=range(qmin, stop=qmax, length=qstps)

# HEISENBERG MODEL -------------------------------------------------------------
#initham=EDutils.inith(0.0, 1.0, true, N)
#evolham=EDutils.inith(0.0, 1.0, true, N)
# ------------------------------------------------------------------------------

# TFIM -------------------------------------------------------------------------
initham=EDutils.initTfim(true, N, 1.0, 0.0)
evolham=EDutils.initTfim(true, N, 1.0, 1.0)
# ------------------------------------------------------------------------------

# initial groundstate
initdiag=eigen(initham)
ψ_0=initdiag.vectors[1,:]

# ------------------------------------------------------------------------------
# INITIALIZE OPERATORS
# ------------------------------------------------------------------------------

# arrays with local Pauli matrices
σzarr=[EDutils.initσz(N, i) for i ∈ 1:N]
σxarr=[EDutils.initσx(N, i) for i ∈ 1:N]
σyarr=[EDutils.initσy(N, i) for i ∈ 1:N]

# arrays with (σ^α)_q operators 
σzqarr=[make_magop(q, N, σzarr) for q ∈ qarr]
σxqarr=[make_magop(q, N, σxarr) for q ∈ qarr]
σyqarr=[make_magop(q, N, σyarr) for q ∈ qarr]

# arrays with ((σ^α)_q)^2 operators
σzqsqarr=σzqarr.^2
σxqsqarr=σxqarr.^2
σyqsqarr=σyqarr.^2


# ------------------------------------------------------------------------------
# TIME EVOLVE
# ------------------------------------------------------------------------------

# first I want to compare how a linearized evolution compares to the full time
# evolution for various operators

σzqtre=zeros(tstps, qstps)
σzqtim=zeros(tstps, qstps)

σzqsqtre=zeros(tstps, qstps)
σzqsqtim=zeros(tstps, qstps)

# Measure echo 

projt=zeros(tstps)

for i = 1:tstps
  projt[i] = adjoint(ψ_0)*exp(-1im*tarr[i]*evolham)*ψ_0
end


#@time begin

for i = 1:tstps
  ψ_t=exp(-1im*tarr[i]*evolham)*ψ_0
  Threads.@threads for j = 1:qstps
    σzqtre[i,j] = real(adjoint(ψ_t)*σzqarr[j]*ψ_t)
    σzqtim[i,j] = imag(adjoint(ψ_t)*σzqarr[j]*ψ_t)
    σzqsqtre[i,j] = real(adjoint(ψ_t)*σzqsqarr[j]*ψ_t)
    σzqsqtim[i,j] = imag(adjoint(ψ_t)*σzqsqarr[j]*ψ_t)
  end
end

#Threads.@threads for i = 1:tstps
#  ψ_t=exp(-1im*tarr[i]*evolham)*ψ_0
#  #push!(magz_t, adjoint(ψ_t)*σzqarr[1]*ψ_t)
#  magztre[i] = real(adjoint(ψ_t)*σzqarr[1]*ψ_t)
#  magztim[i] = imag(adjoint(ψ_t)*σzqarr[1]*ψ_t)
#end

#for i = 1:tstps
#  ψ_t=exp(-1im*tarr[i]*evolham)*ψ_0
#  Threads.@threads for j = 1:qstps
#    σzqtre[i,j] = real(adjoint(ψ_t)*σzqarr[j]*ψ_t)
#    σzqtim[i,j] = imag(adjoint(ψ_t)*σzqarr[j]*ψ_t)
#    σzqsqtre[i,j] = real(adjoint(ψ_t)*σzqsqarr[j]*ψ_t)
#    σzqsqtim[i,j] = imag(adjoint(ψ_t)*σzqsqarr[j]*ψ_t)
#  end
#end

#end # @time

qfit = 4*(σzqsqtre-(σzqtre.^2))



for i = 1:qstps
  qcurr = fill(qarr[i], (1))
  plot(qarr, tarr, qfit[:,i])
end




#qfi_σzqtre = 4 .*(σzqsqtre-(σzqtre.^2))
#
#imshow(qfi_σzqtre, cmap="hot", interpolation="nearest", origin="lower",
#       extent=[0, pi, 0, tarr[tstps]], aspect=(pi/tarr[tstps]))
#axvline(pi/2)
#xlabel("momentum")
#ylabel("time")
#colorbar()
#show()

