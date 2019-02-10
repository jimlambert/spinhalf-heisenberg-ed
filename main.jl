include("EDutils.jl")

using LinearAlgebra

N = 2
q = pi
β = 10
tmin = 0
tmax = 5
tstps = 11
genH = EDutils.inith(true, N) # generator Hamiltonian
evoH = EDutils.inith(-0.5, true, N) # evolution Hamiltonian
tarr = range(tmin, stop=tmax, length=tstps)
szarr = [] # array of sz operators in eigenbasis
sxarr = [] # array of sx operators in eigenbasis
amplitude = []
szvals = []
sz2vals = []
sxvals = []
sx2vals = []
# determine eigenvalues/vectors of H
F = eigen(genH)
U = F.vectors
gstate = U[:,1]

println("exact diagonalization complete")


# generate sz operators
for i::Int=1:(N/2)
  push!(szarr, inv(U)*0.5*EDutils.initσz(N, i)*U)
  push!(sxarr, inv(U)*0.5*EDutils.initσx(N, i)*U)
end

sztot = zeros(2^N, 2^N)
println(szarr)
for (r,sz) ∈ enumerate(szarr)
  global sztot += exp(1im*q*r)*sz
end
sztot = sztot / N
sztot2 = sztot*sztot


sxtot = zeros(2^N, 2^N)
for (r,sx) ∈ enumerate(sxarr)
  global sxtot += exp(1im*q*r)*sx
end
sxtot = sxtot / N
sxtot2=sxtot*sxtot

for t in tarr 
  global amplitude
  global sztot2
  #global sxtot2
  global sztot
  #global sxtot
  gstatenew = exp(-1im*t*evoH)*gstate 
  push!(amplitude, abs(adjoint(gstate)*gstatenew))
  global gstate = gstatenew
  #push!(sztot2, adjoint(gstate)*sztot2*gstate)
  #push!(sxtot2, adjoint(gstate)*sxtot2*gstate)
  #push!(sztot, adjoint(gstate)*sztot*gstate)
  #push!(sxtot, adjoint(gstate)*sxtot*gstate)
end
println(amplitude)
# sz_0 = szarr[1]
# sz_τexpvals = []
# sz_texpvals = []
# for sz in szarr
#   texpval = []
#   for τ in τarr
#     push!(texpval, trace(ρ*expm((-0.1*β+1im*τ)*H)*sz*expm((0.1*β-1im*τ)*H)*sz_0))
#   end
#   push!(sz_texpvals, texpval)
# end

#writedlm("./edtre-mixt.dat", real(sz_texpvals))
#writedlm("./edtim-mixt.dat", imag(sz_texpvals))
