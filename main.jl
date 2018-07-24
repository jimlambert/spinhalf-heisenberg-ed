include("EDutils.jl")

import EDutils

N = 10
β = 100
H = EDutils.inith(true, N)
τstps = 501
τarr = linspace(0.0, β, τstps)
szarr = []

# determine eigenvalues/vectors of H

eigv, U = eig(H)

println("exact diagonalization complete")

# initialize density matrix
H = diagm(eigv, 0)
Z = trace(expm(-β*H))
ρ = expm(-β*H)/Z
# generate sz operators
for i::Int=1:(N/2)
  push!(szarr, inv(U)*0.5*EDutils.initσz(N, i)*U)
end

sz_0 = szarr[1]
sz_τexpvals = []
sz_texpvals = []
for sz in szarr
  τexpval = []
  texpval = []
  for τ in τarr
    push!(τexpval, trace(ρ*expm(τ*H)*sz*expm(-τ*H)*sz_0))
    push!(texpval, trace(ρ*expm(1im*τ*H)*sz*expm(-1im*τ*H)*sz_0))
  end
  push!(sz_τexpvals, τexpval)
  push!(sz_texpvals, texpval)
end

writedlm("./edtau.dat", sz_τexpvals)
writedlm("./edtre.dat", real(sz_texpvals))
writedlm("./edtim.dat", imag(sz_texpvals))
