include("EDutils.jl")

import EDutils

N = 2
β = 10
H = EDutils.inith(true, N)
τstps = 11
τarr = linspace(0.0, β, τstps)
szarr = []

# determine eigenvalues/vectors of H
eigv, U = eig(H)

# initialize density matrix
H = diagm(eigv, 0)
Z = trace(expm(-β*H))
ρ = expm(-β*H)/Z
print(ρ)
# generate sz operators
for i::Int=1:(N/2)
  push!(szarr, inv(U)*0.5*EDutils.initσz(N, i)*U)
end

sz_0 = szarr[1]
sz_expvals = []

for sz in szarr
  expval = []
  println(sz)
  for τ in τarr
    push!(expval, trace(ρ*expm(τ*H)*sz*expm(-τ*H)*sz_0))
  end
  push!(sz_expvals, expval)
end

println(sz_expvals)

