module EDUtils

export flip!
export getid
export getstate


"""
# flip!(Int, Int, Int)

## Arguments

  * s::Int   -   integer corresponding to spins state
  * i::Int   -   index of spin state
  * j::Int   -   index of spin state not equal to i

## Returns

  * b       -   integer id corresponding to flipped state

## Description

Function with one method. Accepts three integers, s, i, j and flips the ith and
jth bits of the binary representation of the integer s. The integer s is
modified in place.  
"""
function flip(s::Int, i::Int, j::Int)
  f::Int64 = 2^(i-1) + 2^(j-1)
  return xor(s-1, f) + 1
end


"""
# getid(Array{Int64, 1})

## Arguments

* state::Array{Int64, 1}   -   array of containing spin state with +1/2 (1) and 
                               -1/2 (0) spins

## Returns

  * id::Int                -  the id number for the input state

## Description

Function with one method. Accepts an array of integers, state. The state of a
spin chain is mapped on to 1's or 0's to ease conversion to and from the state
ids. 

"""
function getid(state::Array{Int64, 1})
  sum = 0
  len = length(state)
  for (site, spin) in enumerate(state) 
    sum += spin*2^(site-1) 
  end
  return sum+1
end


"""
# spin(Int, Int, Int)  

## Arguments

  * id::Int    -   state id
  * n::Int     -   system size
  * i::Int     -   spin of interest

## Returns

  * state[i]::Int   -  the spin at position i of the state corresponding to the
                       id. The return value is +/- one 

## Description

Function with one method that extracts the ith spin of the state corresponding
to the given id. Useful when constructing the hamiltonian matrix and comparing
diagonal/off-diagonal terms.
"""
function spin(id::Int, n::Int, i::Int)
  return digits(id-1, 2, n)[i]
end


"""
# inith(Array{Float64, 1}, Array{Float64, 1}, Bool, Int)

## Argument 

  * λ::Array{Float64, 1}        -   array containing the onsite magnetic fields.
  * Δ::Array{Float64, 1}        -   array of easy-axis anisotropies.
  * bc::Bool                    -   boundary conditions for the chain
    +  True = periodic
    +  False = open
  * n::Int                      -   system size

## Returns

  * H::Array{Float64, 2}    -   Hamiltonian for spin-1/2 Heisenberg chain  

## Description

This function construct the spin-1/2 Heisenberg Hamiltonian. The array passed to
this function is converted in place to the Hamiltonian for a spin-1/2 Heisenberg
chain. The onsite magnetic field is passed as an array to allow for disordered
fields.  
"""
function inith(λ::Array{Float64, 1}, Δ::Array{Float64, 1}, bc::Bool, n::Int)
  H = zeros(n, n) # initialize Hamiltonian
  for a=1:2^n
    bc ? m=n : m=n-1  # if the chain is open, only n-1 bonds. 
    for i=1:m
      j = (i%n)+1 # determine nearest neighbour
      if spin(a, n, i) == spin(a, n, j)
        H[a, a] = H[a, a]+Δ[i]*0.25+0.5*(λ[i]*spin(a, n, i)+λ[j]*spin(a, n, i))
      else
        H[a, a] = H[a, a]-Δ[i]*0.25
        b = flip(a, i, j)
        H[a, b] = H[a, b]+0.5 
      end
    end
  end
  return H
end

function inith(bc::Bool, n::Int)
  H = zeros(2^n, 2^n) # initialize Hamiltonian
  for a=1:2^n
    bc ? m=n : m=n-1  # if the chain is open, only n-1 bonds. 
    for i=1:m
      j = (i%n)+1 # determine nearest neighbour
      if spin(a, n, i) == spin(a, n, j)
        H[a, a] = H[a, a]+0.25
      else
        H[a, a] = H[a, a]-0.25
        b = flip(a, i, j)
        H[a, b] = H[a, b]+0.5 
      end
    end
  end
  return H
end

end # module EDUtils
