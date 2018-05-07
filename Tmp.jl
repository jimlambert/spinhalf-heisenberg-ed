module Tmp

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

  * (none)

## Description

Function with one method. Accepts three integers, s, i, j and flips the ith and
jth bits of the binary representation of the integer s. The integer s is
modified in place.  
"""
function flip!(s::Int, i::Int, j::Int)
  f::Int64 = 2^i + 2^j
  println(f)
  return xor(s, f)
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
  return sum
end

function spin(id::Int64, n::Int64, i::Int64)
  return digits(id, 2, n)[i]
end

function initH!(H::Array{Int64, 2})

end

end # module Tmp
