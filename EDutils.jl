module EDutils

export flip
export getid
export getstate
export inith
export initσz
export initσx

"""
# flip(Int, Int, Int)

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
  return digits(id-1, base=2, pad=n)[i]
end


"""
# inith - 5 method function

## Methods

  * inith(bc::Bool, n::Int)
  * inith(λ::Float64, bc::Bool, n::Int)
  * inith(λ::Float64, Δ::Float64, bc::Bool, n::Int)
  * inith(λ::Array{Float64, 1}, bc::Bool, n::Int)
  * inith(λ::Array{Float64, 1}, Δ::Array{Float64, 1}, bc::Bool, n::Int)

# Arguments 

  * λ::Array{Float64, 1}    -   Uniform magnetic field 
  * λ::Array{Float64, 1}    -   array containing the onsite magnetic fields.
  * Δ::Array{Float64, 1}    -   array of easy-axis anisotropies.
  * bc::Bool                -   boundary conditions for the chain
    +  True = periodic
    +  False = open
  * n::Int                  -   system size

# Returns

  * H::Array{Float64, 2}    -   Hamiltonian for spin-1/2 Heisenberg chain  

## Description

This function construct the spin-1/2 Heisenberg Hamiltonian. The array passed to
this function is converted in place to the Hamiltonian for a spin-1/2 Heisenberg
chain. The onsite magnetic field is passed as an array to allow for disordered
fields.  
"""
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

function inith(λ::Float64, bc::Bool, n::Int)
  H = zeros(2^n, 2^n) # initialize Hamiltonian
  for a=1:2^n
    bc ? m=n : m=n-1  # if the chain is open, only n-1 bonds. 
    for i=1:m
      j = (i%n)+1 # determine nearest neighbour
      if spin(a, n, i) == spin(a, n, j)
        H[a, a] = H[a, a]+0.25-0.5*λ*(spin(a, n, i)+spin(a, n, j))
      else
        H[a, a] = H[a, a]-0.25-0.5*λ*(spin(a, n, i)+spin(a, n, j))
        b = flip(a, i, j)
        H[a, b] = H[a, b]+0.5 
      end
    end
  end
  return H 
end

function inith(λ::Float64, Δ::Float64, bc::Bool, n::Int)
  H = zeros(2^n, 2^n) # initialize Hamiltonian
  for a=1:2^n
    bc ? m=n : m=n-1  # if the chain is open, only n-1 bonds. 
    for i=1:m
      j = (i%n)+1 # determine nearest neighbour
      if spin(a, n, i) == spin(a, n, j)
        H[a, a] = H[a, a]+Δ*0.25-0.5*λ*(spin(a, n, i)+spin(a, n, j))
      else
        H[a, a] = H[a, a]-Δ*0.25-0.5*λ*(spin(a, n, i)+spin(a, n, j))
        b = flip(a, i, j)
        H[a, b] = H[a, b]+0.5 
      end
    end
  end
  return H 
end

function inith(λ::Array{Float64, 1}, bc::Bool, n::Int)
  H = zeros(2^n, 2^n) # initialize Hamiltonian
  for a=1:2^n
    bc ? m=n : m=n-1  # if the chain is open, only n-1 bonds. 
    for i=1:m
      j = (i%n)+1 # determine nearest neighbour
      if spin(a, n, i) == spin(a, n, j)
        H[a, a] = H[a, a]+0.25-0.5*(λ[i]*spin(a, n, i)+λ[j]*spin(a, n, j))
      else
        H[a, a] = H[a, a]-0.25-0.5*(λ[i]*spin(a, n, i)+λ[j]*spin(a, n, j))
        b = flip(a, i, j)
        H[a, b] = H[a, b]+0.5 
      end
    end
  end
  return H
end

function inith(λ::Array{Float64, 1}, Δ::Array{Float64, 1}, bc::Bool, n::Int)
  H = zeros(2^n, 2^n) # initialize Hamiltonian
  for a=1:2^n
    bc ? m=n : m=n-1  # if the chain is open, only n-1 bonds. 
    for i=1:m
      j = (i%n)+1 # determine nearest neighbour
      if spin(a, n, i) == spin(a, n, j)
        H[a, a] = H[a, a]+Δ[i]*0.25-0.5*(λ[i]*spin(a, n, i)+λ[j]*spin(a, n, j))
      else
        H[a, a] = H[a, a]-Δ[i]*0.25-0.5*(λ[i]*spin(a, n, i)+λ[j]*spin(a, n, j))
        b = flip(a, i, j)
        H[a, b] = H[a, b]+0.5 
      end
    end
  end
  return H
end

"""
initσz - 1 method function

## Methods

  * initσz(n::Int, i::Int)

# Arguments

  * n::Int    -   system size.
  * i::Int    -   site on which sz operator acts.

# Returns

  *sz::Array{Float64, 2}    -   The sz array spin-half particle at position i of
                                length n heisenberg chain.

## Description

This function constructs the σ_z operator acting on the ith site of a spin-1/2
chain.

"""
function initσz(n::Int, i::Int)
  
  id = [[1, 0] [0, 1]]
  σz = [[1, 0] [0, -1]]
  tot = Array{Int64, 2}

  if i==1
    tot = σz
    for i=2:n
      tot = kron(tot, id)
    end
  else
    tot = id
    for j=2:(i-1)
      tot = kron(tot, id)
    end
    tot = kron(tot, σz)
    for j=(i+1):n
      tot = kron(tot, id)
    end
  end

  return tot

end

"""
initσx - 1 method function

## Methods

  * initσx(n::Int, i::Int)

# Arguments

  * n::Int    -   system size.
  * i::Int    -   site on which sz operator acts.

# Returns

  *sx::Array{Float64, 2}    -   The sz array spin-half particle at position i of
                                length n heisenberg chain.

## Description

This function constructs the σ_x operator acting on the ith site of a spin-1/2
chain.

"""
function initσx(n::Int, i::Int)
  
  id = [[1, 0] [0, 1]]
  σx = [[0, 1] [1, 0]]
  tot = Array{Int64, 2}

  if i==1
    tot = σx
    for i=2:n
      tot = kron(tot, id)
    end
  else
    tot = id
    for j=2:(i-1)
      tot = kron(tot, id)
    end
    tot = kron(tot, σx)
    for j=(i+1):n
      tot = kron(tot, id)
    end
  end

  return tot

end

"""
initσy - 1 method function

## Methods

  * initσy(n::Int, i::Int)

# Arguments

  * n::Int    -   system size.
  * i::Int    -   site on which sz operator acts.

# Returns

  *sx::Array{Float64, 2}    -   The sz array spin-half particle at position i of
                                length n heisenberg chain.

## Description

This function constructs the σ_y operator acting on the ith site of a spin-1/2
chain.

"""
function initσy(n::Int, i::Int)
  
  id = [[1, 0] [0, 1]]
  σx = [[0, -1im] [im, 0]]
  tot = Array{Int64, 2}

  if i==1
    tot = σx
    for i=2:n
      tot = kron(tot, id)
    end
  else
    tot = id
    for j=2:(i-1)
      tot = kron(tot, id)
    end
    tot = kron(tot, σx)
    for j=(i+1):n
      tot = kron(tot, id)
    end
  end

  return tot

end

"""
initMagBlock - 5 method function

# Methods

  * initMagBlock(bc::Bool, n::Int, mz::Int)
  * initMagBlock(λ::Float64, bc::Bool, n::Int, mz::Int)
  * initMagBlock(λ::Float64, Δ::Float64, bc::Bool, n::Int, mz::Int)
  * initMagBlock(λ::Array{Float64, 1}, bc::Bool, n::Int, mz::Int)
  * initMagBlock(λ::Array{Float64, 1}, Δ::Array{Float64, 1}, bc::Bool, n::Int,
                 mz::Int)

## Arguments 

  * λ::Float64              -   Uniform magnetic field 
  * Δ::Float64              -   Uniform uniaxial anisotropy
  * λ::Array{Float64, 1}    -   array containing the onsite magnetic fields.
  * Δ::Array{Float64, 1}    -   array of easy-axis anisotropies.
  * bc::Bool                -   boundary conditions for the chain
    +  True = periodic
    +  False = open
  * n::Int                  -   system size
  * mz::Int                  -   magnetization subspace

## Returns

  * H::Array{Float64, 2}    -   Hamiltonian for spin-1/2 Heisenberg chain  

# Description

This function initializes a block of the Hamiltonian that leverage the
conservation of magnetization of the Hamiltonian. 
"""

function initMagBlock(bc::Bool, n::Int, mz::Int)
  
  # select all states with given magnetization mz
  mz_states = []
  dm::Int = 0 
  for a=1:2^n
    tot::Int = 0
    for i=1:n
      if 1==spin(a, n, i)
        tot += 1
      else
        tot -= 1
      end
    end
    if tot==mz
      dm += 1
      push!(mz_states, a)  
    end
  end
 
  H = zeros(dm, dm) # initialize Hamiltonian
  for a=1:dm
    bc ? m=n : m=n-1  # if the chain is open, only n-1 bonds. 
    for i=1:m
      j = (i%n)+1 # determine nearest neighbour
      if spin(mz_states[a], n, i) == spin(mz_states[a], n, j)
        H[a, a] = H[a, a]+0.25
      else
        H[a, a] = H[a, a]-0.25
        b = flip(mz_states[a], i, j)
        b = findall(x->x==b, mz_states)[1]
        H[a, b] = H[a, b]+0.5 
      end
    end
  end
  return H, mz_states
end

function initMagBlock(λ::Float64, bc::Bool, n::Int, mz::Int)
  # select all states with given magnetization mz
  mz_states = []
  dm::Int = 0 
  for a=1:2^n
    tot::Int = 0
    for i=1:n
      if 1==spin(a, n, i)
        tot += 1
      else
        tot -= 1
      end
    end
    if tot==mz
      dm += 1
      push!(mz_states, a)  
    end
  end
 
  H = zeros(dm, dm) # initialize Hamiltonian
  for a=1:dm
    bc ? m=n : m=n-1  # if the chain is open, only n-1 bonds. 
    for i=1:m
      j = (i%n)+1 # determine nearest neighbour
      if spin(mz_states[a], n, i) == spin(mz_states[a], n, j)
        H[a, a] = H[a, a]+0.25-0.5*λ*(spin(mz_states[a], n, i)+spin(mz_states[a], n, j))
      else
        H[a, a] = H[a, a]-0.25-0.5*λ*(spin(mz_states[a], n, i)+spin(mz_states[a], n, j))
        b = flip(mz_states[a], i, j)
        b = findall(x->x==b, mz_states)[1]
        H[a, b] = H[a, b]+0.5 
      end
    end
  end
  return H, mz_states
end

function initMagBlock(λ::Float64, Δ::Float64, bc::Bool, n::Int, mz::Int)
  # select all states with given magnetization mz
  mz_states = []
  dm::Int = 0 
  for a=1:2^n
    tot::Int = 0
    for i=1:n
      if 1==spin(a, n, i)
        tot += 1
      else
        tot -= 1
      end
    end
    if tot==mz
      dm += 1
      push!(mz_states, a)  
    end
  end

  H = zeros(dm, dm) # initialize Hamiltonian
  for a=1:dm
    bc ? m=n : m=n-1  # if the chain is open, only n-1 bonds. 
    for i=1:m
      j = (i%n)+1 # determine nearest neighbour
      if spin(mz_states[a], n, i) == spin(mz_states[a], n, j)
        H[a, a] = H[a, a]+Δ*0.25-0.5*λ*(spin(mz_states[a], n, i)+spin(mz_states[a], n, j))
      else
        H[a, a] = H[a, a]-Δ*0.25-0.5*λ*(spin(mz_states[a], n, i)+spin(mz_states[a], n, j))
        b = flip(mz_states[a], i, j)
        b = findall(x->x==b, mz_states)[1]
        H[a, b] = H[a, b]+0.5 
      end
    end
  end
  return H, mz_states
end

function initMagBlock(λ::Array{Float64, 1}, bc::Bool, n::Int, mz::Int)
  # select all states with given magnetization mz
  mz_states = []
  dm::Int = 0 
  for a=1:2^n
    tot::Int = 0
    for i=1:n
      if 1==spin(a, n, i)
        tot += 1
      else
        tot -= 1
      end
    end
    if tot==mz
      dm += 1
      push!(mz_states, a)  
    end
  end

  H = zeros(dm, dm) # initialize Hamiltonian
  for a=1:dm
    bc ? m=n : m=n-1  # if the chain is open, only n-1 bonds. 
    for i=1:m
      j = (i%n)+1 # determine nearest neighbour
      if spin(a, n, i) == spin(a, n, j)
        H[a, a] = H[a, a]+0.25-0.5*(λ[i]*spin(mz_states[a], n, i)
                                   +λ[j]*spin(mz_states[a], n, j))
      else
        H[a, a] = H[a, a]-0.25-0.5*(λ[i]*spin(mz_states[a], n, i)
                                   +λ[j]*spin(mz_states[a], n, j))
        b = flip(mz_states[a], i, j)
        b = findall(x->x==b, mz_states)[1]
        H[a, b] = H[a, b]+0.5 
      end
    end
  end
  return H, mz_states
end

function initMagBlock(λ::Array{Float64, 1}, Δ::Array{Float64, 1}, bc::Bool,
                      n::Int, mz::Int)
  H = zeros(2^n, 2^n) # initialize Hamiltonian
  for a=1:2^n
    bc ? m=n : m=n-1  # if the chain is open, only n-1 bonds. 
    for i=1:m
      j = (i%n)+1 # determine nearest neighbour
      if spin(a, n, i) == spin(a, n, j)
        H[a, a] = H[a, a]+Δ[i]*0.25-0.5*(λ[i]*spin(mz_states[a], n, i)
                                        +λ[j]*spin(mz_states[a], n, j))
      else
        H[a, a] = H[a, a]-Δ[i]*0.25-0.5*(λ[i]*spin(mz_states[a], n, i)
                                        +λ[j]*spin(mz_states[a], n, j))
        b = flip(mz_states[a], i, j)
        b = findall(x->x==b, mz_states)[1]
        H[a, b] = H[a, b]+0.5 
      end
    end
  end
  return H, mz_states
end


end # module EDutils
