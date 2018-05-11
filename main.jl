include("EDUtils.jl")

import EDUtils

n = 4

for a=1:2^n
  for i=1:n
    j=(i%n) + 1
    b = EDUtils.flip(a, i, j)
    println("Initial State:\t", a)
    println("Flip Sites:\t", i, " ", j)
    println("Final State:\t", b)
    println("--------------------")
  end
end


