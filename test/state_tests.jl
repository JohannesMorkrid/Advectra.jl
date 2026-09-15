# Constructors
using CUDA
using Advectra

# For this test will use (n, p, T) and the following domain
domain = Domain(256; MemoryType=CuArray)

# Fieldnames + Domain
S1 = State(:n, :p, :T, domain)
typeof(S1.data) <: CuArray #Should be CuArray
propertynames(S1) # Should include :n, :p, :T
size(S1) # Should be (256,256,3)
S1.domain === domain # Should be true

# Fieldnames (Tuple) + Domain
S2 = State((:n, :p, :T), domain)
typeof(S2.data) <: CuArray #Should be CuArray
propertynames(S2) # Should include :n, :p, :T
size(S2) == (256, 256, 3) # Should be (256,256,3) 
S2.domain === domain # Should be same domain by reference

# NamedTuple of args mixing Array, initial_condition and Field
S3 = State(domain; n=rand(256, 256), p=(gaussian, (; A=2)), T=Field(rand(256, 256), domain))
typeof(S3.data) <: CuArray
length(intersect(propertynames(S3), (:n, :p, :T))) == 3 # Should include :n, :p, :T
size(S3) == (256, 256, 3) # Should be (256,256,3)
S3.domain === domain # Should be same domain by reference

S4 = State((n=rand(256, 256), p=(gaussian, (; A=2)), T=Field(rand(256, 256), domain)),
    domain)
typeof(S4.data) <: CuArray
length(intersect(propertynames(S4), (:n, :p, :T))) == 3 # Should include :n, :p, :T
size(S4) == (256, 256, 3) # Should be (256,256,3)
S4.domain === domain # Should be same domain by reference

# Broadcasting test
S1 .= S2 .+ 0.5S3 .+ exp.(S4)
S1.data == S2.data .+ 0.5S3.data .+ exp.(S4).data

S3 .^ -1 # Check that works!

#zero, copy, similar, copyto!

sum(S3)
minimum(S3)
maximum(S3)
all(S3 .!= 4)
any(S3 .== 2.0)
using LinearAlgebra
norm(S3)
dot(S3, S4) == dot(S3, S4.data)
dot(S3, S4) == dot(S3.data, S4)

# map, mapreduce

S4 .= S3.data

S1 .= 0.0
fill!(S2, 1.0)

rmul!(S4, 2.0)

### Field
F1, F2, F3 = unpack_state(S4)

# Broadcasting test
F1 .= F2 .+ 0.5exp.(F3)
F1.data == F2.data .+ 0.5exp.(F3.data)

F3 .^ -1 # Check that works!

#zero, copy, similar, copyto!

sum(F3)
minimum(F3)
maximum(F3)
all(F3 .!= 4)
any(F3 .== 2.0)
using LinearAlgebra
norm(F3)
dot(F2, F3) == dot(F2, F3.data)
dot(F2, F3) == dot(F2.data, F3)

# map, mapreduce

F2 .= F3.data

F1 .= 0.0
fill!(F2, 1.0)

rmul!(F2, 2.0)
