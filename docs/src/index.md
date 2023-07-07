# AutoVectors.jl Documentation

```@contents
```

## Constructors

```@docs
AutoVector(0.0)
```

## Functions

```@docs
makeAutoVectorOfVecs(veczero::Vector,mini::Integer,maxi::Integer)
mini(v::AutoVector)
maxi(v::AutoVector)
length(v::AutoVector)
arange(v::AutoVector)
olaprange(v::AutoVector,w::AutoVector)
avlocation(v::AutoVector,i)
avlocmin(v::AutoVector)
avlocmax(v::AutoVector)
avvec(v::AutoVector)
fast(v::AutoVector,i)
clear!(v::AutoVector)
avdot(x::AutoVector,y::AutoVector)
avtriple(x::AutoVector,y::AutoVector,z::AutoVector)
doprint([file io thing],v::AutoVector; spacing = 1)
axpy!(y::AutoVector,a::Float64,x::AutoVector) 
axpy!(y::AutoVector,a::Float64,x::AutoVector, cutoff::Float64) 
avtripconv(u::AutoVector,g::AutoVector,v::AutoVector)
convolve(x::AutoVector,y::AutoVector,cut=1.0e-14) 
makeauto(v::Vector{Float64},offset::Integer)
avnorm(v::AutoVector)
applyshift(x::AutoVector,offset::Integer)
shrink!(x::AutoVector,cut)
reverse_ind(x::AutoVector)

```

## Index

```@index
```
