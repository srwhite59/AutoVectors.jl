# AutoVectors.jl Documentation

```@contents
```

## Constructors

```@docs
AutoVector(0.0)
makeAutoVectorOfVecs(veczero::Vector,mini::Integer,maxi::Integer)
```

## Functions

```@docs
mini(v::AutoVector)
maxi(v::AutoVector)
length(v::AutoVector)
arange(v::AutoVector)
avvec(v::AutoVector)
avdot(x::AutoVector,y::AutoVector)
convolve(x::AutoVector,y::AutoVector,cut=1.0e-14) 
doprint(v::AutoVector)
fast(v::AutoVector,i)
clear!(v::AutoVector)
avnorm(v::AutoVector)
avlocation(v::AutoVector,i)
avlocmin(v::AutoVector)
avlocmax(v::AutoVector)
avtriple(x::AutoVector,y::AutoVector,z::AutoVector)
axpy!(y::AutoVector,a::Float64,x::AutoVector) 
axpy!(y::AutoVector,a::Float64,x::AutoVector, cutoff::Float64) 
avtripconv(u::AutoVector,g::AutoVector,v::AutoVector)
makeauto(v::Vector{Float64},offset::Integer)
applyshift(x::AutoVector,offset::Integer)
shrink!(x::AutoVector,cut)
reverse_ind(x::AutoVector)

```

## Index

```@index
```
