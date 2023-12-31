# AutoVectors.jl Documentation

```@contents
```

## Constructors

```@docs
AutoVector(0.0)
makeauto(v::Vector;offset=nothing, firstindex=nothing, cutoff=0.0)
makeAutoVectorOfVecs(veczero::Vector,mini::Integer,maxi::Integer)
```

## Functions

```@docs
mini(v::AutoVector)
maxi(v::AutoVector)
length(v::AutoVector)
arange(v::AutoVector)
avvec(v::AutoVector)
subav(v::AutoVector,i::Integer,j::Integer)
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
applyshift(x::AutoVector,offset::Integer)
shrink!(x::AutoVector,cut)
reverse_ind(x::AutoVector)
fftav(x::AutoVector{Float64},delta)
ifftav(x::AutoVector{ComplexF64},freqspacing::Float64,len::Integer)

```

## Index

```@index
```
