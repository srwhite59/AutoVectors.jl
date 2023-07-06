module AutoVectors
using LinearAlgebra
using DSP		# for conv = convolution

import Base.deepcopy, Base.getindex, Base.setindex!

mutable struct AutoVector{T}  			#   <: AbstractVector{T}   doesn't work
    mini::Int64
    maxi::Int64
    miniloc::Int64		# the location/index of mini-1 in dat
    def::T
    dat::Vector{T}

    function AutoVector{T}(def::T,mini::Integer,maxi::Integer,miniloc::Integer,
			   				dat::Vector{T}) where T
	length(dat) < maxi-mini+1 && error("dat too small in AutoVector constructor")
	new(mini,maxi,miniloc,deepcopy(def),dat)
    end
end

"""
    AutoVector(def::T=0.0,mini::Integer=1,maxi::Integer=0,miniloc::Integer=0) where T
    AutoVector(f::Function,mini::Integer=1,maxi::Integer=0,miniloc::Integer=0)
    AutoVector(v::Vector,mini::Integer=1,maxi::Integer=0,miniloc::Integer=0)

Frequently you just use

    v = AutoVector(), defaulting to Float64, size 0
or
    v = AutoVector(0), defaulting to Int64, size 0

An AutoVector expands when written to outside its range. Reading outside its range 
does not expand the range, and gives def, normally 0.0.
Arguments:

def--default element, usually 0.0. 
mini and maxi give the index range of the created AutoVector (logical indices, not index in data vector)
miniloc is the location of mini-1 within the data vector, default 0 (physical index)

You can initialize an AutoVector with the default, from a function, or by putting in a vector.
Most functions and constructors deal with logical indices,
Physical indices refers to location within the data vector dat, which should not be something to
worry about normally.
"""
function AutoVector(def=0.0,mini::Integer=1,maxi::Integer=0,miniloc::Integer=0)
    T = typeof(def)
    AutoVector{T}(deepcopy(def),mini,maxi,miniloc,T[deepcopy(def) for j=1:maxi-mini+1])
end

function AutoVector(f::Function,mini::Integer=1,maxi::Integer=0,miniloc::Integer=0)
    dat = [f(i) for i=mini:maxi]
    T = typeof(dat[1])
    AutoVector{T}(zero(T),mini,maxi,miniloc,dat)
end

function AutoVector(v::Vector,mini::Integer=1,maxi::Integer=0,miniloc::Integer=0)
    maxi == 0 && (maxi = mini-1+length(v))
    T = typeof(v[1])
    AutoVector{T}(zero(T),mini,maxi,miniloc,v)
end


"""
    makeAutoVectorOfVecs(veczero::Vector,mini::Integer,maxi::Integer)

Create an AutoVector that holds vectors as elements, index from mini to maxi by v where the zero default vector
is veczero
"""
function makeAutoVectorOfVecs(veczero::Vector,mini::Integer,maxi::Integer)
    dat = [deepcopy(veczero)]
    AutoVector(veczero,mini,maxi,0,dat)
end

"minimum index"
mini(v::AutoVector) = v.mini

"maximum index"
maxi(v::AutoVector) = v.maxi

import Base.length

"length of an AutoVector"
length(v::AutoVector) = v.maxi-v.mini+1

"range of an AutoVector given as a:b"
arange(v::AutoVector) = mini(v):maxi(v)

"range of overlapping indices of two AutoVectors, given as a:b"
olaprange(v::AutoVector,w::AutoVector) = max(mini(v),mini(w)):min(maxi(v),maxi(w))

"physical location of logical index i)"
avlocation(v::AutoVector,i) = i-v.mini+v.miniloc+1

"physical location of mini(v)"
avlocmin(v::AutoVector) = v.miniloc+1

"physical location of maxi(v)"
avlocmax(v::AutoVector) = v.maxi-v.mini+v.miniloc+1

"convert to standard vector, all values from mini to maxi"
avvec(v::AutoVector) = v.dat[avlocmin(v):avlocmax(v)]

"get the value, but return def outside the range"
function getindex(v::AutoVector,i::Integer)
	if length(v.dat) == 0 || i < v.mini || i > v.maxi 
		return v.def
	end
	v.dat[avlocation(v,i)]
end

"get the value, but outside bounds throws range exception"
function fast(v::AutoVector,i)
	v.dat[i-v.mini+v.miniloc+1]
end

"Reset an AutoVector to empty"
function clear!(v::AutoVector{T}) where {T}
	v.mini = 1
	v.maxi = 0
	v.miniloc = 0
	v.dat = T[]
end

"Assign the value, but outside bounds forces resizing and adjustment of mini or maxi"
function setindex!(v::AutoVector{T},x,i::Integer) where {T}
	if(length(v.dat) == 0)
		v.dat = T[deepcopy(v.def) for j=1:19]
		v.miniloc = 10
		v.mini = v.maxi = i
	elseif(i < v.mini || i > v.maxi)
		newmini = min(i,v.mini)
		newmaxi = max(i,v.maxi)
		j = i - v.mini + v.miniloc
		if j >= 0 && j < length(v.dat)
			v.miniloc += newmini - v.mini;
		else
			newlen = (newmaxi-newmini+20)*2
			newminiloc=div(newlen-(newmaxi-newmini+1),2)
			newdat = T[deepcopy(v.def) for j=1:newlen]
			oldbegin = v.miniloc
			oldend = v.miniloc+v.maxi-v.mini
			newbegin = newminiloc - newmini + v.mini
			for k = oldbegin:oldend
				newdat[newbegin+k-oldbegin+1] = v.dat[k+1]
			end
			v.dat = newdat
			v.miniloc = newminiloc
		end
		v.mini = newmini
		v.maxi = newmaxi
	end
	v.dat[i-v.mini+v.miniloc+1] = x
end

import Base.copy
function copy(x::AutoVector)
    deepcopy(x)
end

import Base.+,Base.-

"Add AutoVectors"
function +(v::AutoVector,w::AutoVector)
    AutoVector(i->v[i]+w[i],min(v.mini,w.mini),max(v.maxi,w.maxi))
end

"Subtract AutoVectors"
function -(v::AutoVector,w::AutoVector)
    AutoVector(i->v[i]-w[i],min(v.mini,w.mini),max(v.maxi,w.maxi))
end

"Add a Float64 to all elements of a AutoVector"
function +(v::AutoVector,w::Float64)
    AutoVector(i->v[i] .+ w,v.mini,v.maxi)
end
"Subtract a Float64 from all elements of a AutoVector{Float64}"
function -(v::AutoVector,w::Float64)
    AutoVector(i->v[i] .- w,v.mini,v.maxi)
end

import Base.*

"Multiply all elements of a AutoVector by a Float64"
function *(f::Real,v::AutoVector)
    newv = deepcopy(v)
    newv.dat *= f
    newv
end
function *(v::AutoVector,f::Real) f*v end

"""
    *(coef::Vector{Float64},v::Vector{AutoVector{Float64}})

Add a vector of AutoVectors together scaled by coefficients
"""
function *(coef::Vector{Float64},v::Vector{AutoVector{Float64}})
    n = length(v)
    mi = minimum([mini(v[i]) for i=1:n])
    ma = maximum([maxi(v[i]) for i=1:n])
    len = ma-mi+1
    res = zeros(len)
    for i=1:n
	ra = arange(v[i]) .+ (1-mi)
	res[ra] += coef[i] * avvec(v[i])
    end
    AutoVector(res,mi,ma)
end

function pointmult(x::AutoVector,y::AutoVector)
    a = max(mini(x),mini(y))
    b = min(maxi(x),maxi(y))
    @inbounds xy = Float64[fast(x,i)*fast(y,i) for i=a:b]
    AutoVector(xy,a,b)
end
#=
import Base.iterate
function iterate(v::AutoVector)
    length(v) == 0 && return nothing
    (v[mini(v)],mini(v))
end
function iterate(v::AutoVector,i)
    length(v) == 0 && return nothing
    i >= maxi(v) && return nothing
    (v[i+1],i+1)
end
#import Base: .*
#function .*(v::AutoVector,w::AutoVector)
#    pointmult(x,y)
#end

import Base.broadcast
#Base.broadcast(::typeof(*), ...)"
#function broadcast(*,x::AutoVector,y::AutoVector)
function broadcast(::typeof(*),x::AutoVector,y::AutoVector)
    pointmult(x,y)
end
=#

"""
    avdot(x::AutoVector,y::AutoVector)

Dot product (no conjugating)
"""
function avdot(x::AutoVector,y::AutoVector)
    a = max(mini(x),mini(y))
    b = min(maxi(x),maxi(y))
    a > b && return 0.0

    xoff = -x.mini + x.miniloc + 1
    yoff = -y.mini + y.miniloc + 1
    @views rf = dot(x.dat[a+xoff:b+xoff],y.dat[a+yoff:b+yoff])
    rf
end

"""
    avtriple(x::AutoVector,y::AutoVector,z::AutoVector)

Triple dot product (no conjugating)
"""
function avtriple(x::AutoVector,y::AutoVector,z::AutoVector)
    res = 0.0
    for j = max(mini(x),mini(y),mini(z)):min(maxi(x),maxi(y),maxi(z))
@inbounds	res += fast(x,j) * fast(y,j) * fast(z,j)
    end
    res
end

"""
    doprint([file io thing],v::AutoVector; spacing = 1)

Print to standard output or a file all the elements
"""
function doprint(v::AutoVector; spacing = 1)
    for i=mini(v):maxi(v)
        println(i*spacing," ",v[i])
    end
end
function doprint(FI,v::AutoVector; spacing = 1)
    for i=mini(v):maxi(v)
        println(FI,i*spacing," ",v[i])
    end
end

import LinearAlgebra.axpy!
"""
    axpy!(y::AutoVector,a::Float64,x::AutoVector)	# y += a * x
"""
function axpy!(y::AutoVector,a::Float64,x::AutoVector)	# y += a * x
    for k = mini(x):maxi(x)
	y[k] += a * x[k]
    end
end

"""
    axpy!(y::AutoVector,a::Float64,x::AutoVector, cutoff::Float64)	# y += a * x with cutoff for writing
"""
function axpy!(y::AutoVector,a::Float64,x::AutoVector, cutoff::Float64)	# y += a * x
    for k = mini(x):maxi(x)
	r = a * x[k]
	abs(r) > cutoff && (y[k] += a * x[k])
    end
end

function dotrip(ud,gd,vd,ua,ub,ga,va,gb,vb,uoff,goff,voff)
    res = 0.0
		    #j+k >= va  so k >= va - j
		    #j+k <= vb  so k <= vb - j
    for j = ua:ub, k = max(ga,va-j):min(gb,vb-j)
	@inbounds res += ud[j-uoff] * gd[k-goff] * vd[j+k-voff]
    end
    res
end

"""
    avtripconv(u::AutoVector,g::AutoVector,v::AutoVector)
Same as avdot(convolve(u,g),v)
"""
function avtripconv(u::AutoVector,g::AutoVector,v::AutoVector)
    ua,ub = mini(u),maxi(u)
    ga,gb = mini(g),maxi(g)
    va,vb = mini(v),maxi(v)
    ub+gb < va && return 0.0
    ua+ga > vb && return 0.0
    dotrip(u.dat,g.dat,v.dat,ua,ub,ga,va,gb,vb,u.mini-u.miniloc-1,g.mini-g.miniloc-1,v.mini-v.miniloc-1)
end

function tripconv(u::Vector{Float64},g::Vector{Float64},v::Vector{Float64})
    k,m,n = length(u),length(g),length(v)
    res = 0.0
    for i=1:min(k,n), j=1:min(m,n-i+1)
	@inbounds res += u[i] * g[j] * v[i+j-1]
    end
    res
end

function convolvecheck(x::Vector{Float64},y::Vector{Float64})
    mx = length(x)
    my = length(y)
    res = zeros(Float64,mx+my-1)
    for j = 1:mx
@inbounds	for k = 1:my
	    res[j+k-1] += x[j] * y[k]
	end
    end
    res
end

function conv3(u::Vector{Float64},v::Vector{Float64})   # want v shorter
    m,n = length(u),length(v)
    s = m+n-1
    res = zeros(s)
    for j=1:n
	for i=1:m
            @inbounds res[i+j-1] += u[i] * v[j]
        end
    end
    res
end

function myconv(u::Vector{Float64},v::Vector{Float64})
    m,n = length(u),length(v)
    if m*n > 7000		# Crossover seems to depend on total effort (m*n) more than n
        return conv(u,v)
    end
    if m < n
        return conv3(v,u)
    else
	return conv3(u,v)
    end
end

"""
    convolve(x::AutoVector,y::AutoVector,cut=1.0e-14)		# use absolute cutoff
Convolve with cutoff
"""
function convolve(x::AutoVector,y::AutoVector,cut=1.0e-14)		# use absolute cutoff
    res = AutoVector(x.def)
#println(maxi(x)-mini(x),"  ",maxi(y)-mini(y))
    if typeof(x.def) == Float64
	mix = mini(x)
	mx = maxi(x)
	mx-mix < 0 && return res
	#xx = Float64[fast(x,i+mix-1) for i=1:mx-mix+1]
	xx = avvec(x)
	miy = mini(y)
	my = maxi(y)
	my-miy < 0 && return res
	#yy = Float64[fast(y,i+miy-1) for i=1:my-miy+1]
	yy = avvec(y)
	#@time vres2 = convolvecheck(xx,yy)
	vres = myconv(xx,yy)
	#@show norm(vres-vres2)
	mir = mix+miy
	minj = 1
	for j = 1:length(vres)
	    abs(vres[j]) > cut && (minj = j; break)
	end
	maxj = length(vres)
	for j = length(vres):-1:1
	    abs(vres[j]) > cut && (maxj = j; break)
	end
	res.dat = vres[minj:maxj]
	res.mini=minj+mir-1
	res.maxi=maxj+mir-1
	res.miniloc=0
    else
    println("convolving non-Float64 data")
	for j = mini(x):maxi(x)
	    for k = mini(y):maxi(y)
		res[j+k] += x[j] * y[k]
	    end
	end
    end
    res
end

#function makeautotake(v::Vector{Float64},offset::Integer)		# take v as the data; very fast
#    AutoVector(0.0,1-offset,length(v)-offset,0,v)
#end

"""
    makeauto(v::Vector{Float64},offset::Integer)
Convmake AutoVector out of vector by shifting to left by offset
"""
function makeauto(v::Vector{Float64},offset::Integer)
    res = AutoVector(0.0)
    for i=1:length(v)
	res[i-offset] = v[i]
    end
    res
end
function makeauto(v::Vector{Float64},offset::Integer,cutoff::Float64)
    res = AutoVector(0.0)
    for i=1:length(v)
	vi = v[i]
	abs(vi) > cutoff && ( res[i-offset] = v[i] )
    end
    res
end

"""
    avnorm(v::AutoVector)
norm of AutoVector
"""
avnorm(v::AutoVector) = norm(v.dat[avlocmin(v):avlocmax(v)])

"""
    applyshift(x::AutoVector,offset::Integer)
Shift to left by offset, no new data array
"""
function applyshift(x::AutoVector,offset::Integer)
    AutoVector(x.def,x.mini+offset,x.maxi+offset,x.miniloc,x.dat)
end

#=
function symmetrize!(m::Array{Float64,2})
    n = size(m,1)
    for i=1:n
	for j=i:n
	    m[i,j] = m[j,i] = 0.5 * (m[i,j]+m[j,i])
	end
    end
end
=#

function domin(xd,st,ma,cut)
    @inbounds for i=st:ma
	if(abs(xd[i]) > cut)
	    return i
	end
    end
    return length(xd)+1
end
function domax(xd,st,ma,cut)
    @inbounds for i=ma:-1:st
	if(abs(xd[i]) > cut)
	    return i
	end
    end
    return 0
end

"""
    shrink!(x::AutoVector,cut)

new mini and maxi to zero out tails less than cut
"""
function shrink!(x::AutoVector,cut)
    ami,ama = avlocmin(x),avlocmax(x)
    mi = domin(x.dat,ami,ama,cut)
    newmin = mini(x)+mi-ami
    if newmin > maxi(x) 
	clear!(x)
	return
    end
    ma = domax(x.dat,ami,ama,cut)
    x.maxi = maxi(x) + ma - ama 
    x.miniloc += newmin - x.mini
    x.mini = newmin
end
function shrinkold!(x::AutoVector,cut)
    newmin = maxi(x)+1
    for i=arange(x)
	if(abs(fast(x,i)) > cut)
	    newmin = i
	    break
	end
    end
    if newmin == maxi(x)+1 
	clear!(x)
	return
    end
    newmax = newmin-1
    for i=maxi(x):-1:newmin
	if(abs(fast(x,i)) > cut)
	    newmax = i
	    break
	end
    end
    x.maxi = newmax
    x.miniloc += newmin - x.mini
    x.mini = newmin
end

#=
function eigsym(Marg)
    M = copy(Marg)
    symmetrize!(M)
    F = eigen(M)
    F.values,F.vectors
end
=#

"""
    reverse_ind(x::AutoVector)

new AutoVector goes from -maxi to -mini;  reflection really
"""
function reverse_ind(x::AutoVector)
    newminiloc = length(x.dat) - (x.miniloc + x.maxi - x.mini + 1)
    AutoVector(deepcopy(x.def),-x.maxi,-x.mini,newminiloc,vec(reverse(x.dat)))
end

const autovector = AutoVector

export AutoVector, autovector, mini, maxi, clear!, copy, avdot, doprint, axpy!, convolve, 
		makeauto,makeautotake,applyshift,avtriple, fast, arange, #  symmetrize!, 
		avlocation, avlocmin,avlocmax,avvec, shrink!, avnorm, avtripconv,reverse_ind,
		makeAutoVectorOfVecs,pointmult

end
