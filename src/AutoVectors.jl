module AutoVectors

export AutoVector, autovector, mini, maxi, clear!, copy, avdot, doprint, axpy!, convolve, 
		makeauto,makeautotake,applyshift,subav, avtriple, fast, avrange, arange,
		avlocation, avlocmin,avlocmax,avvec, shrink!, avnorm, avtripconv,reverse_ind,
		makeAutoVectorOfVecs,pointmult,fftav,ifftav

using LinearAlgebra
using DSP		# for conv = convolution
using AbstractFFTs,FFTW

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
    AutoVector(def=0.0,mini::Integer=1,maxi::Integer=0,miniloc::Integer=0)
    AutoVector(f::Function,mini::Integer=1,maxi::Integer=0,miniloc::Integer=0)
    AutoVector(v::Vector,mini::Integer=1,maxi::Integer=0,miniloc::Integer=0)

Frequently you just use

    v = AutoVector(), defaulting to Float64, size 0
or
    v = AutoVector(0), defaulting to Int64, size 0

An AutoVector expands when written to outside its range. Reading outside its range 
does not expand the range, and gives def.

Arguments:

def--default element, usually 0.0. Determines the type T of AutoVector{T}

mini and maxi give the index range of the created AutoVector (logical indices, not index in data vector)

miniloc is the location of mini-1 within the data vector, default 0 (data index)

You can initialize an AutoVector with the default, from a function, or by putting in a vector.
Most functions and constructors deal with logical indices,
Data indices refers to location within the data vector dat, which should not be something to
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

Create an AutoVector that holds Vectors as elements, an AutoVecotr of Vectors. The zero default Vector
is veczero.
"""
function makeAutoVectorOfVecs(veczero::Vector,mini::Integer,maxi::Integer)
    dat = [deepcopy(veczero)]
    AutoVector(veczero,mini,maxi,0,dat)
end

"""
    mini(v::AutoVector)
Minimum logical index.
"""
mini(v::AutoVector) = v.mini

"""
    maxi(v::AutoVector)
Maximum logical index.
"""
maxi(v::AutoVector) = v.maxi

import Base.length

"""
    length(v::AutoVector)
Logical length.
"""
length(v::AutoVector) = v.maxi-v.mini+1

"""
    avrange(v::AutoVector)
The logical range as mini:maxi.  A synonym is arange.
"""
avrange(v::AutoVector) = mini(v):maxi(v)

arange = avrange

"""
    olaprange(v::AutoVector,w::AutoVector)
olaprange(v,w) gives range of overlapping indices of AutoVectors v and w, given as a:b.
"""
olaprange(v::AutoVector,w::AutoVector) = max(mini(v),mini(w)):min(maxi(v),maxi(w))

"""
    avlocation(v::AutoVector,i)
The data location of logical index i.
"""
avlocation(v::AutoVector,i) = i-v.mini+v.miniloc+1

"""
    avlocmin(v::AutoVector)
The data location of mini(v).
"""
avlocmin(v::AutoVector) = v.miniloc+1

"""
    avlocmax(v::AutoVector)
The data location of maxi(v).
"""
avlocmax(v::AutoVector) = v.maxi-v.mini+v.miniloc+1

"""
    avvec(v::AutoVector)
avvec(v) creates a standard Vector with all values from mini to maxi; the new Vector does not know mini;
it only knows its length, given by maxi-mini+1.
"""
avvec(v::AutoVector) = v.dat[avlocmin(v):avlocmax(v)]

function getindex(v::AutoVector,i::Integer)
	if length(v.dat) == 0 || i < v.mini || i > v.maxi 
		return v.def
	end
	v.dat[avlocation(v,i)]
end

"""
    fast(v::AutoVector,i)
fast(v,i) is like accessing v[i], but without the check for being outside the logical range. 
If i is outside, fast  may or may not throw a standard range exception depending on whether i
lands outside the data vector's range. For use when i is known to be inside the logical range. 
Combined with @inbounds, there will be no range checking at all, for optimal speed.
"""
function fast(v::AutoVector,i)
	v.dat[i-v.mini+v.miniloc+1]
end

"""
    clear!(v::AutoVector)
Resets v to empty.
"""
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

import Base.*, Base./

"Multiply all elements of a AutoVector by a Float64"
function *(f::Real,v::AutoVector)
    newv = deepcopy(v)
    newv.dat .*= f
    newv
end
function *(v::AutoVector,f::Real) f*v end
function /(v::AutoVector,f::Real) v*(1/f) end

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

"""
    avdot(x::AutoVector,y::AutoVector)

Dot product (with no complex conjugating), sum_i x_i y_i over the intersection of ranges.
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

Triple dot product (with no complex conjugating), sum_i x_i y_i z_i over the intersection of ranges.
"""
function avtriple(x::AutoVector,y::AutoVector,z::AutoVector)
    res = 0.0
    for j = max(mini(x),mini(y),mini(z)):min(maxi(x),maxi(y),maxi(z))
@inbounds	res += fast(x,j) * fast(y,j) * fast(z,j)
    end
    res
end

"""
    doprint([file descriptor],v::AutoVector; spacing = 1)

Print all the elements to standard output or a file, with indices scaled by spacing.
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
Convolve with cutoff; no writing of elements if abs(value) < cut
"""
function convolve(x::AutoVector,y::AutoVector,cut=1.0e-14)		# use absolute cutoff
    res = AutoVector(x.def)
    if typeof(x.def) == Float64
	mix = mini(x)
	mx = maxi(x)
	mx-mix < 0 && return res
	xx = avvec(x)
	miy = mini(y)
	my = maxi(y)
	my-miy < 0 && return res
	yy = avvec(y)
	vres = myconv(xx,yy)
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
    println("Warning: convolving non-Float64 data")
	for j = mini(x):maxi(x)
	    for k = mini(y):maxi(y)
		res[j+k] += x[j] * y[k]
	    end
	end
    end
    res
end

"""
    makeauto(v::Vector{T};offset=nothing, firstindex=nothing, cutoff=0.0) where T
Make an AutoVector out of a Vector, producing a new data vector.
If offset is supplied, shifts data to left by offset.
If firstindex is supplied, it makes mini=firstindex
You can't specify both offset and firstindex.
With cutoff nonzero, elements are only put in if abs(el) > cutoff
To put in a part of a Vector at a particular range, say putting elements 2 to 4 at
positions 5 to 7, do this:
makeauto(v[2:4],firstindex=5)
"""
function makeauto(v::Vector{T};offset=nothing, firstindex=nothing, cutoff=0.0) where T
    res = AutoVector(zero(T))
    if !(offset == nothing || firstindex == nothing)
	error("You can't specify both offset and firstindex in makeauto")
    end
    offset == nothing && firstindex != nothing && (offset = 1-firstindex)
    offset == nothing && firstindex == nothing && (offset = 0)
    if cutoff == 0.0
	for i=1:length(v)
	    res[i-offset] = v[i]
	end
    else
	for i=1:length(v)
	    abs(v[i]) > cutoff && (res[i-offset] = v[i])
	end
    end
    res
end

"""
    avnorm(v::AutoVector)
Norm of AutoVector
"""
avnorm(v::AutoVector) = norm(v.dat[avlocmin(v):avlocmax(v)])

"""
    applyshift(x::AutoVector,offset::Integer)
Create a new AutoVector shifted to left by offset, sharing the same data array.
"""
function applyshift(x::AutoVector,offset::Integer)
    T = typeof(x.def)
    AutoVector{T}(x.def,x.mini+offset,x.maxi+offset,x.miniloc,x.dat)
end

"""
    subav(x::AutoVector,newmini::Integer,newmaxi::Integer)
Create a new AutoVector with the indicated range, sharing the same data array
"""
function subav(x::AutoVector,newmini::Integer,newmaxi::Integer)
    @assert newmini >= mini(x)
    @assert newmaxi <= maxi(x)
    minshift = newmini-mini(x)
    T = typeof(x.def)
    AutoVector{T}(x.def,newmini,newmaxi,x.miniloc+minshift,x.dat)
end

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

Adjust mini and maxi to omit tails where all elements are less than cut in absolute value.
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

"""
    reverse_ind(x::AutoVector)

A reflection about the origin; the new AutoVector goes from -maxi to -mini.
"""
function reverse_ind(x::AutoVector)
    newminiloc = length(x.dat) - (x.miniloc + x.maxi - x.mini + 1)
    T = typeof(x.def)
    AutoVector{T}(deepcopy(x.def),-x.maxi,-x.mini,newminiloc,vec(reverse(x.dat)))
end

"""
    fftav(f::AutoVector{Float64},delta)

For a real function f(x) defined by sampling on a uniform grid with spacing delta, 
returns a tuple (F, freqsp, len), where F is the Fourier Transform as an AutoVector(ComplexF64), 
freqsp is the frequency spacing of F, and len the length of the vector used in
the FFT (to be used in ifftav). Implicitly, the function is assumed to be zero
outside the range of the input. The normalizations are different from the usual discrete
FT conventions, incorporating the spacing, and corresponding to a continuous integral.
The FT is defined as

``F(k) = \\frac{1}{\\sqrt{2 \\pi}} \\int dx e^{-i k x} f(x)``

"""
function fftav(f::AutoVector{Float64},delta)
    maxind = max(abs(maxi(f)),abs(mini(f)))
    len = 2*maxind+1
    v = zeros(len)
    for i=0:maxi(f)
        v[i+1] = f[i]
    end
    for i=-1:-1:mini(f)
        v[len+i+1] = f[i]
    end
    v .*= delta/sqrt(2*pi)
    o = fftshift(fft(v))
    oa = AutoVector(0.0+0.0im)
    ishift = div(length(o),2)+1
    for i=1:length(o)
        abs(o[i]) > 1.0e-12 && (oa[i-ishift] = o[i])
    end
    freqspacing = 2*pi/len/delta
    (oa,freqspacing,len)
end

"""
ifftav(F::AutoVector{ComplexF64},freqspacing,maxind)

For a function F_i defined on a uniform frequency grid with spacing freqspacing, 
returns the Inverse Fourier Transform.

``f(x) = \\frac{1}{\\sqrt{2 \\pi}} \\int dk e^{i k x} F(k)``

"""
function ifftav(F::AutoVector{ComplexF64},freqspacing::Float64,len::Integer)
    v = zeros(ComplexF64,len)
    for i=0:maxi(x)
        v[i+1] = F[i]
    end
    for i=-1:-1:mini(F)
        v[len+i+1] = F[i]
    end
    o = fftshift(ifft(v))
    oa = AutoVector(0.0)
    ishift = div(length(o),2)+1
    for i=1:length(o)
        abs(o[i]) > 1.0e-12 && (oa[i-ishift] = real(o[i]))
    end
    oa.dat *= freqspacing / sqrt(2*pi) * len
    oa
end



const autovector = AutoVector


end
