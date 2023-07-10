# AutoVectors.jl
Julia vectors whose first and last indices are any integers

Most languages make a choice: arrays start at 0 or 1.  For some applications, one would
like arrays that start and end anywhere, e.g. the array runs from -3 to 10. AutoVectors
have this property. They also auto-resize based on writing outside their range.  The range is always
contigous, and implemented through ordinary vectors, so AutoVectors are fast, say, compared to
a Dict{Int64,Float64} which might implement similar featues.  Reading outside the range gives 0.0 
(or whatever the default element is).  Often these properties mean that you
don't have to worry about what the range is; it just works.

For example

```julia
julia> v = AutoVector(0.0)        # default element, determining type

AutoVector{Float64}(1, 0, 0, 0.0, Float64[])

julia> v[-3] = pi

π = 3.1415926535897...

julia> v[10] = exp(1.0)

2.718281828459045

julia> doprint(v)

-3 3.141592653589793

-2 0.0

-1 0.0

0 0.0

1 0.0

2 0.0

3 0.0

4 0.0

5 0.0

6 0.0

7 0.0

8 0.0

9 0.0

10 2.718281828459045

julia> arange(v)

-3:10

julia> v[15]      # No resizing, give default

0.0

julia> arange(v)

-3:10
```

You can add and subtract AutoVectors (e.g. v+w, v-w), add or subtract a constant to all elements (v+a,v-a),  
and multiply by a constant, a\*v or v\*a, and v/a. Broadcasting (v .\* w) mostly doesn't work, since the lengths are
generally different, but there is a function pointmult(v,w) which does the same thing.
A number of other useful functions are also implemented, such as dot products and convolutions. 

For full documentation, see [documentation](https://srwhite59.github.io/AutoVectors.jl/)
