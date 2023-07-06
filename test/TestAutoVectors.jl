using AutoVectors
using DSP

function testAutoVectors()
    v = AutoVector{Float64}(0.0,1,5,0,zeros(5))
    v = [exp(-0.1*i^2) for i=-49:51]
    @time conv(v,v)
    @time conv(v,v)
    v = [exp(-0.1*i^2) for i=-63:64]
    @time conv(v,v)
    @time conv(v,v)
    v = [exp(-0.1*i^2) for i=-64:64]
    @time conv(v,v)
    @time conv(v,v)
    #exit(0)
    v = AutoVector(i->exp(-i*1.0),-1,5)
    doprint(v)
    v *= 2.0
    doprint(v)
    v *= 0.5

    vv = [v,v,v]
    vv0 = [1.0,2.0,3.0] * vv - 6 * v
    test0 = avdot(vv0,vv0)
    @show test0, "should be 0"

    vb = AutoVector(i->exp(-i*i*0.2),-30,30)
    shrink!(vb,1.0e-10)
    println("vb:")
    doprint(vb)
    println()

    vvb1 = pointmult(v,vb)
    println("v:")
    doprint(v)
    println("pointmult(v,vb):")
    doprint(vvb1)
    #vvb2 = broadcast(*,v,vb)
    #vvb3 = v .* vb
    #doprint(vvb3-vvb1)

    w = AutoVector([1.0,2.0],4,5)
    #doprint(w)
    #return
    for test =1:30
	u = AutoVector(0.0)
	g = AutoVector(0.0)
	v = AutoVector(0.0)
	ua = rand(-1000:1000)
	ub = ua+rand(0:5000)
	ga = rand(-20:20)
	gb = ga+rand(0:20)
	va = rand(-500:500)
	vb = va+rand(0:5000)
	for i=ua:ub
	    u[i] = rand()
	end
	for i=ga:gb
	    g[i] = rand()
	end
	for i=va:vb
	    v[i] = rand()
	end
	res = avdot(convolve(u,g),v)
	res2 = avtripconv(u,g,v)
	if abs(res-res2)/abs(res) > 1.0e-10
	    @show ua,ub,ga,gb,va,vb
	    @show res,res2
	end
	if abs(res-res2)/abs(res) > 1.0e-10
	    @show ua,ub,ga,gb,va,vb
	    @show res,res2
	    @time avdot(convolve(u,g),v)
	    #@code_warntype avtripconv(u,g,v)
	    @time avtripconv(u,g,v)
	    println(" ")
	end
    end
end

testAutoVectors()

#using ProfileView
##@code_warntype testAutoVectors()
#testAutoVectors()
#@profile testAutoVectors()
#ProfileView.view()

