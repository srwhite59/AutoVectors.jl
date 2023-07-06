using AutoVectors
using Test

@testset "AutoVectors.jl" begin
    # Write your tests here.
    function checkinit1()
	v = AutoVector{Float64}(0.0,1,5,0,zeros(5))
	v[3]
    end
    function checkinit2()
	v = AutoVector()
	v[3]
    end
    function checkinit3()
	f = x->sin(x/pi)
	v = AutoVector(f,4,8)
	v[5]
    end
    function checkinit4()
	vs = [sin(x/pi) for x in 4:8]
	v = AutoVector(vs,4,8)
	v[5]
    end
    @test checkinit1() == 0.0
    @test checkinit2() == 0.0
    @test checkinit3() == sin(5/pi)
    @test checkinit4() == sin(5/pi)
end
