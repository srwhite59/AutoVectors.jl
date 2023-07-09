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
    function checklots()
	passedall = true
	v = AutoVector(0)
	passedall = passedall && v[1] == 0
	passedall = passedall && v[10] == 0
	passedall = passedall && avrange(v) == 1:0
	passedall = passedall && length(v) == 0
	passedall = passedall && mini(v) == 1
	passedall = passedall && maxi(v) == 0
	r = AutoVector(0.0)
	r[1] = 1.0
	r[2] = 2.0
	passedall = passedall && avrange(r) == 1:2
	passedall = passedall && avnorm(r) == sqrt(5.0)
	w = copy(r)
	z = convolve(r,w)
	passedall = passedall && avrange(z) == 2:4
	passedall = passedall && avvec(z) == [1.0,4.0,4.0]
	shz = applyshift(z,-2)
	passedall = passedall && avrange(shz) == 0:2
	revz = reverse_ind(z)
	passedall = passedall && avrange(revz) == -4:-2
	passedall = passedall && revz[-4] == z[4]
	z[8] = 1e-20
	passedall = passedall && avrange(z) == 2:8
	shrink!(z,1e-14)
	passedall = passedall && avrange(z) == 2:4
	clear!(z)
	passedall = passedall && length(z) == 0
	v = [1.0,2.0]
	vv = makeauto(v,offset=3)
	passedall = passedall && vv[-2] == 1.0
	vv = makeauto(v,firstindex=3)
	passedall = passedall && vv[3] == 1.0
	passedall
    end

    @test checkinit1() == 0.0
    @test checkinit2() == 0.0
    @test checkinit3() == sin(5/pi)
    @test checkinit4() == sin(5/pi)
    @test checklots() == true
end
