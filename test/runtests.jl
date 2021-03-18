using BiotSavart
using Test

r = 0.01
point = Point(0.5,r,0)
line = Line(Point(0,0,0), Point(1,0,0))

(val, err) = BfromLine(point,line,1);
i_approx = μ₀/2pi/r
@testset "BfromLine" begin
  @test abs(val[3]-i_approx)<0.001
end


#Tests for BfromLines
N=200
theta = (1:N) * 2pi/N
x = cos.(theta)
y = sin.(theta)
z = zeros(200)

points = Point.(x,y,z)
lines = Line.(points[1:N-1],points[2:N])
lines = vcat(lines, Line(points[N],points[1]))

x = x/2
y = y/2

points2 = Point.(x,y,z)
lines2 = Line.(points2[1:N-1],points2[2:N])
lines2 = vcat(lines2, Line(points2[N],points2[1]))

B = BfromLines(Point(0,0,0),lines,1)
BtwoAmps = BfromLines(Point(0,0,0),lines,2)
Bhalf = BfromLines(Point(0,0,0),lines2,1)
bCalc = μ₀/2
bCalcHalf = μ₀
error = (bCalc-B[3])/bCalc
error2 = (bCalcHalf-Bhalf[3])/bCalcHalf
error3 = (2*bCalc-BtwoAmps[3])/(2*bCalc)

#more tests to be made from this: http://hyperphysics.phy-astr.gsu.edu/hbase/magnetic/curloo.html

@testset "BfromLines" begin
  @testset "1 meter coil" begin
    @test abs(error)<0.01
  end
  @testset "0.5 meter coil" begin
    @test abs(error2)<0.01
  end
  @testset "Amperage" begin
    @test abs(error3)<0.001
  end
end
