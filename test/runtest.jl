push!(LOAD_PATH, "/home/thatcher/Desktop/Projects")

using BiotSavart
using Base.Test

r = 0.01
point = Point(0.5,r,0)
line = Line(Point(0,0,0), Point(1,0,0))

(val, err) = BfromLine(point,line,1);
i_approx = μ₀/2pi/r

@test abs(val[3]-i_approx)<0.001


#Tests for BfromLines
N=200
theta = linspace(0,2pi,N)
x = cos(theta)
y = sin(theta)
z = zeros(200)

points = Point.(x,y,z)
lines = Line.(points[1:N-1],points[2:N])
lines = vcat(lines, Line(points[N],points[1]))

B = BfromLines(Point(0,0,0),lines,1)
bCalc = μ₀/2
println(B)
println(bCalc)
@test abs(B[3]-bCalc)<0.00001
