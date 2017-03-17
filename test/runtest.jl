push!(LOAD_PATH, "/home/thatcher/Desktop/Projects")

using BiotSavart
using Base.Test

r = 0.01
point = Point(0.5,r,0)
line = Line(Point(0,0,0), Point(1,0,0))

(val, err) = BfromLine(point,line,1);
i_approx = μ₀/2pi/r

@test abs(val[3]-i_approx)<0.001
