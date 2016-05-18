using Cubature

μ₀ = 4pi*1e-7 # Tm/A
I = 1 # A

immutable Point
  x::Float64
  y::Float64
  z::Float64
end

immutable Line
   p0::Point
   p1::Point
end

point = Point(0.5,1,0)
line = Line(Point(0,0,0), Point(1,0,0))

x(t) = (1-t)line.p0.x + t*line.p1.x
y(t) = (1-t)line.p0.y + t*line.p1.y
z(t) = (1-t)line.p0.z + t*line.p1.z

r_(t) = [point.x, point.y, point.z] - [x(t), y(t), z(t)]
ds_(t) = [-t*line.p0.x + t*line.p1.x, -t*line.p0.y + t*line.p1.y, -t*line.p0.z + t*line.p1.z]

function dB(t,v)
  v[:] = μ₀*I/(4pi)*cross(ds_(t),r_(t))/norm(r_(t))^2
end

# (val, err) = hquadrature(3, (t,v) -> v[:] = μ₀*I/(4pi)*cross(ds_(t),r_(t))/norm(r_(t))^2, 0.0, 1.0, abstol=1e-8)
(val, err) = pquadrature(3, dB, 0.0, 1.0, abstol=1e-12)

println(val)

println(μ₀*I/(2pi)*0.5/sqrt(0.5^2 + 1))
