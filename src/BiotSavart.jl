module BiotSavart
export Point, Line, BfromLine, BfromLines, μ₀

using Cubature

μ₀ = 4pi*1e-7 # Tm/A

r = 0.01

immutable Point
  x::Float64
  y::Float64
  z::Float64
end

immutable Line
   p0::Point
   p1::Point
end

function BfromLine(r,line,I)
  x(t) = (1-t)line.p0.x + t*line.p1.x
  y(t) = (1-t)line.p0.y + t*line.p1.y
  z(t) = (1-t)line.p0.z + t*line.p1.z

  r_(t)  = [r.x, r.y, r.z] - [x(t), y(t), z(t)]

  #ds might be incorrect - check the docs
  ds_(t) = [line.p1.x-line.p0.x, line.p1.y-line.p0.y, -line.p1.z-line.p0.z]

  function dB(t,v)
    v[:] = μ₀*I/(4pi)*cross(ds_(t),r_(t)) / norm(r_(t))^3
  end

  # (val, err) = hquadrature(3, (t,v) -> v[:] = μ₀*I/(4pi)*cross(ds_(t),r_(t))/norm(r_(t))^2, 0.0, 1.0, abstol=1e-8)
  pquadrature(3, dB, 0.0, 1.0, abstol=1e-12)
end

function BfromLines(r,lines,I)
  function f(l)
    (val, _) = BfromLine(r,l,I)
    val
  end
  reduce(+, map(f, lines))
end

end
