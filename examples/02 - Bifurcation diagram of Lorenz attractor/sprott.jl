# Sprott, Julien C., et al. "Megastability: Coexistence of a countable infinity of nested attractors in a periodically-forced oscillator with spatially-periodic damping." The European Physical Journal Special Topics 226 (2017): 1979-1985.
function sprott(du, u, p, t)
  x, v = u
  omega, A, Omega = p
  du[1] = v
  du[2] = -omega^2*x+v*cos(x)+A*sin(Omega*t)
end