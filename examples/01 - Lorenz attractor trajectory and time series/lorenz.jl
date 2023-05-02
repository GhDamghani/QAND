function lorenz!(du, u, p, t)
  x, y, z = u
  sigma, rho, beta = p
  du[1] = sigma * (y - x)
  du[2] = x * (rho - z) - y
  du[3] = x * y - beta * z
end