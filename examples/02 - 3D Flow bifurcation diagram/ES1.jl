# Jafari, Sajad, et al. "Simple chaotic 3D flows with surfaces of equilibria." Nonlinear Dynamics 86 (2016): 1349-1358.
function ES1!(du, u, p, t)
  x, y, z = u
  a = p[1]
  f = x
  du[1] = f*y
  du[2] = f*z
  du[3] = f*(-x+a*y^2-x*z)
end