fn main() {
  // Forward Euler Method for example equation:
  // f(t, y) = dy/dt = -y + (2e^(-t) * cos(2t)) where y(0) = 0
  //  -> the exact solution y(t) = e^(-t) * sin(2t)
  let T = 1.0;
  let N = 5;
  let k = T/N;
  let mut t = vec![0.0, N+1];
  let mut y = vec![0.0, N+1];
  let mut sln = vec![0.0, N+1];
  let mut err = vec![0.0, N+1];
  y[0] = 0.0;
  sln[0] = 0.0;
  for i in range 0..N {
    t[i+1] = i*k;
    y[i+1] = y[i]+k*(-1.0*y[i]+2.0*(-1.0*t[i]).exp()*(2.0*t[i]).cos());
    sln[i+1] = ((-1.0*t[i+1]).exp()*(2.0*t[i+1]).sin());
    err[i+1] = (sln[i+1]-y[i+1]).abs();
  }
}
