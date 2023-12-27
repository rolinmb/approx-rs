// Forward Euler: Yn = Yn + k*f(tn, Yn) where f(t,y) = dy/dt and we know f(0, Y0) = initval
fn forward_euler_known(n: usize, initval: f64, tfinal: f64, dydt: impl Fn(f64, f64) -> f64, y: impl Fn(f64) -> f64, logging: bool) {
  let k = tfinal/(n as f64);
  let mut tn: Vec<f64> = vec![0.0; n+1];
  let mut yn: Vec<f64> = vec![0.0; n+1];
  let mut sln: Vec<f64> = vec![0.0; n+1];
  let mut err: Vec<f64> = vec![0.0; n+1];
  yn[0] = initval;
  sln[0] = initval;
  print!("\nBEFORE ITERATING: k = {}; tn[0] = 0.0; yn[0] = 0.0; solution[0] = 0.0; err[0] = 0.0\n", k);
  for i in 0..(n) {
    tn[i+1] = (i as f64)*k;
    let dydtn = dydt(tn[i], yn[i]);
    yn[i+1] = yn[i] + k*dydtn;
    sln[i+1] = y(tn[i+1]);
    err[i+1] = (sln[i+1]-yn[i+1]).abs();
    if logging {
      print!("-> i = {}; dydt(tn[{}], yn[{}]) = {}; tn[{}] = {}; yn[{}] = {}; sln[{}] = {}; err[{}] = {}\n", i, i, i, dydtn, i+1, tn[i+1], i+1, yn[i+1], i+1, sln[i+1], i+1, err[i+1]);
    }
  }
  print!("\nPOST ITERATION RESULT: k = {}; tn[{}] = {}; yn[{}] = {}; sln[{}] = {}; err[{}] = {}\n", k, n-1, tn[n-1], n-1, yn[n-1], n-1, sln[n-1], n-1, err[n-1]);
}

fn forward_euler(n: usize, initval: f64, tfinal: f64, dydt: impl Fn(f64, f64) -> f64, logging: bool) {
  let k = tfinal/(n as f64);
  let mut tn: Vec<f64> = vec![0.0; n+1];
  let mut yn: Vec<f64> = vec![0.0; n+1];
  yn[0] = initval;
  print!("\nBEFORE ITERATING: k = {}; tn[0] = 0.0; yn[0] = 0.0\n", k);
  for i in 0..(n) {
    tn[i+1] = (i as f64)*k;
    let dydtn = dydt(tn[i], tn[i]);
    yn[i+1] = yn[i] + k*dydtn;
    if logging {
      print!("-> i = {}; dydt(tn[{}], yn[{}]) = {}; tn[{}] = {}; yn[{}] = {}\n", i, i, i, dydtn, i+1, tn[i+1], i+1, yn[i+1]);
    }
  }
  print!("\nPOST ITERATION RESULT: k = {} tn[{}] = {}; yn[{}] = {}\n", k, n-1, tn[n-1], n-1, yn[n-1]);
}

fn main() {
  // Example Problem: f(t, y) = dy/dt = -y + (2e^(-t) * cos(2t)) where y(0) = 0
  //  -> the exact solution y(t) = e^(-t) * sin(2t)
  let n = 100;
  let initval = 0.0;
  let tfinal = 1.0;
  let dydt = |t: f64, y: f64| -> f64 {
    (-1.0*y) + (2.0*(-1.0*t).exp() * (2.0*t).cos())
  };
  let yt = |t: f64| -> f64 {
    (-1.0*t).exp() * (2.0*t).sin()
  };
  forward_euler(n, initval, tfinal, dydt);
  forward_euler_known(n, initval, tfinal, dydt, yt);
}
