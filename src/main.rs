// Newton's Method: given a function y(x) and dydx and some starting x0; approximate the value x'
// such that y(x') = 0 (zero-finding algorithm)
fn newton(xinit: f64, n: usize, tolerance: f64, f: impl Fn(f64) -> f64, dfdx: impl Fn(f64) -> f64, logging: bool) -> f64 {
  let mut xi = xinit;
  for i in 1..(n+1) {
    xi = xi - f(xi)/dfdx(xi);
    if logging {
      println!("-> i = {}; xi = {}; f(xi) = {}; error = {}", i, xi, f(xi), f(xi).abs()-tolerance);
    }
  }
  if f(xi).abs() < tolerance {
    println!("\nnewton(): xinit = {}; xi = {}; f(xi) = {}", xinit, xi, f(xi));
    xi
  } else {
    println!("\nnewton(): Could not find an xi such that f(xi) = 0 within desired tolerance of {}", tolerance);
    f64::MAX
  }
}
// START FROM CHATGPT
fn newton_raphson<F>(maxiters: usize, tolerance: f64, mut f: F, mut x0: f64) -> f64
where
    F: FnMut(f64) -> f64,
{
    for _ in 0..maxiters {
        let f_x0 = f(x0);
        if f_x0.abs() < tolerance {
            return x0;
        }
        let df_x0 = (f(x0 + tolerance) - f_x0) / tolerance;
        x0 = x0 - f_x0 / df_x0;
    }
    x0
}
fn backward_euler_known(n: usize, initval: f64, tfinal: f64, maxiters: usize, tolerance: f64, dydt: impl Fn(f64, f64) -> f64, solution: impl Fn(f64) -> f64, logging: bool) -> f64 {
    let k = tfinal / n as f64;
    let mut t = 0.0;
    let mut y = initval;
    for i in 0..n {
        let next_y = newton_raphson(maxiters, tolerance, |x| y + k*dydt(t+k, x) - x, y);
        y = next_y;
        if logging {
          println!("-> i = {}; t = {}; Estimated y(t) = {}; solution(t) = {}", i, t, y, solution(t));
        }
        t += k;
    }
    println!("\nbackward_euler_known(): Estimated y[tfinal] = {}; solution(tfinal)= {}; err = {}\n", y, solution(tfinal), (solution(tfinal)-y).abs());
    y
}
// Backward Euler: Yn+1 = Yn + k*f(tn+1, Yn+1) where f(t,y) = dy/dt for some y(t) and we know Y0 = initval
fn backward_euler(n: usize, initval: f64, tfinal: f64, maxiters: usize, tolerance: f64, dydt: impl Fn(f64, f64) -> f64, logging: bool) -> f64 {
    let k = tfinal / n as f64;
    let mut t = 0.0;
    let mut y = initval;
    for i in 0..n {
        let next_y = newton_raphson(maxiters, tolerance, |x| y + k*dydt(t+k, x) - x, y);
        y = next_y;
        if logging {
          println!("-> i = {}; t = {}; Estimated y(t) = {}", i, t, y);
        }
        t += k;
    }
    println!("\nbackward_euler(): Estimated y[tfinal] = {}\n", y);
    y
}
// END FROM CHATGPT
// Forward Euler: Yn = Yn + k*f(tn, Yn) where f(t,y) = dy/dt for some y(t) and we know Y0 = initval
fn forward_euler_known(n: usize, initval: f64, tfinal: f64, dydt: impl Fn(f64, f64) -> f64, y: impl Fn(f64) -> f64, logging: bool) -> f64 {
  let k = tfinal/(n as f64);
  let mut tn: Vec<f64> = vec![0.0; n+1];
  let mut yn: Vec<f64> = vec![0.0; n+1];
  let mut sln: Vec<f64> = vec![0.0; n+1];
  let mut err: Vec<f64> = vec![0.0; n+1];
  yn[0] = initval;
  sln[0] = initval;
  for i in 0..(n) {
    tn[i+1] = (i as f64)*k;
    let dydtval = dydt(tn[i], yn[i]);
    yn[i+1] = yn[i] + k*dydtval;
    sln[i+1] = y(tn[i+1]);
    err[i+1] = (sln[i+1]-yn[i+1]).abs();
    if logging {
      println!("-> i = {}; dydt(tn[{}], yn[{}]) = {}; tn[{}] = {}; yn[{}] = {}; sln[{}] = {}; err[{}] = {}", i, i, i, dydtval, i+1, tn[i+1], i+1, yn[i+1], i+1, sln[i+1], i+1, err[i+1]);
    }
  }
  println!("forward_euler_known(): POST ITERATION RESULT: k = {}; tn[{}] = {}; yn[{}] = {}; sln[{}] = {}; err[{}] = {}", k, n-1, tn[n-1], n-1, yn[n-1], n-1, sln[n-1], n-1, err[n-1]);
  yn[n]
}

fn forward_euler(n: usize, initval: f64, tfinal: f64, dydt: impl Fn(f64, f64) -> f64, logging: bool) -> f64 {
  let k = tfinal/(n as f64);
  let mut tn: Vec<f64> = vec![0.0; n+1];
  let mut yn: Vec<f64> = vec![0.0; n+1];
  yn[0] = initval;
  for i in 0..(n) {
    tn[i+1] = (i as f64)*k;
    let dydtval = dydt(tn[i], tn[i]);
    yn[i+1] = yn[i] + k*dydtval;
    if logging {
      println!("-> i = {}; dydt(tn[{}], yn[{}]) = {}; tn[{}] = {}; yn[{}] = {}", i, i, i, dydtval, i+1, tn[i+1], i+1, yn[i+1]);
    }
  }
  println!("forward_euler(): POST ITERATION RESULT: k = {} tn[{}] = {}; yn[{}] = {}", k, n-1, tn[n-1], n-1, yn[n-1]);
  yn[n]
}

fn main() {
  // Example Problem: f(t, y) = dy/dt = -y + (2e^(-t) * cos(2t)) where y(0) = 0
  //  -> the exact solution to this ODE/IVP is y(t) = e^(-t) * sin(2t)
  let n = 100;
  let initval = 0.0;
  let tfinal = 1.0;
  let dydt = |t: f64, y: f64| -> f64 {
    (-1.0*y) + (2.0*(-1.0*t).exp() * (2.0*t).cos())
  };
  let yt = |t: f64| -> f64 {
    (-1.0*t).exp() * (2.0*t).sin()
  };
  let fe = forward_euler(n, initval, tfinal, dydt, true);
  let fek = forward_euler_known(n, initval, tfinal, dydt, yt, true);
  let tolerance = 1e-12;
  let be = backward_euler(n, initval, tfinal, n, tolerance, dydt, true);
  let bek = backward_euler_known(n, initval, tfinal, n, tolerance, dydt, yt, true);
  println!("* Forward Euler trials:\n\n-> {} and {} (known)\n\n* Backward Euler trials:\n\n-> {} and {} (known)\n\n* Expected Solution = {}", fe, fek, be, bek, yt(tfinal));
  let fx = |x: f64| -> f64 {
    x.exp() - (x*x) + (3.0*x) - 2.0
  };
  let dfdx = |x: f64| -> f64 {
    x.exp() - (2.0*x) + 3.0
  };
  let xinit = -0.5;
  let xend = newton(xinit, n, tolerance, fx, dfdx, true);
  let result = fx(xend);
  println!("* Newton Method Trial:\n\n-> xinit = {}; xend = {}; f(xend) = {}; error = {}; tolerance = {}", xinit, xend, result, tolerance-result.abs(), tolerance);
}
