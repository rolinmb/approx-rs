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
    if logging {
      println!("\nnewton(): xinit = {}; xi = {}; f(xi) = {}", xinit, xi, f(xi));
    }
    xi
  } else {
    println!("\nnewton(): Could not find an xi such that f(xi) = 0 within desired tolerance of {}", tolerance);
    f64::MAX
  }
}
// START FROM CHATGPT
// Newton-Raphson Method: another zero-finding algorithm likle Newton's method
fn newtonraphson<F>(maxiters: usize, tolerance: f64, mut f: F, mut x0: f64) -> f64
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
// Backward Euler's Method: Yn+1 = Yn + k*f(tn+1, Yn+1) where f(t,y) = dy/dt for some y(t) and we know Y0 = initval
fn backwardeuler(n: usize, initval: f64, tfinal: f64, maxiters: usize, tolerance: f64, dydt: impl Fn(f64, f64) -> f64, logging: bool) -> f64 {
  let k = tfinal / n as f64;
  let mut t = 0.0;
  let mut y = initval;
  for i in 0..n {
      let next_y = newtonraphson(maxiters, tolerance, |x| y + k*dydt(t+k, x) - x, y);
      y = next_y;
      if logging {
        println!("-> i = {}; t = {}; Estimated y(t) = {}", i, t, y);
      }
      t += k;
  }
  if logging {
    println!("backwardeuler(): Estimated y[tfinal] = {}", y);
  }
  y
}
fn backwardeuler_known(n: usize, initval: f64, tfinal: f64, maxiters: usize, tolerance: f64, dydt: impl Fn(f64, f64) -> f64, solution: impl Fn(f64) -> f64, logging: bool) -> f64 {
    let k = tfinal / n as f64;
    let mut t = 0.0;
    let mut y = initval;
    for i in 0..n {
        let next_y = newtonraphson(maxiters, tolerance, |x| y + k*dydt(t+k, x) - x, y);
        y = next_y;
        if logging {
          println!("-> i = {}; t = {}; Estimated y(t) = {}; solution(t) = {}", i, t, y, solution(t));
        }
        t += k;
    }
    if logging {
      println!("backwardeuler_known(): Estimated y[tfinal] = {}; solution(tfinal)= {}; err = {}", y, solution(tfinal), (solution(tfinal)-y).abs());
    }
    y
}
// END FROM CHATGPT
// Forward Euler's Method: Yn+1 = Yn + k*f(tn, Yn) where f(t,y) = dy/dt for some y(t) and we know Y0 = initval
fn forwardeuler(n: usize, initval: f64, tfinal: f64, dydt: impl Fn(f64, f64) -> f64, logging: bool) -> f64 {
  let k = tfinal/(n as f64);
  let mut tn: Vec<f64> = vec![0.0; n+1];
  let mut yn: Vec<f64> = vec![0.0; n+1];
  let mut dydtval: f64;
  yn[0] = initval;
  for i in 0..n {
    tn[i+1] = (i as f64)*k;
    dydtval = dydt(tn[i], yn[i]);
    yn[i+1] = yn[i] + k*dydtval;
    if logging {
      println!("-> i = {}; dydt(tn[{}], yn[{}]) = {}; tn[{}] = {}; yn[{}] = {}", i, i, i, dydtval, i+1, tn[i+1], i+1, yn[i+1]);
    }
  }
  if logging {
    println!("forwardeuler(): POST ITERATION RESULT: k = {} tn[{}] = {}; yn[{}] = {}", k, n-1, tn[n-1], n-1, yn[n-1]);
  }
  yn[n]
}

fn forwardeuler_known(n: usize, initval: f64, tfinal: f64, dydt: impl Fn(f64, f64) -> f64, y: impl Fn(f64) -> f64, logging: bool) -> f64 {
  let k = tfinal/(n as f64);
  let mut tn: Vec<f64> = vec![0.0; n+1];
  let mut yn: Vec<f64> = vec![0.0; n+1];
  let mut sln: Vec<f64> = vec![0.0; n+1];
  let mut err: Vec<f64> = vec![0.0; n+1];
  let mut dydtval: f64;
  yn[0] = initval;
  sln[0] = initval;
  for i in 0..n {
    tn[i+1] = (i as f64)*k;
    dydtval = dydt(tn[i], yn[i]);
    yn[i+1] = yn[i] + k*dydtval;
    sln[i+1] = y(tn[i+1]);
    err[i+1] = (sln[i+1]-yn[i+1]).abs();
    if logging {
      println!("-> i = {}; dydt(tn[{}], yn[{}]) = {}; tn[{}] = {}; yn[{}] = {}; sln[{}] = {}; err[{}] = {}", i, i, i, dydtval, i+1, tn[i+1], i+1, yn[i+1], i+1, sln[i+1], i+1, err[i+1]);
    }
  }
  if logging {
    println!("forwardeuler_known(): POST ITERATION RESULT: k = {}; tn[{}] = {}; yn[{}] = {}; sln[{}] = {}; err[{}] = {}", k, n-1, tn[n-1], n-1, yn[n-1], n-1, sln[n-1], n-1, err[n-1]);
  }
  yn[n]
}
// Modified Euler's Method (Heun's Method): Yn+1 = Yn + k/2 * (f(tn,Yn) + f(tn + k, Yn + k*f(tn,Yn))) 
fn modifiedeuler(n: usize, initval: f64, tfinal: f64, dydt: impl Fn(f64, f64) -> f64, logging: bool) -> f64 {
  let h = tfinal/(n as f64);
  let mut tn: Vec<f64> = vec![0.0; n+1];
  let mut yn: Vec<f64> = vec![0.0; n+1];
  let mut k1: f64;
  let mut k2: f64;
  yn[0] = initval;
  for i in 0..n {
    tn[i+1] = (i as f64)*h;
    k1 = dydt(tn[i], yn[i]);
    k2 = dydt(tn[i] + h, yn[i] + h*dydt(tn[i], yn[i]));
    yn[i+1] = yn[i] + (h/2.0)*(k1+k2);
    if logging {
      println!("-> i = {};  tn[{}] = {}; yn[{}] = {}", i, i+1, tn[i+1], i+1, yn[i+1]);
    }
  }
  if logging {
    println!("modifiedeuler(): POST ITERATION RESULT: h = {}; tn[{}] = {}; yn[{}] = {}", h, n-1, tn[n-1], n-1, yn[n-1]);
  }
  yn[n]
}

fn modifiedeuler_known(n: usize, initval: f64, tfinal: f64, dydt: impl Fn(f64, f64) -> f64, y: impl Fn(f64) -> f64, logging: bool) -> f64 {
  let h = tfinal/(n as f64);
  let mut tn: Vec<f64> = vec![0.0; n+1];
  let mut yn: Vec<f64> = vec![0.0; n+1];
  let mut sln: Vec<f64> = vec![0.0; n+1];
  let mut err: Vec<f64> = vec![0.0; n+1];
  let mut k1: f64;
  let mut k2: f64;
  yn[0] = initval;
  sln[0] = initval;
  for i in 0..n {
    tn[i+1] = (i as f64)*h;
    k1 = dydt(tn[i], yn[i]);
    k2 = dydt(tn[i] + h, yn[i] + h*dydt(tn[i], yn[i]));
    yn[i+1] = yn[i] + (h/2.0)*(k1+k2);
    sln[i+1] = y(tn[i+1]);
    err[i+1] = (sln[i+1]-yn[i+1]).abs();
    if logging {
      println!("-> i = {}; tn[{}] = {}; yn[{}] = {}; sln[{}] = {}; err[{}] = {}", i, i+1, tn[i+1], i+1, yn[i+1], i+1, sln[i+1], i+1, err[i+1]);
    }
  }
  if logging {
    println!("modifiedeuler_known(): POST ITERATION RESULT: h = {}; tn[{}] = {}; yn[{}] = {}; sln[{}] = {}; err[{}] = {}", h, n-1, tn[n-1], n-1, yn[n-1], n-1, sln[n-1], n-1, err[n-1]);
  }
  yn[n]
}
//Runge-Kutta Method(s): given a function f(t,y) = y'(t) for some unknown y(t) and y(tstart) = initval,
//we can approximate y(t) with even greater precision by 'modifying' the euler method further
/*fn rkfour(n: usize, initval: f64, tfinal: f64, dydt: impl Fn(f64, f64) -> f64, logging: bool) -> f64 {
  let h = tfinal/(n as f64);
  let mut tn: Vec<f64> = vec![0.0; n];
  let mut yn: Vec<f64> = vec![0.0; n];
  let mut k1: f64;
  let mut k2: f64;
  let mut k3: f64;
  let mut k4: f64;
  yn[0] = initval;
  for i in 0..n-1 {
    tn[i+1] = (i as f64)*h;
    k1 = dydt(tn[i], yn[i]);
    k2 = dydt(tn[i] + (h/2.0), yn[i] + (h*k1/2.0));
    k3 = dydt(tn[i] + (h/2.0), yn[i] + (h*k2/2.0));
    k4 = dydt(tn[i] + h, yn[i] + (h*k3));
    yn[i+1] = (h/6.0)*(k1 + (2.0*k2) + (2.0*k3) + k4);
    if logging {
      println!("-> i = {}; tn[{}] = {}; yn[{}] = {}", i, i+1, tn[i+1], i+1, yn[i+1]);
    }
  }
  if logging {
    println!("rkfour(): POST ITERATION RESULT: h = {}; tn[{}] = {}; yn[{}] = {}",  h, n-1, tn[n-1], n-1, yn[n-1]);
  }
  yn[n-1]
}

fn rkfour_known(n: usize, initval: f64, tfinal: f64, dydt: impl Fn(f64, f64) -> f64, y: impl Fn(f64) -> f64, logging: bool) -> f64 {
  let h = tfinal/(n as f64);
  let mut tn: Vec<f64> = vec![0.0; n];
  let mut yn: Vec<f64> = vec![0.0; n];
  let mut sln: Vec<f64> = vec![0.0; n];
  let mut err: Vec<f64> = vec![0.0; n];
  let mut k1: f64;
  let mut k2: f64;
  let mut k3: f64;
  let mut k4: f64;
  yn[0] = initval;
  sln[0] = initval;
  for i in 0..n-1 {
    tn[i+1] = (i as f64)*h;
    k1 = dydt(tn[i], yn[i]);
    k2 = dydt(tn[i] + (h/2.0), yn[i] + (h*k1/2.0));
    k3 = dydt(tn[i] + (h/2.0), yn[i] + (h*k2/2.0));
    k4 = dydt(tn[i] + h, yn[i] + (h*k3));
    yn[i+1] = (h/6.0)*(k1 + (2.0*k2) + (2.0*k3) + k4);
    sln[i+1] = y(tn[i+1]);
    err[i+1] = (sln[i+1]-yn[i+1]).abs();
    if logging {
      println!("-> i = {}; tn[{}] = {}; yn[{}] = {}; sln[{}] = {}; err[{}] = {}", i, i+1, tn[i+1], i+1, yn[i+1], i+1, sln[i+1], i+1, err[i+1]);
    }
  }
  if logging {
    println!("rkfour_known(): POST ITERATION RESULT: h = {}; tn[{}] = {}; yn[{}] = {}; sln[{}] = {}; err[{}] = {}",  h, n-1, tn[n-1], n-1, yn[n-1], n-1, sln[n-1], n-1, err[n-1]);
  }
  yn[n-1]
}*/

fn lnfactorial(n: usize) -> f64 {
	if n == 0 { return 0.0; }
	let mut result = 0.0;
	for i in 1..n+1 {
		result += (i as f64).log(std::f64::consts::E);
	}
	result
}

// Forward Finite Difference for approximation of desired order (nth) derivative evaluated at some point xinit
fn forwardfinitediff(order: usize, h: f64, xinit: f64, f: impl Fn(f64) -> f64, logging: bool) -> f64 {
  let mut ord = order;
  if order < 1 {
    println!("\nforwardfinitediff(): Order is less than 1 (order = {}); defaulting to order = 1\n", order);
    ord = 1;
  }
  let mut est = 0.0;
  for i in 0..=ord {
    let fi = i as f64;
    let coef = lnfactorial(ord) / (lnfactorial(i) * lnfactorial(ord - i)).exp();
    est += coef * f(xinit + fi * h);
    if logging {
      println!("-> i = {}; f(xinit+{}*h) = {}; est = {}", i, i, f(xinit+fi*h), est);
    }
  }
  if logging {
    println!("\ncentralfinitediff(): POST ITERATION RESULT: order = {}, h = {}, xinit = {}, est = {}\n", order, h, xinit, est);
  }
  est / f64::powf(h, ord as f64)
}

fn forwardfinitediff_known(order: usize, h: f64, xinit: f64, f: impl Fn(f64) -> f64, df: impl Fn(f64) -> f64, logging: bool) -> f64 {
  let mut ord = order;
  if order < 1 {
    println!("\nforwardfinitediff_known(): Order is less than 1 (order = {}); defaulting to order = 1\n", order);
    ord = 1;
  }
  let mut est = 0.0;
  for i in 0..=ord {
    let fi = i as f64;
    let coef = lnfactorial(ord) / (lnfactorial(i) * lnfactorial(ord - i)).exp();
    est += coef * f(xinit + fi * h);
    if logging {
      println!("-> i = {}; f(xinit+{}*h) = {}; est = {}; df(xinit+{}*h) = {}; err = {}", i, i, f(xinit+fi*h), est, i, df(xinit+fi*h), (df(xinit+fi*h)-est).abs());
    }
  }
  if logging {
    println!("\ncentralfinitediff_known(): POST ITERATION RESULT: order = {}, h = {}, xinit = {}, est = {}; df(xinit) = {}; err = {}\n", order, h, xinit, est, df(xinit), (df(xinit)-est).abs());
  }
  est / f64::powf(h, ord as f64)
}

// Central Finite Difference for approximation of desired order (nth) derivative evaluated at some point xinit
fn centralfinitediff(order: usize, h: f64, xinit: f64, f: impl Fn(f64) -> f64, logging: bool) -> f64 {
  let mut ord = order;
  if order < 1 {
    println!("\ncentralfinitediff(): Order is less than 1 (order = {}); defaulting to order = 1\n", order);
    ord = 1;
  }
  let mut est = 0.0;
  for i in 0..=ord {
    let fi = i as f64;
    let coef = lnfactorial(ord) / (lnfactorial(i) * lnfactorial(ord - i)).exp();
    est += coef * f(xinit + fi * h / 2.0) * if i % 2 == 0 { 1.0 } else { -1.0 };
    if logging {
      println!("-> i = {}; f(xinit+{}*h/2) = {}; est = {}", i, i, f(xinit + fi * h / 2.0), est);
    }
  }
  if logging {
    println!("\ncentralfinitediff(): POST ITERATION RESULT: order = {}, h = {}, xinit = {}, est = {}\n", order, h, xinit, est);
  }
  est / f64::powf(h, ord as f64)
}

fn centralfinitediff_known(order: usize, h: f64, xinit: f64, f: impl Fn(f64) -> f64, df: impl Fn(f64) -> f64, logging: bool) -> f64 {
  let mut ord = order;
  if order < 1 {
    println!("\ncentralfinitediff_known(): Order is less than 1 (order = {}); defaulting to order = 1\n", order);
    ord = 1;
  }
  let mut est = 0.0;
  for i in 0..=ord {
    let fi = i as f64;
    let coef = lnfactorial(ord) / (lnfactorial(i) * lnfactorial(ord - i)).exp();
    est += coef * f(xinit + fi * h / 2.0) * if i % 2 == 0 { 1.0 } else { -1.0 };
    if logging {
      println!("-> i = {}; f(xinit+{}*h/2) = {}; est = {}; df(xinit+{}*h/2) = {}; err = {}", i, i, f(xinit+fi*h/2.0), est, i, df(xinit+fi*h/2.0), (df(xinit+fi*h/2.0)-est).abs());
    }
  }
  if logging {
    println!("\ncentralfinitediff_known(): POST ITERATION RESULT: order = {}, h = {}, xinit = {}, est = {}; df(xinit) = {}; err = {}\n", order, h, xinit, est, df(xinit), (df(xinit)-est).abs());
  }
  est / f64::powf(h, ord as f64)
}

// Backwards Finite Difference for approximation of desired order (nth) derivative evaluated at some point xinit
fn backfinitediff(order: usize, h: f64, xinit: f64, f: impl Fn(f64) -> f64, logging: bool) -> f64 {
  let mut ord = order;
  if order < 1 {
    print!("\nbackfinitediff(): Order is less than 1 (order = {}); defaulting to order = 1\n", order);
    ord = 1;
  }
  let mut est = 0.0;
  for i in 0..=ord {
    let fi = i as f64;
    let coef = f64::powf(-1.0, fi) * (lnfactorial(ord) * lnfactorial(i) * lnfactorial(ord-i)).exp();
		est += coef * f(xinit - fi*h);
    if logging {
      print!("-> i = {}; f(xinit-{}*h) = {}; est = {}", i, i, f(xinit-fi*h), est);
    } 
  }
  if logging {
    print!("\nbackfinitediff(): POST ITERATION RESULT: order = {}, h = {}, xinit = {}, est = {}\n", order, h, xinit, est);
  }
  est / f64::powf(h, ord as f64)
}

fn backfinitediff_known(order: usize, h: f64, xinit: f64, f: impl Fn(f64) -> f64, df: impl Fn(f64) -> f64, logging: bool) -> f64 {
  let mut ord = order;
  if order < 1 {
    print!("\nbackfinitediff_known(): Order is less than 1 (order = {}); defaulting to order = 1\n", order);
    ord = 1;
  }
  let mut est = 0.0;
  for i in 0..=ord {
    let fi = i as f64;
    let coef = f64::powf(-1.0, fi) * (lnfactorial(ord) * lnfactorial(i) * lnfactorial(ord-i)).exp();
		est += coef * f(xinit-fi*h);
    if logging {
      print!("-> i = {}; f(xinit-{}*h) = {}; est = {}; df(xinit-{}*h) = {}; err = {}", i, i, f(xinit-fi*h), est, i, df(xinit-fi*h), (df(xinit-fi*h)-est).abs());
    }  
  }
  if logging {
    print!("\nbackfinitediff_known(): POST ITERATION RESULT: order = {}, h = {}, xinit = {}, est = {}; df(xinit) = {}; err = {}\n", order, h, xinit, est, df(xinit), (df(xinit)-est).abs());
  }
  est / f64::powf(h, ord as f64)
}
// START CHATGPT
struct ButcherTableau {
  a: Vec<Vec<f64>>,
  b: Vec<f64>,
  stages: usize,
}

fn general_rungekutta(f: impl Fn(f64, f64) -> f64, y0: f64, mut t0: f64, tfinal: f64, n: usize, tableau: &ButcherTableau) -> f64 {
  let h = (tfinal - t0) / (n as f64);
  let a = &tableau.a;
  let b = &tableau.b;
  let stages = tableau.stages;
  let mut y = y0;
  for _ in 0..n {
    let t = t0 + h;
    let mut k = 0.0;
    for i in 0..stages {
      let ti = t0 + a[i][0] * h;
      let yi = y;
      let mut tempyi = yi;
      for j in 1..stages {
        tempyi = tempyi + a[i][j] * h * k;
      }
      k = f(ti, tempyi);
    }
    y = y + h * b.iter().zip(std::iter::repeat(k)).map(|(bi, ki)| bi * ki).sum::<f64>();
    t0 = t;
  }
  y
}
/*END CHATGPT
Example ODE / IVP Problem: f(t, y) = dy/dt = -y + (2e^(-t) * cos(2t)) where y(t0 = 0.0) = 0.0
  -> the exact solution to this ODE/IVP is y(t) = e^(-t) * sin(2t)
*/
fn main() {
  let n = 1000000; // h (step size) = (tfinal - t0) / n
  let initval = 0.0;
  let tfinal = 1.0;
  let tolerance = 1e-10;
  println!("\nmain(): n = {}; initval = {}; tfinal = {}; tolerance = {}", n, initval, tfinal, tolerance);
  let dydt = |t: f64, y: f64| -> f64 {
    (-1.0*y) + (2.0*(-1.0*t).exp() * (2.0*t).cos())
  };
  let yt = |t: f64| -> f64 {
    (-1.0*t).exp() * (2.0*t).sin()
  };
  let fx = |x: f64| -> f64 {
    x.exp() - (x*x) + (3.0*x) - 2.0
  };
  let dfdx = |x: f64| -> f64 {
    x.exp() - (2.0*x) + 3.0
  };
  let xinit = -0.5;
  let nrend = newton(xinit, n, tolerance, fx, dfdx, false);
  let nresult = fx(nrend);
  println!("\n* Newton Method trial:\n\n > xinit = {}\n\n=> nrend = {}; f(nrend) = {}; error = {}; tolerance = {}", xinit, nrend, nresult, (tolerance-nresult).abs(), tolerance);
  let nrrend = newtonraphson(n, tolerance, fx, xinit);
  let nrresult = fx(nrrend);
  println!("\n* Newton-Raphson Method trial:\n\n > xinit = {}\n\n=> nrrend = {}; f(nrrend) = {}; error = {}; tolerance = {}\n", xinit, nrrend, nrresult, (tolerance-nrresult).abs(), tolerance);
  let fe = forwardeuler(n, initval, tfinal, dydt, false);
  let fek = forwardeuler_known(n, initval, tfinal, dydt, yt, false);
  let me = modifiedeuler(n, initval, tfinal, dydt, false);
  let mek = modifiedeuler_known(n, initval, tfinal, dydt, yt, false);
  let be = backwardeuler(n, initval, tfinal, n, tolerance, dydt, false);
  let bek = backwardeuler_known(n, initval, tfinal, n, tolerance, dydt, yt, false);
  println!("* Forward Euler Method trials:\n\n > {} and {} (known)\n\n* Backward Euler Method trials:\n\n > {} and {} (known)\n\n* Modified Euler Method trials:\n\n > {} and {} (known)\n\n=> Expected Solution = {}", fe, fek, be, bek, me, mek, yt(tfinal));
  //let rkf = rkfour(n, initval, tfinal, dydt, false);
  //let rkfe = rkfour_known(n, initval, tfinal , dydt, yt, false);
  //println!("\n* Runge-Kutta (4) trials:\n\n > {} and {} (known)\n\n=> Expected Solution = {}; error = {}", rkf, rkfe, yt(tfinal), (yt(tfinal)).abs());
  let step = 0.01;
  let ffd = forwardfinitediff(1, step, xinit, fx, false);
  let ffde = forwardfinitediff_known(1, step, xinit, fx, dfdx, false);
  println!("\n* Forward Finite Difference trials:\n\n > {} and {} (known); xinit = {}; Expected Solution = df(xinit) = {}; Error = {}\n", ffd, ffde, xinit, dfdx(xinit), (dfdx(xinit)-ffd).abs());
  let cfd = centralfinitediff(1, step, xinit, fx, false);
  let cfde = centralfinitediff_known(1, step, xinit, fx, dfdx, false);
  println!("\n* Central Finite Difference trials:\n\n > {} and {} (known); xinit = {}; Expected Solution = df(xinit) = {}; Error = {}\n", cfd, cfde, xinit, dfdx(xinit), (dfdx(xinit)-cfd).abs());
  let bfd = backfinitediff(1, step, xinit, fx, false);
  let bfde = backfinitediff_known(1, step, xinit, fx, dfdx, false);
  println!("\n* Backwards Finite Difference trials:\n\n > {} and {} (known); xinit = {}; Expected Solution = df(xinit) = {}; Error = {}\n", bfd, bfde, xinit, dfdx(xinit), (dfdx(xinit)-bfd).abs());
  let tableau_rk4 = ButcherTableau {
    a: vec![
      vec![0.0, 0.0, 0.0, 0.0],
      vec![0.5, 0.0, 0.0, 0.0],
      vec![0.0, 0.5, 0.0, 0.0],
      vec![0.0, 0.0, 1.0, 0.0],
    ],
    b: vec![1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0],
    stages: 4,
  };
  let tstart = 0.0;
  let mut yresult = general_rungekutta(dydt, initval, tstart, tfinal, n, &tableau_rk4);
  println!("\n* general_rungekutta (4):\n\n > Final value: {}; Expected Soultion = yt(tfinal) = {}; Error = {}", yresult, yt(tfinal), (yt(tfinal) - yresult).abs());
  let tableau_dp = ButcherTableau {
    a: vec![
      vec![0.0, 0.0, 0.0, 0.0, 0.0],
      vec![0.2, 0.0, 0.0, 0.0, 0.0],
      vec![3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0],
      vec![0.3, -0.9, 1.2, 0.0, 0.0],
      vec![-11.0 / 54.0, 2.5, -70.0 / 27.0, 35.0 / 27.0, 0.0],
      vec![1631.0 / 55296.0, 175.0 / 512.0, 575.0 / 13824.0, 44275.0 / 110592.0, 253.0 / 4096.0],
    ],
    b: vec![37.0 / 378.0, 0.0, 250.0 / 621.0, 125.0 / 594.0, 0.0, 512.0 / 1771.0],
    stages: 5,
  };
  yresult = general_rungekutta(dydt, initval, tstart, tfinal, n, &tableau_dp);
  println!("\n* general_rungekutta (Dormand-Prince):\n\n > Final value: {}; Expected Soultion = yt(tfinal) = {}; Error = {}", yresult, yt(tfinal), (yt(tfinal) - yresult).abs());
}
