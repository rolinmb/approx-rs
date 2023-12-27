fn main() {
  // Forward Euler Method for example equation:
  // f(t, y) = dy/dt = -y + (2e^(-t) * cos(2t)) where y(0) = 0
  //  -> the exact solution y(t) = e^(-t) * sin(2t)
  let t = 1.0;
  let n = 50;
  let k = t/(n as f64);
  let mut tn: Vec<f64> = vec![0.0; n];
  let mut yn: Vec<f64> = vec![0.0; n];
  let mut solution: Vec<f64> = vec![0.0; n];
  let mut err: Vec<f64> = vec![0.0; n];
  //yn[0] = 0.0; // technically are already 0.0
  //solution[0] = 0.0;
  print!("\nBEFORE ITERATING: tn[0] = 0.0; yn[0] = 0.0; solution[0] = 0.0; err[0] = 0.0\n\n");
  for i in 0..(n-1) {
    tn[i+1] = (i as f64)*k;
    yn[i+1] = yn[i] + k*(-1.0*yn[i] + 2.0*(-1.0*tn[i]).exp() * (2.0*tn[i]).cos() );
    solution[i+1] = (-1.0*tn[i+1]).exp()*(2.0*tn[i+1]).sin();
    err[i+1] = (solution[i+1]-yn[i+1]).abs();
    print!("> i = {}; tn[{}] = {}; yn[{}] = {}; solution[{}] = {}; err[{}] = {}\n", i, i+1, tn[i+1], i+1, yn[i+1], i+1, solution[i+1], i+1, err[i+1]);
  }
  print!("\nPOST ITERATION RESULT: tn[{}] = {}; yn[{}] = {}; solution[{}] = {}; err[{}] = {}\n", n-1, tn[n-1], n-1, yn[n-1], n-1, solution[n-1], n-1, err[n-1]);
}
