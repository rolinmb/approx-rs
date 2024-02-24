Implementing various numerical methods for mathematical approximations in rust.

* Central and Forward Finite Differences do not work due to numerical instability (floating point precision problems?); maybe need to use rust 'special' crate for
ln_gamma() function to replace current lnfactorial() function

* Dormand-Prince RK5 method performs worse than standard RK4 at the moment... may need to introduce adaptive step sizing?