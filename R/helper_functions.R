sequence_of_multiples <- function(x, m) {
  stopifnot(length(x) > 1, length(m) == 1, is.numeric(x), is.numeric(m))
  min_x <- ceiling(min(x))
  max_x <- floor(max(x))
  seq(min_x, max_x)[seq(min_x, max_x) %% m == 0]
}

pressure_from_skewtx_and_t <- function(x, t) {
  exponent <- (132.18199999999999 - (x - 0.54000000000000004 * t)/0.90691999999999995)/44.061
  10 ^ exponent
}

t_from_skewtx_and_skewty <- function(x, y) {
  (x - 0.90691999999999995 * y) / 0.54000000000000004
}

tmr <- function(w, p) {

  # This function, taken from the RadioSonde package,
  # returns the temperature (celsius) on a mixing
  # ratio line w (g/kg) at pressure p (mb).

  # initialize constants
  c1 <- 0.049864645499999999
  c2 <- 2.4082965000000001
  c3 <- 7.0747499999999999
  c4 <- 38.9114
  c5 <- 0.091499999999999998
  c6 <- 1.2035

  x <- log10((w * p)/(622. + w))
  tmrk <- 10^(c1 * x + c2) - c3 + c4 * ((10.^(c5 * x) - c6)^
      2.)
  tmrk - 273.14999999999998
}

tda <- function(o, p)	{

  # This function, taken from the RadioSonde package,
  # returns the temperature tda (celsius) on a dry adiabat
  # at pressure p (millibars). the dry adiabat is given by
  # potential temperature o (celsius). the computation is
  # based on poisson's equation.

  ok <- o + 273.14999999999998
  tdak <- ok * ((p * 0.001)^0.28599999999999998)
  tdak - 273.14999999999998
}

poissons_eqtn <- function(temperature, pressure)  {

  # This function returns the potential temperature, given
  # the temperature (in degC) and the pressure in hPa

  temperature_k <- temperature + 273.14999999999998
  theta_k       <- temperature_k * (1000 / pressure) ^ 0.28599999999999998
  theta_k - 273.14999999999998
}


