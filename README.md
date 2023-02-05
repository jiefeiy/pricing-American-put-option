# prcing-American-put-option
Pricing American option with a single asset.

Run '/pathGen/test_pathBS1d.m' to 
generate asset price paths (or anti paths for variance reduction) under the Black Scholes model.

Run '/InterpValue1d/main3chebfun.m' to implement the dynamic Chebyshev method. This needs the toolbox 'chebfun', see https://www.chebfun.org/.

Run '/interp1d/AmerPut1D.m' to price American put by interpolation and CC quadrature.

Run '/interp1d/AmerCall1D.m' to price American call with dividend yield by interpolation and CC quadrature.

Run '/LSMC/test_LSM1d.m' to price American put by least square Monte Carlo (LSMC), see Longstaff and Schwartz. 
