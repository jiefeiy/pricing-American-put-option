function V = bsPutDividend(Price,Strike,Rate,Time,Expiration,Volatility,DivYield)
theta = Expiration - Time;
d1 = ( log(Price./Strike) + (Rate - DivYield +Volatility.^2).*theta )./(Volatility.*sqrt(theta));
d2 = d1 - Volatility.*sqrt(theta);
V = Strike.*exp(-Rate.*theta).*normcdf(-d2) - Price.*exp(-DivYield.*theta).*normcdf(-d1);