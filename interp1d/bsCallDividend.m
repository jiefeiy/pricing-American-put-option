function V = bsCallDividend(Price,Strike,Rate,Time,Expiration,Volatility,DivYield)
theta = Expiration - Time;
d1 = ( log(Price./Strike) + (Rate - DivYield +Volatility.^2).*theta )./(Volatility.*sqrt(theta));
d2 = d1 - Volatility.*sqrt(theta);
V = Price.*exp(-DivYield.*theta).*normcdf(d1) - Strike.*exp(-Rate.*theta).*normcdf(d2);