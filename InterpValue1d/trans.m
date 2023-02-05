function s = trans(x, range1, range2)
s = range2(1,1) + (range2(1,2) - range2(1,1))/(range1(1,2) - range1(1,1))*(x-range1(1,1));
