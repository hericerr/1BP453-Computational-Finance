function [fmin, xmin] = goldenRaioMin(f, a, b, TOL)

r = (sqrt(5)-1)/2;
leftMid = b - r*(b-a);
rightMid = a + r*(b-a);

if (abs(b-a) <=  TOL)
   xmin = (a+b)/2;
   fmin = f(xmin);
else
   if (f(leftMid) > f(rightMid))
       [fmin, xmin] = goldenRaioMin(f, leftMid, b, TOL);       
   else
       [fmin, xmin] = goldenRaioMin(f, a, rightMid, TOL);
   end
end    

end