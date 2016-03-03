function [ T ] = dhtrafo( a, alpha, q, d )
% produce dh trafo
  
  cq=cos(q);
  sq=sin(q);
  ca=cos(alpha);
  sa=sin(alpha);
  
  T=[cq,    -sq,    0,      a;
     sq*ca, cq*ca,  -sa,    -sa*d;
     sq*sa, cq*sa,  ca,     ca*d;
     0,     0,      0,      1];
  
end