function [ T ] = calctrafos( dhparams, ndof, q )
% calculate trafos from DH parameters
  
  T=cell(1);
  T{1}=eye(4);
  for i=1:ndof
      if dhparams.types(i)== 0
         qakt=q(i);
         dakt=dhparams.d(i);
      else
         qakt=dhparams.q(i);
         dakt=q(i);
      end
      T{i+1}=T{i}*dhtrafo(dhparams.a(i),dhparams.alpha(i),qakt,dakt);
  end
  
end