function coef = nf_DH(odefile,jacobian,hessians,der3,x,p,nphase)
%
% coef = nf_DH(odefile,jacobian,hessians,der3,x,p,nphase)
% compute normal form coefficients for double-hopf.
%
global cds
  jac=cjac(odefile,jacobian,x,p,ones(length(p),1));
  [X,D] = eig(jac);
  index=find(real(diag(D))<1e-6 & sign(imag(diag(D)))==1);% This should give a 1x2 vector
if(size(index)==1)					  % Otherwise there is a neutral saddle involved.
  debug('Neutral saddle\n');
  coef =[0 0;0 0];
  return;
end
  ev0 = diag(D(index(1),index(1)));
  ev1 = diag(D(index(2),index(2)));
  q0 = X(:,index(1));
  q1 = X(:,index(2));
  [X,DD] = eig(jac');
  K = diag(D);
  index2=find(conjugate(diag(DD))==K(index));
  qad0 = X(:,index2(1));
  qad1 = X(:,index2(2));
  p0=qad0/(q0'*qad0);
  p1=qad1/(q1'*qad1);
  hessIncrement = (cds.options.Increment)^(3.0/4.0);
  ten3Increment = (cds.options.Increment)^(3.0/5.0);
if (cds.options.SymDerivative >= 3)
  hess = chess(odefile,jacobian,hessians,x,p,ones(length(p),1));
  tens = ctens3(odefile,jacobian,hessians,der3,x,p,ones(length(p),1));
else
  hess = [];
  tens = [];
end
%2nd order vectors
  h2000 = (2*ev0*eye(nphase)-jac)\multilinear2(odefile,hess,q0,q0,x,p,hessIncrement);		% (2iw_0-A)\B(q0,q0)
  h1100 = -jac\multilinear2(odefile,hess,q0,conj(q0),x,p,hessIncrement);			% -A\B(q0,conj(q0))
  h1010 = ((ev0+ev1)*eye(nphase)-jac)\multilinear2(odefile,hess,q0,q1,x,p,hessIncrement);	% (i(w_0+w_1)-A)\B(q0,q1)
  h1001 = ((ev0-ev1)*eye(nphase)-jac)\multilinear2(odefile,hess,q0,conj(q1),x,p,hessIncrement);	% (i(w_0-w_1)-A)\B(q0,conj(q1))
  h0020 = (2*ev1*eye(nphase)-jac)\multilinear2(odefile,hess,q1,q1,x,p,hessIncrement);		% (2iw_1-A)\B(q1,q1)
  h0011 = -jac\multilinear2(odefile,hess,q1,conj(q1),x,p,hessIncrement);			% -A\B(q1,conj(q1))
%3rd order vectors
  h2100 = multilinear3(odefile,tens,q0,q0,conj(q0),x,p,ten3Increment);				%  C(q0,q0,conj(q0))
  h2100 = h2100 + 2*multilinear2(odefile,hess,q0,h1100,x,p,hessIncrement);			%+2B(h1100,q0)
  h2100 = h2100 + multilinear2(odefile,hess,h2000,conj(q0),x,p,hessIncrement);			%+ B(h2000,conj(q0))  
  h1011 = multilinear3(odefile,tens,q0,q1,conj(q1),x,p,ten3Increment);				%  C(q0,q1,conj(q1))
  h1011 = h1011 + multilinear2(odefile,hess,h0011,q0,x,p,hessIncrement);			%+ B(q0,h0011)
  h1011 = h1011 + multilinear2(odefile,hess,h1001,q1,x,p,hessIncrement);			%+ B(q1,h1001)
  h1011 = h1011 + multilinear2(odefile,hess,h1010,conj(q1),x,p,hessIncrement);			%+ B(conj(q1),h1010)
  h1110 = multilinear3(odefile,tens,q1,q0,conj(q0),x,p,ten3Increment);				%  C(q1,q0,conj(q0))
  h1110 = h1011 + multilinear2(odefile,hess,h1100,q1,x,p,hessIncrement);			%+ B(q1,h1100)
  h1110 = h1011 + multilinear2(odefile,hess,conj(h1001),q0,x,p,hessIncrement);			%+ B(q0,conj(h1001))
  h1110 = h1011 + multilinear2(odefile,hess,h1010,conj(q0),x,p,hessIncrement);			%+ B(conj(q0),h1010)
  h0021 = multilinear3(odefile,tens,q1,q1,conj(q1),x,p,ten3Increment);				%  C(q1,q1,conj(q1))
  h0021 = h0021 + 2*multilinear2(odefile,hess,q1,h0011,x,p,hessIncrement);			%+2B(h0011,q1)
  h0021 = h0021 + multilinear2(odefile,hess,h0020,conj(q1),x,p,hessIncrement);			%+ B(h0020,conj(q1))
%coefficients
  coef = [ p0'*h2100/2.0 p0'*h1011; p1'*h1101 p1'*h0021/2.0 ];
