function [x,y,z,s,w,info] = miip(A,b,c,ubounds)
% MIIP        - Main program for MIIP algorithm.

% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global x y z s w
global Ubounds_exist nub nt
global Dense_cols_exist idense ispars
global MONITOR_ON NNZA NNZL
global phi0 tau0
global Hist
global DISPLAY_LIPSOL
global MAXITER_LIPSOL

if DISPLAY_LIPSOL >= 1,
    disp('<<<<< This is MIIP algorithm >>>>>');
end

parameters;
[m,n] = size(A); nt = n + nub;
bnrm = norm(b); cnrm = norm(c); unrm = []; Hist = [];
if (Ubounds_exist) unrm = norm(ubounds); end

if Dense_cols_exist
    Aabs = abs(A(:,ispars));
else Aabs = abs(A); end;
symbfct(Aabs*Aabs'); clear Aabs;

[x,y,z,s,w] = initpoint(A,b,c,ubounds,bnrm,cnrm);
[Rxz,Rsw,dgap] = complementy(x,z,s,w);
[Rb,Rc,Ru,rb,rc,ru] = feasibility(A,b,c,ubounds);
[trerror,rrb,rrc,rru,rdgap,objp,objd] = ...
    errornobj(b,rb,bnrm,c,rc,cnrm,ubounds,ru,unrm,dgap);

iter = 0; converged = 0; mzeros = sparse(m,1); nzeros = sparse(n,1);

if DISPLAY_LIPSOL >= 1,
    fprintf('\n  Residuals:   Primal    Dual   U-bounds   Gap    TR_error\n');
    fprintf('  ---------------------------------------------------------\n');
end

Vd = [];
while iter <= MAXITER_LIPSOL  %%% Loop begins %%%

    if DISPLAY_LIPSOL == 2,
        monitor(iter,rb,rc,ru,dgap); 
    end;

    if DISPLAY_LIPSOL >= 1,
        fprintf('  Iter %4i:  ', iter);
        fprintf('%8.2e %8.2e %8.2e %8.2e %8.2e\n',rb,rc,ru,dgap,trerror);
    end

    if iter > 0
        [stop, converged] = stopping(tol);
        if stop break; end;
    end;

    [vmid,xn1,sn1] = getmisc(x,z,s,w);
    [P,U] = getpu(A,vmid);
    Vd = [Vd vmid];

    [dx,dy,dz,ds,dw] = ...
        direction(A,P,U,Rb,Rc,Ru,Rxz,Rsw,vmid,xn1,sn1,z,w,0,bnrm,1);
    [ap,ad] = ratiotest(dx,dz,ds,dw);

    if (tau0*ap < 1 | tau0*ad < 1)
        mu = centering(dx,dz,ds,dw,ap,ad,dgap,trerror);
        Rxz = dx.*dz; Rsw = ds.*dw;
        [dx2,dy2,dz2,ds2,dw2] = direction(A,P,U,mzeros,nzeros,nzeros,...
            Rxz,Rsw,vmid,xn1,sn1,z,w,mu,bnrm,2);
        dx = dx + dx2; dy = dy + dy2; dz = dz + dz2;
        ds = ds + ds2; dw = dw + dw2;
        [ap,ad] = ratiotest(dx,dz,ds,dw);
    end;

    [Rxz,Rsw,dgap] = update(ap,ad,dx,dy,dz,ds,dw,trerror,tol);
    [Rb,Rc,Ru,rb,rc,ru] = feasibility(A,b,c,ubounds);
    [trerror,rrb,rrc,rru,rdgap,objp,objd] = ...
        errornobj(b,rb,bnrm,c,rc,cnrm,ubounds,ru,unrm,dgap);

    Hist = [Hist [trerror rrb rrc rru rdgap objp objd]'];
    iter = iter + 1;
end;                  %%% Loop ends %%%
if iter >= MAXITER_LIPSOL,
    disp('Maximum number of iterations reached.');
end

info(1) = converged;
info(2) = iter;
info(3) = trerror;
%save vd Vd;
