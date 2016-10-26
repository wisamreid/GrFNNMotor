function [A,b,c,lbounds,ubounds,FEASIBLE] = preprocess(A,b,c,lbounds,ubounds)
% PREPROCESS  - Preprocessing input data.
% Usage:  [A,b,c,lbounds,ubounds] = preprocess(A,b,c,lbounds,ubounds)

% Yin Zhang, last updated April, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global NNZA
global b_orig c_orig
global data_changed
global Lbounds_non0
global Ubounds_exist nub
global Fixed_exist ifix infx xfix
global Zrcols_exist izrcol inzcol xzrcol
global Sgtons_exist isolved insolved xsolved
global DISPLAY_LIPSOL


% Big number
BIG = 1.0000e+032;

if DISPLAY_LIPSOL >= 1,
    disp('Preprocessing ...');
end
[m, n] = size(A); FEASIBLE = 1;

if any(lbounds > ubounds)
    error('Preprocessor: Lower bound exceeds upper bound');
end;

if ~issparse(A) A = sparse(A); end;
b = sparse(b); c = sparse(c);
lbounds = sparse(lbounds);
b_orig = b; c_orig = c;
data_changed = 0;

%----- delete fixed variables -----
fixed = lbounds == ubounds;
Fixed_exist = any(fixed);
if (Fixed_exist)
    ifix = find(fixed);
    infx = find(1 - fixed);
    xfix = lbounds(ifix);
    c    = c(infx);
    b    = b - A(:,ifix)*sparse(xfix);
    A    = A(:,infx);
    lbounds = lbounds(infx);
    ubounds = ubounds(infx);
    data_changed = 1;
end

%----- delete zero rows -----
rnnzct = sum(spones(A'));
if any(rnnzct == 0)
    izrows = find(rnnzct == 0);
    if any(b(izrows) ~= 0)
        disp('Preprocessor: problem infeasible');
        FEASIBLE = 0;
        return;
    end;
    inzrows = find(rnnzct > 0);
    A = A(inzrows,:); b = b(inzrows);
    rnnzct = rnnzct(inzrows);

    if DISPLAY_LIPSOL >= 1,
        disp(sprintf(' (%i 0-rows)', length(izrows)));
    end
    data_changed = 1; [m, n] = size(A);
end

%----- make A structurally "full rank" -----
sprk = sprank(A');
if (sprk < m)
    Row_deleted = 1;
    [dmp, tmp] = dmperm(A);
    irow = dmp(1:sprk);
    A = A(irow,:); b = b(irow);
    rnnzct = rnnzct(irow);
    if DISPLAY_LIPSOL >= 1,
        fprintf(' (%i dep-rows)', m-sprk);
    end
    data_changed = 1;
end

%----- delete zero columns -----
Zrcols_exist = 0;
zrcol = (max(abs(A)) == 0)';
if any(zrcol == 1)
    Zrcols_exist = 1;
    izrcol = find(zrcol);
    if any(c(izrcol) < 0 & ubounds(izrcol) > BIG-1)
        disp('Preprocessor: problem unbounded below');
        FEASIBLE = 0; 
        return;
    end
    xzrcol = zeros(size(izrcol))...
        + (c(izrcol) < 0).*ubounds(izrcol)...
        + (c(izrcol) > 0).*lbounds(izrcol);
    inzcol = find(1 - zrcol);
    A = A(:,inzcol);
    c = c(inzcol);
    lbounds = lbounds(inzcol);
    ubounds = ubounds(inzcol);
    if DISPLAY_LIPSOL >= 1,
        fprintf(' (%i 0-columns)', nnz(zrcol));
    end
    data_changed = 1;
end

%----- solve singleton rows -----
Sgtons_exist = 0;
singleton = (rnnzct == 1);
nsgrows = nnz(singleton);
if nsgrows >= max(1, .01*size(A,1))
    Sgtons_exist = 1;
    isgrows = find(singleton);
    iothers = find(1 - singleton);
    if DISPLAY_LIPSOL >= 1,
        fprintf(' (%i singletons)',nsgrows);
    end

    Atmp = A(isgrows,:); Atmp1 = spones(Atmp); btmp = b(isgrows);
    if nsgrows == 1
        isolved  = find(Atmp1);
        insolved = find(Atmp1 == 0);
        xsolved  = b(isgrows)/Atmp(isolved);
    else
        colnnzct = sum(Atmp1);
        isolved  = find(colnnzct);
        insolved = find(colnnzct == 0);
        [ii, jj] = find(Atmp);
        Atmp = Atmp(ii,jj); btmp = btmp(ii);
        xsolved  = btmp./diag(Atmp);
        if any(colnnzct >  1)
            repeat = diff([0; jj]) == 0;
            for i = 1:length(xsolved) - 1
                if repeat(i+1) & xsolved(i+1) ~= xsolved(i)
                    disp('Preprocessor: problem infeasible');
                    FEASIBLE = 0; 
                    return;
                end;
            end;
            ii = find(~repeat); jj = ii;
            Atmp = Atmp(ii,jj); btmp = btmp(ii);
            xsolved  = btmp./diag(Atmp);
        end;
    end;

    if any(xsolved < lbounds(isolved)) | any(xsolved > ubounds(isolved))
        disp('Preprocessor: problem infeasible');
        FEASIBLE = 0; 
        return;
    end;

    b = b(iothers) - A(iothers,isolved)*xsolved;
    A = A(iothers, insolved);
    c = c(insolved);
    lbounds = lbounds(insolved);
    ubounds = ubounds(insolved);
    data_changed = 1;
end;

%----- shift nonzero lower bounds -----
Lbounds_non0 = any(lbounds ~= 0);
if (Lbounds_non0)
    b = b - A*lbounds;
    data_changed = 1;
end

%----- find upper bounds -----
nub = 0; iubounds = ubounds < BIG - 1;
Ubounds_exist = full(any(iubounds));
if (Ubounds_exist)
    ubounds = sparse(iubounds.*(ubounds-lbounds));
    nub = nnz(ubounds);
end

[m, n] = size(A); NNZA = nnz(A);
if DISPLAY_LIPSOL >= 1,
    fprintf(' (m=%i, n=%i)\n',m,n);
end
