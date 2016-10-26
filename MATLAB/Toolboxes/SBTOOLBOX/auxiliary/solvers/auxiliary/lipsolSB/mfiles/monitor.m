function monitor(iter,rb,rc,ru,dgap)
% MONITOR     - Graphic monitor for iteration progress.

% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global rb_prev rc_prev ru_prev dg_prev

if iter == 0
   clf;
   drawnow;
   ymax = min(16,ceil(1+log10(max([rb,rc,ru,dgap]))));
   axis([0 50 -16 ymax]); axis(axis);  hold on;
   title('LIPSOL Monitor'); 
   xlabel('Iteration'); 
   ylabel('Log Residuals');
   text(25,ymax-1.5,'solid line: duality gap');
   text(25,ymax-3.0,'dotted line: p-infeasibility');
   text(25,ymax-4.5,'dashdot line: d-infeasibility');
   text(25,ymax-6.0,'(dashed line: upper bounds)');
   rb   = max(1.e-32,rb);   rb_prev = rb;
   rc   = max(1.e-32,rc);   rc_prev = rc;
   ru   = max(1.e-32,ru);   ru_prev = ru;
   dgap = max(1.e-32,dgap); dg_prev = dgap;
else
   rb   = max(1.e-32,rb);
   rc   = max(1.e-32,rc);
   ru   = max(1.e-32,ru);
   dgap = max(1.e-32,dgap);
   plot([iter-1 iter], log10([rb_prev rb]),   ':', ...
        [iter-1 iter], log10([rc_prev rc]),   '-.', ...
        [iter-1 iter], log10([ru_prev ru]),   '--', ...
        [iter-1 iter], log10([dg_prev dgap]), '-' );
   drawnow;
   rb_prev = rb;
   rc_prev = rc;
   ru_prev = ru;
   dg_prev = dgap;
end;
