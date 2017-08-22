function [fNew,step,rNew,iter,err, localProdA] = spgLine(f,d,gtd,x,fMax,funForward, funPenalty, params, b,feasSrchIt,linear)
% Nonmonotone linesearch.

localProdA = 0;
EXIT_CONVERGED  = 0;
EXIT_LINEITERATIONS = 1;
maxIts = feasSrchIt;
step   = 1;
iter   = 0;
gamma  = 1e-4;
gtd    = -abs(undist(gtd)); % 03 Aug 07: If gtd is complex,
% then should be looking at -abs(gtd).

%Ad = Aprod(d,1);

if(linear)
Ad = funForward(d, [], params); 
localProdA = localProdA + 1;
r = x; %CAREFUL HERE: we are passing in rOld if linear. 
end

while 1

% Evaluate trial point and function value.
if(linear)
rNew = r - step*Ad;
else
rNew = b - funForward(x + step*d, [], params);
localProdA = localProdA + 1;
end


fNew = funPenalty(rNew, params);
%fNew = funComposite(x + step*d, b, funForward, funPenalty, params);
%    rNew = r - step*Ad;
%    fNew = norm(rNew)^2 / 2;

% Check exit conditions.
if fNew < fMax + gamma*step*gtd  % Sufficient descent condition.
err = EXIT_CONVERGED;
break
elseif  iter >= maxIts           % Too many linesearch iterations.
err = EXIT_LINEITERATIONS;
break
end

% New linesearch iteration.
iter = iter + 1;

% Safeguarded quadratic interpolation.
if step <= 0.1
step  = step / 2;
else
tmp = (-gtd*step^2) / (2*(fNew-f-step*gtd));
if tmp < 0.1 || tmp > 0.9*step || isnan(tmp)
tmp = step / 2;
end
step = tmp;
end

end % while 1

end % function spgLine
