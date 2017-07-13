function [x] = afun(x,Ind,params)

if params.logical==1
    x = vec(x);
    x(~Ind)=0;
else
    x = vec(x);
    x(Ind)=0;
end
end