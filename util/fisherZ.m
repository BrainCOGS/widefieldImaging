function z = fisherZ(r)

%  z = fisherZ(r)
% fisher' z tranformation

r(r==1)  =  .999;
r(r==-1) = -.999;
z        = .5 .* log((1+r) ./ (1-r));