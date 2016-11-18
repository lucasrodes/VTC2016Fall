function fzv = integrateOverSecondPair(fun,x,y,zmin,zmax,vmin,vmax)

% DBLQUAD evaluates FUN for a vector x and scalar y. The backwards loop is
% Matt Fig's dynamic preallocation.
for ii=numel(x):-1:1
  g = @(z,v) fun(x(ii),y,z,v);
  fzv(ii) = dblquad(g,zmin,zmax,vmin,vmax);
end
