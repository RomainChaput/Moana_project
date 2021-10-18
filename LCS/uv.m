function out = uv(t, xy, Fu, Fv)
% t : time of the interpolation
% xy: list of coordinates xy = [x1,...,xn,y1,...,yn];
% veloicities: out = [u1,...,un,v1,...,vn];

n = length(xy)/2;

x = xy(1:n);
y = xy(n+1:2*n);
t = repmat(t, [n 1]);

% output velocities
out = zeros(2*n,1);
out(1:n) = Fu(x, y, t);
out(n+1:2*n) = Fv(x, y, t);

end