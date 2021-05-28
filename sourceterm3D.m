function [ s ] = sourceterm3D
%
global options element fracture

num_volumes = size(element.vertices,1) + size(fracture.coord,1);
s = zeros(num_volumes,1);
 
syms x y z
reg = unique(element.region);
for i=1:size(reg,1)
    g = [diff(options.solanalit{reg(i)},x);diff(options.solanalit{reg(i)},y);diff(options.solanalit{reg(i)},z)];
    flux = options.tensor{reg(i)}*g;
    q = divergence([flux(1), flux(2), flux(3)], [x, y, z]);
    Q{i} = matlabFunction(q);    
end 

for i=1:size(element.vertices,1)
    
    if nargin(Q{element.region(i)})==0
        S = Q{element.region(i)}();
    elseif nargin(Q{element.region(i)})==3
        S = Q{element.region(i)}(element.centroid(i,1),element.centroid(i,2),element.centroid(i,3));
    elseif nargin(Q{element.region(i)})==2
        S = Q{elem(i,5)}(element.centroid(i,1),element.centroid(i,2));
    elseif options.caso==4
        S = Q{element.region(i)}(element.centroid(i,3));
    end

    s(i,1) = -S*element.volume(i);

end

end


