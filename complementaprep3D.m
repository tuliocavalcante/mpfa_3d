function [ Keq, GV ] = complementaprep3D
%
global face

num_faces = size(face.inner.vertices,1)+size(face.bound.vertices,1);
Keq = zeros(num_faces,1); GV = zeros(num_faces,3);

for ti = 1:size(face.inner.vertices,1)
    [ Keq, GV ] = calc_M_GV3D( ti, Keq, GV );
end

for tb = 1:size(face.bound.vertices,1)
    calc_M_B_bound3D( tb );
end

end