function [ M, B ] = contribuvertex3D( M, B, Keq, GV, w, a, i )
%
global face vertex

inicesurn = vertex.elsurvertpointer(i);
if i==size(vertex.elsurvertpointer,2)
    fimesurn = size(vertex.elsurvertex,2);
else
    fimesurn = vertex.elsurvertpointer(i+1)-1;
end
col = vertex.elsurvertex(inicesurn:fimesurn);
wei = w(inicesurn:fimesurn);
for j=1:3
    faces = find(face.inner.vertices(:,j)==i);
    if isempty(faces)==0
        elemfaces = [face.inner.montelem(faces); face.inner.juselem(faces)];
        contribM = [-(Keq(faces).*GV(faces,j))*wei; (Keq(faces).*GV(faces,j))*wei];
        contribB = [Keq(faces).*GV(faces,j); -Keq(faces).*GV(faces,j)];
        [elemfaces,idx] = sort(elemfaces);
        contribM = contribM(idx,:);
        contribB = contribB(idx,:);
        facevec = unique(elemfaces);
        C = zeros(size(facevec,1),size(wei,2));
        I = zeros(size(facevec,1),1);
        t = 1;
        for k = facevec'
            lin = find(elemfaces==k);
            C(t,:) = sum(contribM(lin,:),1);
            I(t,1) = sum(contribB(lin,:));
            t = t +1;
        end
        M(facevec,col) = M(facevec,col) + C;
        if vertex.flag(i)<300
            B(facevec) = B(facevec) + I*a(i);
        end
    end
end

end

