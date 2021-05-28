function [ q, influx, bflux, G, NG ] = flowbalance3D( p, Keq, GV )
%
global options vertex element face

syms x y z
reg = unique(element.region);
for i=1:size(reg,1)
    g = [diff(options.solanalit{reg(i)},x);diff(options.solanalit{reg(i)},y);diff(options.solanalit{reg(i)},z)];
    grad{i} = matlabFunction(g);
end 

q = zeros(size(p,1),1); NG = zeros(size(p,1),3); regel = zeros(size(p,1),1); G = zeros(size(p,1),3);
pno = zeros(size(vertex.elsurvertpointer,2),1); regno = zeros(size(vertex.elsurvertpointer,2),1);

for i=1:size(face.inner.vertices,1)
    
    for j=1:4
        if regno(element.vertices(face.inner.montelem(i),j))==0
            [ pno ] = recuperapressnos3D( p, pno, element.vertices(face.inner.montelem(i),j) );
            regno(element.vertices(face.inner.montelem(i),j)) = 1;
        end
        if regno(element.vertices(face.inner.juselem(i),j))==0
            [ pno ] = recuperapressnos3D( p, pno, element.vertices(face.inner.juselem(i),j) );
            regno(element.vertices(face.inner.juselem(i),j)) = 1;
        end
    end
    
    [ G, regel ] = calcG3D( G, regel, face.inner.montelem(i), grad );
    [ G, regel ] = calcG3D( G, regel, face.inner.juselem(i), grad );
        
    nosL = element.vertices(face.inner.montelem(i),:); nosR = element.vertices(face.inner.juselem(i),:);
    IJK = face.inner.vertices(i,1:3); noQL = setdiff(nosL,IJK); noQR = setdiff(nosR,IJK);
    
    NG(face.inner.montelem(i),:) = NG(face.inner.montelem(i),:)-(face.inner.normal(i,:)*pno(noQL))*(1/(3*element.volume(face.inner.montelem(i))));
    NG(face.inner.juselem(i),:) = NG(face.inner.juselem(i),:)+(face.inner.normal(i,:)*pno(noQR))*(1/(3*element.volume(face.inner.juselem(i))));
    
    influx(i,1) = -Keq(i)*((p(face.inner.juselem(i))-p(face.inner.montelem(i)))...
                  +GV(i,1)*pno(face.inner.vertices(i,1))+GV(i,2)*pno(face.inner.vertices(i,2))...
                  +GV(i,3)*pno(face.inner.vertices(i,3)));
    
    q(face.inner.montelem(i)) = q(face.inner.montelem(i)) + influx(i,1);
    q(face.inner.juselem(i)) = q(face.inner.juselem(i)) - influx(i,1);
    
end

for i=1:size(face.bound.vertices,1)

if face.bound.flag(i)<200
    
    for j=1:4
        if regno(element.vertices(face.bound.montelem(i),j))==0
            [ pno ] = recuperapressnos3D( p, pno, element.vertices(face.bound.montelem(i),j) );
            regno(element.vertices(face.bound.montelem(i),j)) = 1;
        end
    end
    
    I = vertex.coord(face.bound.vertices(i,1),:)'; 
    J = vertex.coord(face.bound.vertices(i,2),:)'; 
    K = vertex.coord(face.bound.vertices(i,3),:)';

    L = element.centroid(face.bound.montelem(i),:)';
    
    reg4 = element.region(face.bound.montelem(i));
    KkL = options.tensor{reg4}(L(1),L(2),L(3));

    JI = I - J; JK = K - J; LJ = J - L;     
    N = face.bound.normal(i,:)'; TJI=cross(N,JI); TJK=cross(N,JK);
    
    nosL = element.vertices(face.bound.montelem(i),:);
    IJK = face.bound.vertices(i,:); noQL = setdiff(nosL,IJK);
    
    NG(face.bound.montelem(i),:) = NG(face.bound.montelem(i),:)-(face.bound.normal(i,:)*pno(noQL))*(1/(3*element.volume(face.bound.montelem(i))));

    hL = abs(dot(LJ,N)/norm(N));

    KnL = (N'*KkL*N)/(norm(N)^2);

    KTJIL = (N'*KkL*TJI)/(norm(N)^2);
    KTJKL = (N'*KkL*TJK)/(norm(N)^2);

    bflux(i,1) = (1/(2*hL*norm(N)))*((dot(-TJK,LJ)*KnL+hL*norm(N)*KTJKL)*...
                 (pno(face.bound.vertices(i,1))-pno(face.bound.vertices(i,2)))+2*(norm(N)^2)*KnL*...
                 (p(face.bound.montelem(i))-pno(face.bound.vertices(i,2)))+(dot(-TJI,LJ)*KnL+hL*...
                 norm(N)*KTJIL)*(pno(face.bound.vertices(i,2))-pno(face.bound.vertices(i,3))));
                               
    q(face.bound.montelem(i)) = q(face.bound.montelem(i)) + bflux(i,1);

else
    
    bflux(i,1) = 0;
    
end

end

q = q - element.sourceterm;

end

%-------------------------------------------------------------------------%
function [ pno ] = recuperapressnos3D( p, pno, no )
%
global vertex options

if vertex.flag(no)<200
    reg4 = vertex.flag(no) - 100;
    I = vertex.coord(no,:)';
    pno(no) = options.solanalit{reg4}(I(1),I(2),I(3));
else
    inicesurn = vertex.elsurvertpointer(no);
    if no==size(vertex.elsurvertpointer,2)
        fimesurn = size(vertex.elsurvertex,2);
    else
        fimesurn = vertex.elsurvertpointer(no+1)-1;
    end
    pre = p(vertex.elsurvertex(inicesurn:fimesurn));
    wei = vertex.weights(inicesurn:fimesurn);
    pno(no) = wei*pre;
end

end

%-------------------------------------------------------------------------%
function [ G, regel ] = calcG3D( G, regel, adjelem, grad )
%
global element

reg4 = element.region(adjelem);

if regel(adjelem)==0
    if nargin(grad{reg4})==0
        G(adjelem,:) = grad{reg4}()';
    elseif nargin(grad{reg4})==3
        G(adjelem,:) = grad{reg4}(element.centroid(adjelem,1),element.centroid(adjelem,2),element.centroid(adjelem,3))';
    elseif nargin(grad{reg4})==2
        G(adjelem,:) = grad{reg4}(element.centroid(adjelem,1),element.centroid(adjelem,2))';
    elseif caso==4
        G(adjelem,:) = grad{reg4}(element.centroid(adjelem,3))';
    elseif caso==2
        G(adjelem,:) = grad{reg4}(element.centroid(adjelem,2))';
    end
    regel(adjelem) = 1;
end

end

