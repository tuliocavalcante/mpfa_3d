function calc_M_B_bound3D(i)
%
global vertex element face sist options

for j=1:3
   vertex.flag(face.bound.vertices(i,j)) = min(face.bound.flag(i),vertex.flag(face.bound.vertices(i,j)));
end    

if face.bound.flag(i)<200
    
    I = vertex.coord(face.bound.vertices(i,1),:)';
    J = vertex.coord(face.bound.vertices(i,2),:)';
    K = vertex.coord(face.bound.vertices(i,3),:)';
    nL = face.bound.montelem(i); L = element.centroid(nL,:)';
    
    reg4 = element.region(nL); KkL = options.tensor{reg4}(L(1),L(2),L(3));

    JI = I - J; JK = K - J; LJ = J - L;
    N = face.bound.normal(i,:)'; TJI = cross(N,JI); TJK = cross(N,JK);

    reg4 = face.bound.flag(i) - 100;
    pI = options.solanalit{reg4}(I(1),I(2),I(3));
    pJ = options.solanalit{reg4}(J(1),J(2),J(3));
    pK = options.solanalit{reg4}(K(1),K(2),K(3));   
    
    hL = abs(dot(LJ,N)/norm(N)); KnL = (N'*KkL*N)/(norm(N)^2);
    KTJIL = (N'*KkL*TJI)/(norm(N)^2); KTJKL = (N'*KkL*TJK)/(norm(N)^2);

    sist.Btpfa(nL) = sist.Btpfa(nL) - ((1/(2*hL*norm(N)))*((dot(-TJK,LJ)*KnL+hL*norm(N)*KTJKL)*...
                                      (pI-pJ)-2*(norm(N)^2)*KnL*pJ+(dot(-TJI,LJ)*KnL+hL*norm(N)*KTJIL)*(pJ-pK)));

    sist.Mtpfa(nL,nL) = sist.Mtpfa(nL,nL) + (1/(hL*norm(N)))*((norm(N)^2)*KnL);

elseif face.bound.flag(i)>200
    
    nL = face.bound.montelem(i); bcflag = options.flagcorresp;
    
    c1 = find(bcflag(:,1)==face.bound.flag(i));
    if bcflag(c1,1)>200, g = bcflag(c1,2); end
    
    N = face.bound.normal(i,:)';
    
    sist.Btpfa(nL) = sist.Btpfa(nL) - norm(N)*g;
    
end

end

