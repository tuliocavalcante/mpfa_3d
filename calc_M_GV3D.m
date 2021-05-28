function [ Keq, GV ] = calc_M_GV3D( i, Keq, GV )
%
global vertex element face sist options

nL = face.inner.montelem(i); reg4 = element.region(nL);
nR = face.inner.juselem(i); reg5 = element.region(nR);
L = element.centroid(nL,:)'; R = element.centroid(nR,:)';

KkL = options.tensor{reg4}(L(1),L(2),L(3));
KkR = options.tensor{reg5}(R(1),R(2),R(3));

I = vertex.coord(face.inner.vertices(i,1),:)';
J = vertex.coord(face.inner.vertices(i,2),:)';
K = vertex.coord(face.inner.vertices(i,3),:)';

JI = I - J; JK = K - J; LR = R - L; LJ = J - L; RJ = J - R;

N = face.inner.normal(i,:)'; TJI=cross(N,JI); TJK=cross(N,JK); 

hL = abs(dot(LJ,N)/norm(N)); hR = abs(dot(RJ,N)/norm(N));

KnL = (N'*KkL*N)/(norm(N)^2); KnR = (N'*KkR*N)/(norm(N)^2);

KTJIL = (N'*KkL*TJI)/(norm(N)^2); KTJKR = (N'*KkR*TJK)/(norm(N)^2);
KTJIR = (N'*KkR*TJI)/(norm(N)^2); KTJKL = (N'*KkL*TJK)/(norm(N)^2);

DJI = (dot(TJI,LR)/(norm(N)^2))-(1/norm(N))*(hL*(KTJIL/KnL)+hR*(KTJIR/KnR));
DJK = (dot(TJK,LR)/(norm(N)^2))-(1/norm(N))*(hL*(KTJKL/KnL)+hR*(KTJKR/KnR));

Keq(i,1) = ((KnR*KnL)/(KnR*hL+KnL*hR))*norm(N);

% G=0.5*(-DJK*(pJ-pI)+DJI*(pJ-pK)); --------------------------------------%
GI = 0.5*DJK; GJ = 0.5*DJI-0.5*DJK; GK = -0.5*DJI; GV(i,:)=[GI GJ GK];

sist.Mtpfa(nL,nL) = sist.Mtpfa(nL,nL) + Keq(i);
sist.Mtpfa(nL,nR) = sist.Mtpfa(nL,nR) - Keq(i);
sist.Mtpfa(nR,nR) = sist.Mtpfa(nR,nR) + Keq(i);
sist.Mtpfa(nR,nL) = sist.Mtpfa(nR,nL) - Keq(i);

end

