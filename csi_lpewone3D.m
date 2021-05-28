function [ CSI, csik ] = csi_lpewone3D( t1,t2,t3,t4,t5,t6,Kk )
% CSI calculation for LPEW1 of the MPFA-3D --------------------------------
%  t1 = [T2;T1;k;Q]; t2 = [T3;T2;k;Q]; t3 = [T4;T3;k;Q]; t4 = [T5;T4;k;Q];
%  t5 = [T6;T5;k;Q]; t6 = [T1;T6;k;Q];

T(:,:,1)=t1;
T(:,:,2)=t2;
T(:,:,3)=t3;
T(:,:,4)=t4;
T(:,:,5)=t5;
T(:,:,6)=t6;

for i=1:6
    N=zeros(3,6);
    I=T(1,:,i)'; J=T(2,:,i)'; K=T(3,:,i)'; Q=T(4,:,i)';
    % Oposto a Q:
    JI=I-J; JK=K-J; JQ = Q-J; N(:,4) = 0.5*cross(JI,JK); r = dot(JQ,N(:,4));
    if r > 0, N(:,4)=-N(:,4); end
    % Oposto a I:
    KQ=Q-K; KJ=J-K; KI = I-K; N(:,1) = 0.5*cross(KQ,KJ); r = dot(KI,N(:,1));
    if r > 0, N(:,1)=-N(:,1); end
    % Oposto a J:
    KQ=Q-K; KJ=J-K; KI = I-K; N(:,2) = 0.5*cross(KQ,KI); r = dot(KJ,N(:,2));
    if r > 0, N(:,2)=-N(:,2); end
    % Oposto a K:
    JI=I-J; JK=K-J; JQ = Q-J; N(:,3) = 0.5*cross(JI,JQ); r = dot(JK,N(:,3));
    if r > 0, N(:,3)=-N(:,3); end
    v=abs(dot(JQ,cross(JI,JK))/6);
    for j=1:3
        csi(i,j)=(N(:,4)'*Kk*N(:,j))/(3*v);
    end
end

CSI = [csi(1,2)+csi(6,1) csi(2,2)+csi(1,1) csi(3,2)+csi(2,1) csi(4,2)+csi(3,1) csi(5,2)+csi(4,1) csi(6,2)+csi(5,1)];

csik = csi(1,3)+csi(2,3)+csi(3,3)+csi(4,3)+csi(5,3)+csi(6,3);

end

