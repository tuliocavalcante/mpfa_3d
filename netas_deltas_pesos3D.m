function [ NETA,SIGMA,NBF ] = netas_deltas_pesos3D( t1,t2,t3,t4,t5,t6,...
                            Kk,Ko1,Ko2,Ko3,g1,g2,g3,o1,o2,o3,centelem,k )

% t1 = [T1;T2;k;Q]; t2 = [T2;T3;k;Q]; t3 = [T3;T4;k;Q]; t4 = [T4;T5;k;Q];
% t5 = [T5;T6;k;Q]; t6 = [T6;T1;k;Q];
                                         
NETA = zeros(6,6);
SIGMA = zeros(6,4);
NBF = zeros(6,1);

% eps=1e-10;
T(:,:,1)=t1;
T(:,:,2)=t2;
T(:,:,3)=t3;
T(:,:,4)=t4;
T(:,:,5)=t5;
T(:,:,6)=t6;
% theta=0.3*pi;
% R=[cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0;0 0 1];
if o1~=0, O1=centelem(o1,:)'; end
if o2~=0, O2=centelem(o2,:)'; end
if o3~=0, O3=centelem(o3,:)'; end
% if o1~=0, O1=centelem(o1,:)'; kO1=O1'-k; k=k+(1e-5)*(R*kO1')'; end
% if o2~=0, O2=centelem(o2,:)'; kO2=O2'-k; k=k+(1e-5)*(R*kO2')'; end
% if o3~=0, O3=centelem(o3,:)'; kO3=O3'-k; k=k+(1e-5)*(R*kO3')'; end
% O = [o1 o2 o3];

for j=1:6
    % Oposto a k:
    Q=T(4,:,j)';
    T1=T(1,:,j)';
    T2=T(2,:,j)';
    T1Q=Q-T1;
    T2Q=Q-T2;
    kQ=Q-k';
    N=0.5*cross(-T2Q,-T1Q); 
    Nf(:,j)=N;
    if (j==1 || j==2) && o1~=0
        TO=O1-T1;
        Ko=Ko1;
        QO=O1-Q; 
        ho(j)=abs(dot(TO,N))/norm(N);
    elseif (j==3 || j==4) && o2~=0 
        TO=O2-T1;
        Ko=Ko2;
        QO=O2-Q;
        ho(j)=abs(dot(TO,N))/norm(N);
    elseif (j==5 || j==6) && o3~=0
        TO=O3-T1;
        Ko=Ko3;
        QO=O3-Q;
        ho(j)=abs(dot(TO,N))/norm(N);
    else
        Ko=0;
        QO=[0 0 0]'; 
        ho(j)=Inf;
    end
    hk(j)=abs(dot(kQ,N))/norm(N);
    Knk(j)=(N'*Kk*N)/(norm(N)^2);
    Kno(j)=(N'*Ko*N)/(norm(N)^2);
    sigmak(j)=norm(N)*Knk(j)/hk(j);
    sigmao(j)=norm(N)*Kno(j)/ho(j);
    for i=1:2
        QT(:,i)=T(i,:,j)'-Q;
        tQT(:,i)=cross(N,QT(:,i));
        Ktk(j,i)=(N'*Kk*tQT(:,i))/(norm(N)^2);
        Kto(j,i)=(N'*Ko*tQT(:,i))/(norm(N)^2); 
        neta(j,i)=-((dot(tQT(:,i),kQ))/(2*hk(j)*norm(N)))*Knk(j)-...
        ((dot(tQT(:,i),-QO))/(2*ho(j)*norm(N)))*Kno(j)+0.5*Ktk(j,i)-0.5*Kto(j,i);
    end
end

for i=1:6      
    NETA(i,i)=-neta(i,2);
    if i<6 
        NETA(i,i+1)=neta(i,1);
    else
        NETA(i,1)=neta(i,1);
    end
end

for i=1:6      
   for j=1:4
      SIGMA(i,1)=-sigmak(i);
      if i==1 || i==2
         SIGMA(i,2)=-sigmao(i);
         if o1==0
             NBF(i,1)=norm(Nf(:,i))*g1; 
         end
      end
      if i==3 || i==4
         SIGMA(i,3)=-sigmao(i);
          if o2==0
             NBF(i,1)=norm(Nf(:,i))*g2; 
          end
      end
      if i==5 || i==6
         SIGMA(i,4)=-sigmao(i);
          if o3==0
             NBF(i,1)=norm(Nf(:,i))*g3;  
          end
      end
   end  
end

end

