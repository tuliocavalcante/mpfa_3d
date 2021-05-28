function [ tinf, fint ] = insideintnodes( polyt, polyf )
%
tol = 4;
centpolyf = mean(polyf,1);
centpolyt = mean(polyt,1);

if size(polyt,1)>size(polyf,1)
    v1 = polyt(2,:) - polyt(1,:); v2 = polyt(3,:) - polyt(1,:);
else
    v1 = polyf(2,:) - polyf(1,:); v2 = polyf(3,:) - polyf(1,:);
end
n = cross(v1,v2)/norm(cross(v1,v2));

if size(polyt,1)>=3
    c = [1:size(polyt,1) 1];
    for i=1:size(polyt,1)
        Vt(i,:) = polyt(c(i+1),:) - polyt(c(i),:); Vt(i,:) = Vt(i,:)/norm(Vt(i,:)); 
        Nt(i,:) = axisrotation3D(n,pi/2,Vt(i,:))*sign(dot(polyt(i,:)-centpolyt,axisrotation3D(n,pi/2,Vt(i,:))));
        Nt(i,:) = Nt(i,:)/norm(Nt(i,:)); 
        for j=1:size(polyf,1)
            verifint(j,i) = round((polyf(j,:) - polyt(i,:))*Nt(i,:)',tol)<=0;
        end
    end
    posfint = sum(verifint,2)==size(polyt,1);
    fint = polyf(posfint==1,:);
else
    fint = 0;
end

if size(polyf,1)>=3
    c = [1:size(polyf,1) 1];
    for i=1:size(polyf,1)
        Vf(i,:) = polyf(c(i+1),:) - polyf(c(i),:); Vf(i,:) = Vf(i,:)/norm(Vf(i,:)); 
        Nf(i,:) = axisrotation3D(n,pi/2,Vf(i,:))*sign(dot(polyf(i,:)-centpolyf,axisrotation3D(n,pi/2,Vf(i,:))));
        Nf(i,:) = Nf(i,:)/norm(Nf(i,:));
        for j=1:size(polyt,1)
            veritinf(j,i) = round((polyt(j,:) - polyf(i,:))*Nf(i,:)',tol)<=0;
        end
    end
    postinf = sum(veritinf,2)==size(polyf,1);
    tinf = polyt(postinf==1,:);
else
    tinf = 0;
end

end

%-------------------------------------------------------------------------%
function [ vrot ] = axisrotation3D( n, theta, v )
%

nx = [0 -n(3) n(2);
      n(3) 0 -n(1);
      -n(2) n(1) 0];
  
nnT = [n(1)^2 n(1)*n(2) n(1)*n(3);
       n(1)*n(2) n(2)^2 n(2)*n(3);
       n(1)*n(3) n(2)*n(3) n(3)^2];

R = cos(theta)*eye(3) + sin(theta)*nx + (1 - cos(theta))*nnT;

vrot = (R*v')';

end

