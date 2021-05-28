function pesos3D( Keq, GV )
%
global options vertex face sist

w = zeros(size(vertex.elsurvertex)); a = zeros(size(vertex.coord,1),1);

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
for i=1:size(vertex.elsurvertpointer,2)
    
    if vertex.flag(i)>200
    
        if i==size(vertex.elsurvertpointer,2)
            posvelem = vertex.elsurvertpointer(i):size(vertex.elsurvertex,2);
        else
            posvelem = vertex.elsurvertpointer(i):vertex.elsurvertpointer(i+1)-1;
        end 
        
        velem = vertex.elsurvertex(posvelem);
        
        [ W, A ] = calcWA3D( velem, i );

        if sum(W)~=0
            a(i) = sum(A)/sum(W); W = W/sum(W); w(posvelem) = W;
        end

        clear velem W A

        [ sist.Mcdt, sist.Bcdt ] = contribuvertex3D( sist.Mcdt, sist.Bcdt, Keq, GV, w, a, i );
    
    else
        
        reg4 = vertex.flag(i)-100; 
        pV = options.solanalit{reg4}(vertex.coord(i,1),vertex.coord(i,2),vertex.coord(i,3));
        for j=1:3
            faces = find(face.inner.vertices(:,j)==i);
            if isempty(faces)==0              
                leftofaces = face.inner.montelem(faces);
                rigtofaces = face.inner.juselem(faces);
                for t=1:size(leftofaces,1)
                    sist.Bcdt(leftofaces(t)) = sist.Bcdt(leftofaces(t)) + Keq(faces(t))*GV(faces(t),j)*pV;
                    sist.Bcdt(rigtofaces(t)) = sist.Bcdt(rigtofaces(t)) - Keq(faces(t))*GV(faces(t),j)*pV;
                end
            end
        end
        
    end
    
end

vertex.weights = w; vertex.flowcontribs = a;

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [ W, A ] = calcWA3D( velem, i )
%
global options vertex element

    Q = vertex.coord(i,:); W = zeros(size(velem)); A = zeros(size(velem));

    if strcmp(options.tipopeso,'IDW')==1
        for j=1:size(velem,2)
            W(j) = 1/norm(element.centroid(velem(j),:)-Q);
        end

    elseif strcmp(options.tipopeso,'LSW')==1
        xk = element.centroid(velem,1)-Q(1); xm = mean(xk)*ones(size(velem,2),1);
        yk = element.centroid(velem,2)-Q(2); ym = mean(yk)*ones(size(velem,2),1);
        zk = element.centroid(velem,3)-Q(3); zm = mean(zk)*ones(size(velem,2),1);
        sigmax2 = (1/size(velem,2))*sum((xk-xm).^2);
        sigmay2 = (1/size(velem,2))*sum((yk-ym).^2);
        sigmaz2 = (1/size(velem,2))*sum((zk-zm).^2);
        W = (1/size(velem,2))*(1+(((xm-xk).*xm)./sigmax2)+...
                                 (((ym-yk).*ym)./sigmay2)+...
                                 (((zm-zk).*zm)./sigmaz2));

    elseif strcmp(options.tipopeso,'YG2019')==1
        for j=1:size(velem,2)
            o1=0;o2=0;o3=0; 
            Ko1 = [1 0 0;0 1 0;0 0 1]; Ko2 = [1 0 0;0 1 0;0 0 1]; Ko3 = [1 0 0;0 1 0;0 0 1];
            k = element.centroid(velem(j),:); reg = element.region(velem(j));
            Kk = options.tensor{reg}(k(1),k(2),k(3));
            faceQ = element.vertices(velem(j),find(element.vertices(velem(j),1:4)~=i));
            I = vertex.coord(faceQ(1),:); J = vertex.coord(faceQ(2),:); K = vertex.coord(faceQ(3),:);
            JI = I-J; JK = K-J; JQ = Q-J; N = 0.5*cross(JI,JK); r = dot(JQ,N);
            if r > 0, J1 = K; f1 = faceQ(3); K = J; faceQ(3) = faceQ(2); faceQ(2) = f1; J = J1; end
            for t=1:size(velem,2)
                if (velem(t)~=velem(j))&&(isempty(find(element.vertices(velem(t),1:4)==faceQ(1)))==0)...
                   &&(isempty(find(element.vertices(velem(t),1:4)==faceQ(3)))==0)
                   o1 = velem(t); reg = element.region(o1); cent1 = element.centroid(o1,:);
                   nk1 = -cross(I-Q,K-Q); nk1 = nk1/norm(nk1); no1 = -nk1; h1 = abs(dot(cent1-I,nk1));
                   Ko1 = options.tensor{reg}(element.centroid(o1,1),element.centroid(o1,2),element.centroid(o1,3));
                   Lk1 = (nk1*Kk*nk1')/h1; Lo1 = (no1*Ko1*no1')/h1;
                elseif (velem(t)~=velem(j))&&(isempty(find(element.vertices(velem(t),1:4)==faceQ(2)))==0)...
                   &&(isempty(find(element.vertices(velem(t),1:4)==faceQ(3)))==0)
                   o2 = velem(t); reg = element.region(o2); cent2 = element.centroid(o2,:);
                   nk2 = -cross(K-Q,J-Q); nk2 = nk2/norm(nk2); no2 = -nk2; h2 = abs(dot(cent2-J,nk2));
                   Ko2 = options.tensor{reg}(element.centroid(o2,1),element.centroid(o2,2),element.centroid(o2,3));
                   Lk2 = (nk2*Kk*nk2')/h2; Lo2 = (no2*Ko2*no2')/h2;
                elseif (velem(t)~=velem(j))&&(isempty(find(element.vertices(velem(t),1:4)==faceQ(1)))==0)...
                   &&(isempty(find(element.vertices(velem(t),1:4)==faceQ(2)))==0)
                   o3 = velem(t); reg = element.region(o3); cent3 = element.centroid(o3,:);
                   nk3 = -cross(J-Q,I-Q); nk3 = nk3/norm(nk3); no3 = -nk3; h3 = abs(dot(cent3-J,nk3)); 
                   Ko3 = options.tensor{reg}(element.centroid(o3,1),element.centroid(o3,2),element.centroid(o3,3));
                   Lk3 = (nk3*Kk*nk3')/h3; Lo3 = (no3*Ko3*no3')/h3;
                end
            end
            
            mi = [Lk1/(Lk1+Lo1) Lk2/(Lk2+Lo2) Lk3/(Lk3+Lo3)];
            
            T1 = ( Lk1*k' + Lo1*cent1' + ( Kk' - Ko1' )*nk1' )/( Lk1 + Lo1 ); 
            T2 = ( Lk2*k' + Lo2*cent2' + ( Kk' - Ko2' )*nk2' )/( Lk2 + Lo2 ); 
            T3 = ( Lk3*k' + Lo3*cent3' + ( Kk' - Ko3' )*nk3' )/( Lk3 + Lo3 );             
            
            denphi = ones(4,4); denphi(1:3,1:4) = [k' T1 T2 T3];
            
            for ty=1:4
                numphi = denphi;
                numphi(1:3,ty) = Q';
                phi(ty) = det(numphi)/det(denphi);
            end

            C = (1/size(velem,2))*[phi(1) + mi(1)*phi(2) + mi(2)*phi(3) + mi(3)*phi(4) (1 - mi(1))*phi(2) (1 - mi(2))*phi(3) (1 - mi(3))*phi(4)];
            
            W(j) = W(j) + C(1);
            if o1~=0, tp=find(velem==o1); W(tp) = W(tp) + C(2); end
            if o2~=0, tp=find(velem==o2); W(tp) = W(tp) + C(3); end
            if o3~=0, tp=find(velem==o3); W(tp) = W(tp) + C(4); end
            
        end
            
    else
        for j=1:size(velem,2)
            o1=0;o2=0;o3=0; 
            Ko1 = [1 0 0;0 1 0;0 0 1]; Ko2 = [1 0 0;0 1 0;0 0 1]; Ko3 = [1 0 0;0 1 0;0 0 1];
            k = element.centroid(velem(j),:); reg = element.region(velem(j));
            Kk = options.tensor{reg}(k(1),k(2),k(3));
            faceQ = element.vertices(velem(j),find(element.vertices(velem(j),1:4)~=i));
            I = vertex.coord(faceQ(1),:); J = vertex.coord(faceQ(2),:); K = vertex.coord(faceQ(3),:);
            JI = I-J; JK = K-J; JQ = Q-J; N = 0.5*cross(JI,JK); r = dot(JQ,N);
            if r > 0, J1 = K; f1 = faceQ(3); K = J; faceQ(3) = faceQ(2); faceQ(2) = f1; J = J1; end
            [ f1,f2,f3 ] = identelemface( velem(j),i,faceQ(1),faceQ(2),faceQ(3),element.faces );
            for t=1:size(velem,2)
                if (velem(t)~=velem(j))&&(isempty(find(element.vertices(velem(t),1:4)==faceQ(1)))==0)...
                   &&(isempty(find(element.vertices(velem(t),1:4)==faceQ(3)))==0)
                   o1 = velem(t); reg = element.region(o1);
                   Ko1 = options.tensor{reg}(element.centroid(o1,1),element.centroid(o1,2),element.centroid(o1,3));
                elseif (velem(t)~=velem(j))&&(isempty(find(element.vertices(velem(t),1:4)==faceQ(2)))==0)...
                   &&(isempty(find(element.vertices(velem(t),1:4)==faceQ(3)))==0)
                   o2 = velem(t); reg = element.region(o2);
                   Ko2 = options.tensor{reg}(element.centroid(o2,1),element.centroid(o2,2),element.centroid(o2,3));
                elseif (velem(t)~=velem(j))&&(isempty(find(element.vertices(velem(t),1:4)==faceQ(1)))==0)...
                   &&(isempty(find(element.vertices(velem(t),1:4)==faceQ(2)))==0)
                   o3 = velem(t); reg = element.region(o3);
                   Ko3 = options.tensor{reg}(element.centroid(o3,1),element.centroid(o3,2),element.centroid(o3,3));
                end
            end

            [ CSI,csik,T1,T2,T3,T4,T5,T6 ] = determinecsi3D( I,J,K,Q,Kk,k );

            t1 = [T1;T2;k;Q]; t2 = [T2;T3;k;Q]; t3 = [T3;T4;k;Q]; t4 = [T4;T5;k;Q]; t5 = [T5;T6;k;Q]; t6 = [T6;T1;k;Q];

            [ g1,g2,g3 ] = defineg1g2g3( f1,f2,f3 );

            [ NETA,SIGMA,NBF ] = netas_deltas_pesos3D( t1,t2,t3,t4,t5,t6,Kk,Ko1,Ko2,Ko3,g1,g2,g3,o1,o2,o3,element.centroid,k );

            if strcmp(options.calcpeso,'semiexp')==1
                C = CSI*inv(NETA)*SIGMA;
                D = sum(NBF)+CSI*inv(NETA)*NBF;
            elseif strcmp(options.calcpeso,'compexp')==1
                [ Bb ] = Binvneta( NETA );
                [ C, D ] = calcCneta( CSI, SIGMA, Bb, NBF );
            end

            A(j) = A(j) + D;
            W(j) = W(j) + C(1) + csik;
            if o1~=0, tp=find(velem==o1); W(tp)=W(tp)+C(2); end
            if o2~=0, tp=find(velem==o2); W(tp)=W(tp)+C(3); end
            if o3~=0, tp=find(velem==o3); W(tp)=W(tp)+C(4); end

        end
    end
    
end

%-------------------------------------------------------------------------%
function [ f1,f2,f3 ] = identelemface( k,Q,I,J,K,elemface )
% Verifica faces no contorno.
global face

f = [elemface(k,1) elemface(k,2) elemface(k,3) elemface(k,4)];
f1 = 0; f2 = 0; f3 = 0;

for i=1:4
    if f(i)>size(face.inner.vertices,1)
        a = f(i) - size(face.inner.vertices,1);
        v = face.bound.vertices(a,1:3);
        if (v(1)==Q||v(1)==K||v(1)==I)&&(v(2)==Q||v(2)==K||v(2)==I)&&(v(3)==Q||v(3)==K||v(3)==I)
            f1 = f(i);
        elseif (v(1)==Q||v(1)==K||v(1)==J)&&(v(2)==Q||v(2)==K||v(2)==J)&&(v(3)==Q||v(3)==K||v(3)==J)
            f2 = f(i);
        elseif (v(1)==Q||v(1)==J||v(1)==I)&&(v(2)==Q||v(2)==J||v(2)==I)&&(v(3)==Q||v(3)==J||v(3)==I)
            f3 = f(i);
        end
    end
end

end

%-------------------------------------------------------------------------%
function [ g1,g2,g3 ] = defineg1g2g3( f1,f2,f3 )
%
global face options

n = size(face.inner.vertices,1);
g1 = 0; g2 = 0; g3 = 0; 

if f1~=0
    b1 = face.bound.flag(f1-n);
    c1 = find(options.flagcorresp(:,1)==b1);
    if options.flagcorresp(c1,1)>200
        g1 = options.flagcorresp(c1,2);
    end
end

if f2~=0
    b2 = face.bound.flag(f2-n);
    c2 = find(options.flagcorresp(:,1)==b2);
    if options.flagcorresp(c2,1)>200
        g2 = options.flagcorresp(c2,2);
    end
end

if f3~=0
    b3 = face.bound.flag(f3-n);
    c3 = find(options.flagcorresp(:,1)==b3);
    if options.flagcorresp(c3,1)>200
        g3 = options.flagcorresp(c3,2);
    end
end

end

%-------------------------------------------------------------------------%
function [ B ] = Binvneta( NETA )
%
r = [1 2 3 4 5 6 1 2 3 4 5 6];
B = zeros(6,6);
p = NETA(1,2)*NETA(2,3)*NETA(3,4)*NETA(4,5)*NETA(5,6)*NETA(6,1);
P = NETA(1,1)*NETA(2,2)*NETA(3,3)*NETA(4,4)*NETA(5,5)*NETA(6,6);
p1=1;
for i=1:6
    if NETA(r(i),r(i+1))~=0
        p1=p1*NETA(r(i),r(i+1));
    end
end
C = p1/(p-P);

for i=1:6
    for j=1:6
        if i>j
            N = i-j-1;
        elseif i==j
            N = 5;
        else
            N = 5-j+i;
        end
        num=1;den=1;
        for n=j+1:j+N
            num = num*(-NETA(r(n),r(n)));
        end
        for n=j:j+N
            if NETA(r(n),r(n+1))~=0
                den = den*NETA(r(n),r(n+1));
            end
        end
        B(i,j)=C*(num/den);
    end
end

end

%-------------------------------------------------------------------------%
function [ C,D ] = calcCneta( CSI,SIGMA,B,NBF )
%
D=0;
for j=1:6
    S=0;
    for i=1:6
        S=S+CSI(i)*B(i,j);
    end
    D=D+SIGMA(j,1)*S;
end
C(1)=D;

D=0;E=0;
for j=1:2
    S=0;
    for i=1:6
        S=S+CSI(i)*B(i,j);
    end
    D=D+SIGMA(j,2)*S;
    E=E+NBF(j)*(1+S);
end
C(2)=D;
F(1)=E;

D=0;E=0;
for j=3:4
    S=0;
    for i=1:6
        S=S+CSI(i)*B(i,j);
    end
    D=D+SIGMA(j,3)*S;
    E=E+NBF(j)*(1+S);
end
C(3)=D;
F(2)=E;

D=0;E=0;
for j=5:6
    S=0;
    for i=1:6
        S=S+CSI(i)*B(i,j);
    end
    D=D+SIGMA(j,4)*S;
    E=E+NBF(j)*(1+S);
end
C(4)=D;
F(3)=E;

D=sum(F); 

end

%-------------------------------------------------------------------------%
function [CSI,csik,T1,T2,T3,T4,T5,T6] = determinecsi3D(I,J,K,Q,Kk,k)
%
global options

if strcmp(options.tipopeso,'LPEW1')==1
    
    T1 = 0.5*(Q+I); T2 = (1/3)*(Q+K+I); T3 = 0.5*(Q+K); T4 = (1/3)*(Q+K+J);
    T5 = 0.5*(Q+J); T6 = (1/3)*(Q+J+I);

    t1 = [T2;T1;k;Q]; t2 = [T3;T2;k;Q]; t3 = [T4;T3;k;Q]; t4 = [T5;T4;k;Q];
    t5 = [T6;T5;k;Q]; t6 = [T1;T6;k;Q];

    [ CSI, csik ] = csi_lpewone3D( t1,t2,t3,t4,t5,t6,Kk );

elseif strcmp(options.tipopeso,'LPEW2')==1

    T1 = 0.5*(Q+I); T2 = (1/3)*(Q+K+I); T3 = 0.5*(Q+K); T4 = (1/3)*(Q+K+J); T5 = 0.5*(Q+J); T6 = (1/3)*(Q+J+I);
    t1 = [T1;T6;T2;Q]; t2 = [T2;T4;T3;Q]; t3 = [T4;T6;T5;Q]; t4 = [T2;T6;T4;Q];

    [ CSI ] = csi_pesos3D( t1,t2,t3,t4,Kk );
    
    csik = 0;
    
elseif strcmp(options.tipopeso,'LPEW3')==1
    
    wQ = 0; wK = 1;  wI = 1;  wJ = 1;
    T1 = (1/(wQ+wI))*(wQ*Q+wI*I); T2 = (1/(wQ+wK+wI))*(wQ*Q+wI*I+wK*K); 
    T3 = (1/(wQ+wK))*(wQ*Q+wK*K); T4 = (1/(wQ+wK+wJ))*(wQ*Q+wJ*J+wK*K);
    T5 = (1/(wQ+wJ))*(wQ*Q+wJ*J); T6 = (1/(wQ+wI+wJ))*(wQ*Q+wJ*J+wI*I);
    t1 = [T1;T6;T2;Q]; t2 = [T2;T4;T3;Q]; t3 = [T4;T6;T5;Q]; t4 = [T2;T6;T4;Q];

    [ CSI ] = csi_pesos3D( t1,t2,t3,t4,Kk );
    
    csik = 0;

end

end


