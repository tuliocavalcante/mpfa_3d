function [ face, element, fracture ] = portionvolumes( coordtetra, N1, N2, N3, N4, Nf, ...
                                       intpolynodes, area, tetraj, frati, face, element, ...
                                       fracture, vertex, options )
%

for i=1:4   
    fractetra(i,:) = coordtetra(i,:) - mean(intpolynodes,1);
    veri(1,i) = (Nf(frati,:)/norm(Nf(frati,:)))*fractetra(i,:)'>0; 
end
posmaior = find(veri>0); posmenor = find(veri==0);

cenf = mean(intpolynodes,1); dist = 0;

poly1 = [intpolynodes; coordtetra(posmaior,:)]; 
poly2 = [intpolynodes; coordtetra(posmenor,:)];

DT1 = delaunayTriangulation(poly1(:,1),poly1(:,2),poly1(:,3));
DT2 = delaunayTriangulation(poly2(:,1),poly2(:,2),poly2(:,3));
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
portionvol1 = 0;
for i=1:size(DT1.ConnectivityList,1)
    
    tetra1 = poly1(DT1.ConnectivityList(i,:),:); 
    
    cent1 = mean(tetra1,1); 
    
    noI = tetra1(1,:); noJ = tetra1(2,:);
    noK = tetra1(3,:); noQ = tetra1(4,:);

    IQ = noQ - noI; IJ = noJ - noI; IK = noK - noI;
    
    vol1 = abs(dot(cross(IJ,IK),IQ)/6);
    
    dist = dist + abs(dot((Nf(frati,:)/norm(Nf(frati,:))),cent1-cenf)*vol1);

    portionvol1 = portionvol1 + vol1;
    
end

portionvol2 = 0;
for i=1:size(DT2.ConnectivityList,1)
    
    tetra2 = poly2(DT2.ConnectivityList(i,:),:);
    
    cent2 = mean(tetra2,1);
    
    noI = tetra2(1,:); noJ = tetra2(2,:);
    noK = tetra2(3,:); noQ = tetra2(4,:);

    IQ = noQ - noI; IJ = noJ - noI; IK = noK - noI;
    
    vol2 = abs(dot(cross(IJ,IK),IQ)/6);
    
    dist = dist + abs(dot((Nf(frati,:)/norm(Nf(frati,:))),cent2-cenf)*vol2);

    portionvol2 = portionvol2 + vol2;
    
end

if portionvol2>portionvol1, volmenor = portionvol1; volmaior = portionvol2; end
if portionvol1>portionvol2, volmenor = portionvol2; volmaior = portionvol1; end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
n1 = Nf(frati,:)/norm(Nf(frati,:)); n2 = -Nf(frati,:)/norm(Nf(frati,:));
voltotal = volmaior + volmenor; dist = dist/voltotal;

OF = mean(intpolynodes,1) - mean(coordtetra,1);

if dot(n1,OF)>=0, n = n1; elseif dot(n2,OF)>0, n = n2; end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
if dot(n,N1(tetraj,:))>=0
    factor(1) = area*abs(dot(n,N1(tetraj,:))/norm(N1(tetraj,:)))*(1 - (volmenor/voltotal)); 
else
    factor(1) = area*abs(dot(n,N1(tetraj,:))/norm(N1(tetraj,:)))*(1 - (volmaior/voltotal)); 
end

if dot(n,N2(tetraj,:))>=0
    factor(2) = area*abs(dot(n,N2(tetraj,:))/norm(N2(tetraj,:)))*(1 - (volmenor/voltotal)); 
else
    factor(2) = area*abs(dot(n,N2(tetraj,:))/norm(N2(tetraj,:)))*(1 - (volmaior/voltotal)); 
end

if dot(n,N3(tetraj,:))>=0
    factor(3) = area*abs(dot(n,N3(tetraj,:))/norm(N3(tetraj,:)))*(1 - (volmenor/voltotal)); 
else
    factor(3) = area*abs(dot(n,N3(tetraj,:))/norm(N3(tetraj,:)))*(1 - (volmaior/voltotal)); 
end

if dot(n,N4(tetraj,:))>=0
    factor(4) = area*abs(dot(n,N4(tetraj,:))/norm(N4(tetraj,:)))*(1 - (volmenor/voltotal)); 
else
    factor(4) = area*abs(dot(n,N4(tetraj,:))/norm(N4(tetraj,:)))*(1 - (volmaior/voltotal)); 
end
%-------------------------------------------------------------------------%

numface = element.faces(tetraj,:);
for y=1:4
    if numface(y)<=size(face.inner.vertices,1)
        if face.inner.montelem(numface(y))==tetraj
            fracture.connecwithcell(frati,face.inner.juselem(numface(y))) = factor(y);
        else
            fracture.connecwithcell(frati,face.inner.montelem(numface(y))) = factor(y);
        end          
    else
        boundface = numface(y)-size(face.inner.vertices,1);
        fracture.connecwithbound(frati,boundface) = factor(y);
        onbound = zeros(2,3); tw = 1;
        for w=1:size(intpolynodes,1)
            verif = verifypointinplpane(intpolynodes(w,:),face.bound.vertices(boundface,:),vertex.coord);
            if verif==1
                onbound(tw,:) = intpolynodes(w,:); tw = tw +1;
            end
        end 
        fracture.slicingbound(frati,boundface) = 0;
        for w=2:size(onbound,1)
            fracture.slicingbound(frati,boundface) = fracture.slicingbound(frati,boundface) + norm(onbound(w,:) - onbound(w-1,:)); 
        end
        if size(onbound,1)>=2
            reg = face.bound.flag(boundface) - 100; centb = mean(onbound,1); centf = mean(intpolynodes,1);
            veb = onbound(2,:) - onbound(1,:); vec = centf - onbound(1,:); 
            vec1 = vec/(norm(vec)+eps*(norm(vec)==0)); veb1 = veb/norm(norm(veb)+eps*(norm(veb)==0)); 
            theta = acos(dot(vec1,veb1)); h = norm(vec)*sin(theta); 
            if face.bound.flag(boundface)<200
                fracture.valueonbound(frati,boundface) = options.solanalit{reg}(centb(1),centb(2),centb(3));
            end
            fracture.slicingbound(frati,boundface) = fracture.slicingbound(frati,boundface);
            fracture.disttobound(frati,boundface) = h;
        end
    end    
end

element.areafrac(tetraj,frati) = area;
element.distfrac(tetraj,frati) = dist;

end

%-------------------------------------------------------------------------%
function verif = verifypointinplpane(node,face,coord)

A = ones(4,4);

A(1:3,1:3) = coord(face,:)';
A(1:3,4) = node';

verif = det(A)==0;

end
