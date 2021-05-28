function [ vertex, element, face, options, sist, fracture, wells ] = preprocessor3D

[ options.caso, options.malha, options.tipopeso, options.calcpeso, ...
  options.limitadora, options.fractarch, options.flagcorresp ] = readopt3D;

filepath = sprintf('Malhas\\%s',options.malha);

[~,nameso,~] = fileparts(options.malha);
folderpath = sprintf('Malhas\\%s',nameso);

if exist(sprintf('%s\\vertex',folderpath),'file') == 2
    
    load(sprintf('%s\\vertex',folderpath));
    load(sprintf('%s\\element',folderpath));
    load(sprintf('%s\\face',folderpath));
    
else
    
    if exist(folderpath,'dir')==0
        mkdir(folderpath);
    end
    
    %"coord" matrix - The coordinates (x,y,z) of each point in a carts: 
    [ coord, nnode, nflag ] = getcoord( filepath );

    %"elem" matrix - The points that constitute each element:
    [ elem, bedgeref, nflag ] = getelem( filepath, nnode, nflag );

    %"bedge" and "inedge" - It get inform. about boundary and internal edges:
    [ bedge, inedge, centelem, elemvolume ] = getinfo( coord, elem, bedgeref );

    %"esurn1" and "esurn2" - Elements surrounding nodes:
    [ esurn1, esurn2 ] = getesurn( elem, coord );

    %"elemface" and "normals" - Faces in elements with normals:
    [ elemface, innormals, bounnormals ] = idelemface3D( elem, inedge, bedge, coord, centelem );

    vertex.coord = coord; vertex.elsurvertex = esurn1; 
    vertex.elsurvertpointer = esurn2;
    vertex.flag = nflag;

    element.vertices = elem(:,1:4); element.centroid = centelem;
    element.volume = elemvolume; element.faces = elemface;
    element.region = elem(:,5);   

    face.bound.vertices = bedge(:,1:3); face.bound.flag = bedge(:,5);
    face.bound.montelem = bedge(:,4); face.bound.normal = bounnormals; 
    face.inner.vertices = inedge(:,1:3); face.inner.montelem = inedge(:,4); 
    face.inner.juselem = inedge(:,5); face.inner.normal = innormals; 

    save(sprintf('%s\\vertex',folderpath),'vertex');
    save(sprintf('%s\\element',folderpath),'element');
    save(sprintf('%s\\face',folderpath),'face');
    
end
%--------------------------------------------------------------------------
[ options.solanalit, options.tensor, wells ] = caso3D( options.caso, element, vertex );

[ fracture ] = gerafrac3D(options);

[ element, fracture, face ] = polygoninpolyhedron ( element, fracture, face, vertex, options );

num_volumes = size(elem,1) + size(fracture.coord,1);
sist.Mtpfa = sparse(num_volumes,num_volumes);
sist.Mcdt = sparse(num_volumes,num_volumes);
sist.Btpfa = sparse(num_volumes,1);
sist.Bcdt = sparse(num_volumes,1);

if isempty(fracture.coord)==0
    sist = assemblyfract( sist, fracture, element, options, face );
end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%FUNCTION DEFINITION
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%FUNCTION "getcoord"
%--------------------------------------------------------------------------
function [coord,nnode,nflag] = getcoord(filepath)
%Open the *.msh file
readmsh = fopen(filepath);

%"nnode" is the number of nodes in the discrete domain
getmshdata = textscan(readmsh,'%u',1,'HeaderLines',4);
%Attribute the data to "nnode"
nnode = getmshdata{1};

%"coord" is a matrix which contain the coordinate of each point of domain
getmshdata = textscan(readmsh,'%*u %f64 %f64 %f64',nnode,'HeaderLines',1);
%Fill the matrix "coord" with the x, y and z coordinates.
coord = cell2mat(getmshdata);

%Close the *.msh file
fclose(readmsh);

nflag = 5000*ones(nnode,1);

%--------------------------------------------------------------------------
%FUNCTION "getelem"
%--------------------------------------------------------------------------
function [elem,bedgeref,nflag] = getelem(filepath,nnode,nflag)
      
%Open the *.msh file
readmsh = fopen(filepath);        
%"nent" is the number of all entities understud as element in *.msh file
getmshdata = textscan(readmsh,'%u',1,'HeaderLines',7 + nnode);
%Attribute the data to "nnode"
nent = getmshdata{1};
%Swept the "elements" of "*.msh" in order to fill the parameters above. 
getmshdata = textscan(readmsh,'%*n %n %*n %n %*[^\n]',nent);
%Attribute the data to "auxmat"
auxmat = cell2mat(getmshdata);
%Fill "entitype"
entitype = auxmat(:,1);
%"nbe" verifies how many entities are edges which constitute the 
nbedge = sum(entitype == 2);
nelem = sum(entitype == 4);
%Close the *.msh file
fclose(readmsh);

%---------------------------
%Define the matrix ("elem"):
%Open again the *.msh file
readmsh = fopen(filepath);
%Create and initialize the parameter "elem"
elem = zeros(nelem,5);  %nelem rows and 5 columns
if nbedge > 0  
    H1 = 8 + nnode + sum(entitype == 15) + sum(entitype == 1);
    getmshdata = textscan(readmsh,'%*u %*u %*u %u %*u %u %u %u',nbedge,...
        'HeaderLines',H1);
    bedgeaux = cell2mat(getmshdata);
    bedgeref(1:nbedge,5) = bedgeaux(:,1); 
    bedgeref(1:nbedge,1:3) = bedgeaux(:,2:4); 
end  %End of second IF

if nelem > 0
    getmshdata = textscan(readmsh,'%*u %*u %*u %u %*u %u %u %u %u',nelem,...
        'HeaderLines',1);
    elemaux = cell2mat(getmshdata);
    %Attribute to "elem" (fifth column)
    elem(1:nelem,5) = elemaux(:,1); 
    %Attribute to "elem" (another column)
    elem(1:nelem,1:4) = elemaux(:,2:5); 
end  %End of second IF

for bi=1:size(bedgeref,1)
    vert = bedgeref(bi,1:3);
    nflag(vert) = bedgeref(bi,5);
end

%Close the *.msh file
fclose(readmsh);

%Open the *.msh file
readmsh = fopen(filepath);        
%"nent" is the number of all entities understud as element in *.msh file
getmshdata = textscan(readmsh,'%u',1,'HeaderLines',7 + nnode);
%Attribute the data to "nnode"
nent = getmshdata{1};
%Swept the "elements" of "*.msh" in order to fill the parameters above. 
getmshdata = textscan(readmsh,'%*n %n %*n %n %*[^\n]',nent);
%Attribute the data to "auxmat"
auxmat = cell2mat(getmshdata);
%Fill "entitype"
entitype = auxmat(:,1);
%"nbe" verifies how many entities are edges which constitute the 
espnode = sum(entitype == 15);
nflag(1:espnode) = auxmat(1:espnode,2);
%Close the *.msh file
fclose(readmsh);

%--------------------------------------------------------------------------
%FUNCTION "getinfoedge"
%--------------------------------------------------------------------------
function [bedge,inedge,centelem,elemvolume] = getinfo(coord,elem,bedgeref)
%
ti=1;
centelem = zeros(size(elem,1),3);
elemvolume = zeros(size(elem,1),1);
for ielem = 1:size(elem,1)
   
   %"elemnode" receives three or four nodes which constitute the element
    elemnode = elem(ielem,1:sum(elem(ielem,1:4) ~= 0));
    %Calculate the centroid
    centelem(ielem,:) = mean(coord(elemnode,:));

    I = coord(elem(ielem,1),:);
    J = coord(elem(ielem,2),:);
    K = coord(elem(ielem,3),:);
    L = coord(elem(ielem,4),:);

    IL = L - I;
    IJ = J - I;
    IK = K - I;

    elemvolume(ielem) = abs(dot(cross(IJ,IK),IL)/6);
    
    centelem(ielem,:) = centelem(ielem,:) + (elemvolume(ielem)^(1/3))*(1e-5)*rand(1,3);
    
   if ti>1
      [ vi ] = verifyface3D( inedgeaux, elem, ielem );
   else
      vi = zeros(4,1); 
   end
   ir = [1 2 3 4 1 2];
   for j=1:4
       if vi(j)==0
           inedgeaux(ti,:)=[0 0 0 0 0];
           ve=sort([elem(ielem,ir(j)) elem(ielem,ir(j+1)) elem(ielem,ir(j+2))]);
           inedgeaux(ti,1:3)=ve;
           inedgeaux(ti,4)=ielem;
           ti=ti+1;
       else
           inedgeaux(vi(j),5)=ielem;
           invec=inedgeaux(vi(j),:);
           [ invec ]=reorface3D( invec, coord, centelem, invec(4), invec(5) );
           inedgeaux(vi(j),:)=invec;   
       end
   end
   
end

[~,idx] = sort(inedgeaux(:,1));
inedgeaux = inedgeaux(idx,:);
bedgeref(:,1:3) = sort(bedgeref(:,1:3),2);
[~,idx] = sort(bedgeref(:,1));
bedgeref = bedgeref(idx,:);
inedge = inedgeaux(inedgeaux(:,5)~=0,:);
bedge = inedgeaux(inedgeaux(:,5)==0,:);

[ bedge ] = bedgeflag3D( bedge, bedgeref, coord, centelem );

%--------------------------------------------------------------------------
%FUNCTION "reorface3D"
%--------------------------------------------------------------------------
function [ matrixvec ] = reorface3D( matrixvec, coord, centelem, elem1, elem2 )
%

I=coord(matrixvec(1),:);
J=coord(matrixvec(2),:);
K=coord(matrixvec(3),:);
M=(1/3)*(I+J+K);   

JI=I-J;
JK=K-J;

v1=M-centelem(elem1,:);
c1 = dot(v1,cross(JK,JI));
if c1>0 && elem2==0
    ip=matrixvec(3);
    matrixvec(3)=matrixvec(1);
    matrixvec(1)=ip;
    matrixvec(4)=elem1;
elseif c1<0 && elem2==0
    matrixvec(4)=elem1;
elseif c1>0 && elem2~=0
    matrixvec(5)=elem1;
    matrixvec(4)=elem2;
elseif c1<0 && elem2~=0
    matrixvec(4)=elem1;
    matrixvec(5)=elem2;
end

%--------------------------------------------------------------------------
%FUNCTION "getesurn"
%--------------------------------------------------------------------------
function [ esurn1, esurn2 ] = getesurn( elem,coord )
%

t=1;
for i=1:size(coord,1)
    e1 = find(elem(:,1)==i)';
    e2 = find(elem(:,2)==i)';
    e3 = find(elem(:,3)==i)';
    e4 = find(elem(:,4)==i)';
    e = [e1 e2 e3 e4];
    k = 1;
    for j=t:t-1+size(e,2)
      esurn1(j) = e(k);
      k = k+1;
    end
    esurn2(i) = t;
    t=t+size(e,2);   
end

%--------------------------------------------------------------------------
%FUNCTION "idelemface3D"
%--------------------------------------------------------------------------
function [ elemface, innormals, bounnormals ] = idelemface3D( elem, inedge, bedge, coord, centelem )
%

elemface = zeros(size(elem,1),4);
innormals = zeros(size(inedge,1),3);
bounnormals = zeros(size(bedge,1),3);


for i=1:size(inedge,1)
    if elemface(inedge(i,4),1)==0
        j=1;
    elseif elemface(inedge(i,4),2)==0
        j=2;
    elseif elemface(inedge(i,4),3)==0
        j=3;
    else
        j=4;
    end
    elemface(inedge(i,4),j)=i; 
    if elemface(inedge(i,5),1)==0
        j=1;
    elseif elemface(inedge(i,5),2)==0
        j=2;
    elseif elemface(inedge(i,5),3)==0
        j=3;
    else
        j=4;
    end
    elemface(inedge(i,5),j)=i; 
    
    L = centelem(inedge(i,4),:); R = centelem(inedge(i,5),:);
    v1 = coord(inedge(i,2),:)-coord(inedge(i,1),:);
    v2 = coord(inedge(i,3),:)-coord(inedge(i,1),:);
    n = 0.5*cross(v1,v2); LR = R - L;
    innormals(i,:) = sign(dot(n,LR))*n;
    
end

ni = size(inedge,1);

for i=1:size(bedge,1)
    if elemface(bedge(i,4),1)==0
        j=1;
    elseif elemface(bedge(i,4),2)==0
        j=2;
    elseif elemface(bedge(i,4),3)==0
        j=3;
    else
        j=4;
    end
    elemface(bedge(i,4),j)=i+ni;
    
    L = centelem(bedge(i,4),:); I = coord(bedge(i,1),:);
    v1 = coord(bedge(i,2),:)-coord(bedge(i,1),:);
    v2 = coord(bedge(i,3),:)-coord(bedge(i,1),:);
    n = 0.5*cross(v1,v2); LI = I - L;
    bounnormals(i,:) = sign(dot(n,LI))*n;
    
end

%--------------------------------------------------------------------------
%FUNCTION "bedgeflag3D"
%--------------------------------------------------------------------------
function [ bedge ] = bedgeflag3D( bedge, bedgeref, coord, centelem )
%

for i=1:size(bedge,1)
    
    invec = bedge(i,:);
    [ invec ] = reorface3D( invec, coord, centelem, invec(4), invec(5) );
    bedge(i,:) = invec;
    
    edge = sort(bedge(i,1:3));
    
    linhasr = find(bedgeref(:,1)==edge(1));

    for j=linhasr'
        if bedgeref(j,1)==edge(1) && bedgeref(j,2)==edge(2) && bedgeref(j,3)==edge(3)
            bedge(i,5)=bedgeref(j,5);
        end
    end
    
end

%--------------------------------------------------------------------------
%FUNCTION "verifyface3D"
%--------------------------------------------------------------------------
function [ vi ] = verifyface3D( inedgeaux, elem, i )
%
vi=zeros(4,1);

v1=sort([elem(i,1) elem(i,2) elem(i,3)]);
v2=sort([elem(i,2) elem(i,3) elem(i,4)]);
v3=sort([elem(i,3) elem(i,4) elem(i,1)]);
v4=sort([elem(i,4) elem(i,1) elem(i,2)]);

linhas1 = find(inedgeaux(:,1)==v1(1));
linhas2 = find(inedgeaux(:,1)==v2(1));
linhas3 = find(inedgeaux(:,1)==v3(1));
linhas4 = find(inedgeaux(:,1)==v4(1));

if isempty(linhas1)==0
    for j=linhas1'
        if inedgeaux(j,1)==v1(1) && inedgeaux(j,2)==v1(2) && inedgeaux(j,3)==v1(3)
            vi(1)=j;
        end
    end
end
if isempty(linhas2)==0
    for j=linhas2'
        if inedgeaux(j,1)==v2(1) && inedgeaux(j,2)==v2(2) && inedgeaux(j,3)==v2(3)
            vi(2)=j;
        end
    end
end
if isempty(linhas3)==0
    for j=linhas3'
        if inedgeaux(j,1)==v3(1) && inedgeaux(j,2)==v3(2) && inedgeaux(j,3)==v3(3)
            vi(3)=j;
        end
    end
end
if isempty(linhas4)==0
    for j=linhas4'
        if inedgeaux(j,1)==v4(1) && inedgeaux(j,2)==v4(2) && inedgeaux(j,3)==v4(3)
            vi(4)=j;
        end
    end
end

%--------------------------------------------------------------------------
%FUNCTION "readopt3D"
%--------------------------------------------------------------------------
function [ caso, malha, tipopeso, calcpeso, limitadora, fractarch, bcflag ] = readopt3D
%
opt = fopen('Options3D.txt','r'); 
celcaso = textscan(opt,'%u',1,'headerlines',4); caso = celcaso{1};
celmalha = textscan(opt,'%s',1,'headerlines',3); malha = char(celmalha{1});
celtipopeso = textscan(opt,'%s',1,'headerlines',3); tipopeso = char(celtipopeso{1});
celcalcpeso = textscan(opt,'%s',1,'headerlines',3); calcpeso = char(celcalcpeso{1});
cellimitadora = textscan(opt,'%s',1,'headerlines',3); limitadora = char(cellimitadora{1});
celfractarch = textscan(opt,'%s',1,'headerlines',6); fractarch = char(celfractarch{1});
bcflag = [101 0; 102 1; 103 2; 201 0; 202 -1; 203 1];
fclose(opt);

%--------------------------------------------------------------------------
%FUNCTION "assemblyfract"
%--------------------------------------------------------------------------
function sist = assemblyfract( sist, fracture, element, options, face )
%
for col=1:size(element.areafrac,2)
    v1 = fracture.coord{col}(2,:) - fracture.coord{col}(1,:);
    v2 = fracture.coord{col}(3,:) - fracture.coord{col}(1,:);
    n = cross(v1,v2); n = n/norm(n);
    linhas = find(element.areafrac(:,col)~=0);
    for lin = linhas'
        reg = element.region(lin);
        cent = element.centroid(lin,:);
        Km = options.tensor{reg}(cent(1),cent(2),cent(3));
        Kf = fracture.permeability{fracture.type(col)};
        Tmf = 2*(n*Km*n')*element.areafrac(lin,col)/element.distfrac(lin,col);
        Tfm = 2*(n*Kf*n')*element.areafrac(lin,col)/fracture.aperture(fracture.type(col));
        Mcol = col + size(element.vertices,1); TF = ( Tmf^-1 + Tfm^-1 )^-1;
        sist.Mtpfa(lin,lin) = sist.Mtpfa(lin,lin) + TF;
        sist.Mtpfa(lin,Mcol) = sist.Mtpfa(lin,Mcol) - TF;
        sist.Mtpfa(Mcol,lin) = sist.Mtpfa(Mcol,lin) - TF;
        sist.Mtpfa(Mcol,Mcol) = sist.Mtpfa(Mcol,Mcol) + TF;
    end 
end

for lin=1:size(fracture.connecwithfract,1)
    v11 = fracture.coord{lin}(2,:) - fracture.coord{lin}(1,:);
    v12 = fracture.coord{lin}(3,:) - fracture.coord{lin}(1,:);
    n1 = cross(v11,v12); n1 = n1/norm(n1);
    Kf1 = fracture.permeability{fracture.type(lin)};
    colunas = find(fracture.connecwithfract(lin,:)~=0);
    for col = colunas
        v21 = fracture.coord{col}(2,:) - fracture.coord{col}(1,:);
        v22 = fracture.coord{col}(3,:) - fracture.coord{col}(1,:);
        n2 = cross(v21,v22); n2 = n2/norm(n2);
        Kf2 = fracture.permeability{fracture.type(col)};
        Tf1 = (n1*Kf1*n1')*fracture.connecwithfract(lin,col)*fracture.aperture(fracture.type(lin));
        Tf2 = (n2*Kf2*n2')*fracture.connecwithfract(col,lin)*fracture.aperture(fracture.type(col));
        TF = ( Tf1^-1 + Tf2^-1 )^-1;
        Mcol = col + size(element.vertices,1);
        Mlin = lin + size(element.vertices,1);
        sist.Mtpfa(Mcol,Mcol) = sist.Mtpfa(Mcol,Mcol) + TF;
        sist.Mtpfa(Mcol,Mlin) = sist.Mtpfa(Mcol,Mlin) - TF;
        sist.Mtpfa(Mlin,Mcol) = sist.Mtpfa(Mlin,Mcol) - TF;
        sist.Mtpfa(Mlin,Mlin) = sist.Mtpfa(Mlin,Mlin) + TF;
    end 
end

for lin=1:size(fracture.slicingbound,1)
    v11 = fracture.coord{lin}(2,:) - fracture.coord{lin}(1,:);
    v12 = fracture.coord{lin}(3,:) - fracture.coord{lin}(1,:);
    n1 = cross(v11,v12); n1 = n1/norm(n1);
    Kf1 = fracture.permeability{fracture.type(lin)};
    colunas = find(fracture.slicingbound(lin,:)~=0);
    for col = colunas
        if face.bound.flag(col)>200
            bcflag = options.flagcorresp;   
            c1 = find(bcflag(:,1)==face.bound.flag(col));
            if bcflag(c1,1)>200, g = bcflag(c1,2); end
            Mlin = lin + size(element.vertices,1);
            sist.Btpfa(Mlin,1) = sist.Btpfa(Mlin,1) + g*fracture.slicingbound(lin,col)*fracture.aperture(fracture.type(lin));
        else
            Tf1 = (n1*Kf1*n1')*fracture.slicingbound(lin,col)*fracture.aperture(fracture.type(lin))/fracture.disttobound(lin,col);
            Mlin = lin + size(element.vertices,1);
            sist.Mtpfa(Mlin,Mlin) = sist.Mtpfa(Mlin,Mlin) + Tf1;
            sist.Btpfa(Mlin,1) = sist.Btpfa(Mlin,1) + Tf1*fracture.valueonbound(lin,col);
        end 
    end
end


    
    








