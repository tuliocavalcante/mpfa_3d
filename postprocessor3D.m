function postprocessor3D( p, idxvpi )
% Create a vtk file for solution visualization ---------------------------% 
global options element vertex fracture

%--------------------------------------------------------------------------
fname = 'Output'; ext = '.vtk';

nnode = size(vertex.coord,1); nelem = size(element.vertices,1);
type_elem = zeros(nelem,1); ivpi = num2str(idxvpi);
sumnodebyelem = sum(sum(element.vertices(1:nelem,1:4) ~= 0));
nfrac = size(fracture.coord,1);
%--------------------------------------------------------------------------

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

fname_vtk = [fname '_' ivpi ext]; [~,nameso,~] = fileparts(options.malha);

if strcmp(options.limitadora,'sim')
    typewei = [options.tipopeso '_' 'Limited'];
else
    typewei = options.tipopeso;
end

if exist(sprintf('3D - Caso %u\\%s\\%s\\Results',options.caso,typewei,nameso),'dir')==0
    mkdir(sprintf('3D - Caso %u\\%s\\%s\\Results',options.caso,typewei,nameso));
end
fid = fopen(sprintf('3D - Caso %u\\%s\\%s\\Results\\%s',options.caso,typewei,nameso,fname_vtk),'w'); 
%Write head informations
fprintf(fid,'# vtk DataFile Version 2.0\n');
fprintf(fid,'Pressure Field Data \r\n');
fprintf(fid,'ASCII \r\n');
%Information about grid type
fprintf(fid,'DATASET UNSTRUCTURED_GRID \r\n\r\n');

%Write the POINT informations:
%Head (POINT)
fprintf(fid,'POINTS %i float \r\n\r\n',nnode);
%Distribution (POINT)
data1 = [vertex.coord(:,1:3)'];
%Print the distribution
fprintf(fid,'%26.16E %26.16E %26.16E \r\n',data1);
%Jump a line
fprintf(fid,'\r\n');

%Write the CELL informations:
%Head (CELL)
fprintf(fid,'CELLS %i %i \r\n\r\n',nelem,sumnodebyelem+nelem);

nodebyelem = 4;
%"numdata" is the first column of data to be printed in the CELL sect.
numdata = nodebyelem*ones(nelem,1);

%Considering which the initial data (counter) begining from zero, the 
%node number is diminish of one (-1) in the node definition (command 
%below). The decision command below is used in order to define if the 
%elemnts ploted is a triangle ("nodebyelem" == 3) or quadrangle 
%("nodebyelem" == 4)
data2 = [numdata';(element.vertices(1:nelem,1)-1)';...
        (element.vertices(1:nelem,2)-1)';...
        (element.vertices(1:nelem,3)-1)';...
        (element.vertices(1:nelem,4)-1)'];
fprintf(fid,'%i %i %i %i %i \r\n',data2);
%Definition of element type (10 is to tetrahedro)
type_elem(1:nelem) = 10;

%Jump a line
fprintf(fid,'\r\n');

%Cell type information
fprintf(fid,'CELL_TYPES %i \r\n\r\n',nelem);
fprintf(fid,'%i \r\n',type_elem);
%Jump a line
fprintf(fid,'\r\n');

%Write the CELL_DATA informations:
%Head (CELL_DATA)
fprintf(fid,'CELL_DATA %i \r\n',nelem);
%Write data related to PRESSURE
fprintf(fid,'SCALARS Pressure float 1 \r\n');
fprintf(fid,'LOOKUP_TABLE default \r\n\r\n');
fprintf(fid,'%26.16E \r\n',full(p));

%Jump a line
fprintf(fid,'\r\n');

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

if size(fracture.coord,1)>0

fname_vtk = [fname '_Fract_' ivpi ext]; [~,nameso,~] = fileparts(options.malha);

if strcmp(options.limitadora,'sim')
    typewei = [options.tipopeso '_' 'Limited'];
else
    typewei = options.tipopeso;
end

if exist(sprintf('3D - Caso %u\\%s\\%s\\Results',options.caso,typewei,nameso),'dir')==0
    mkdir(sprintf('3D - Caso %u\\%s\\%s\\Results',options.caso,typewei,nameso));
end
fid = fopen(sprintf('3D - Caso %u\\%s\\%s\\Results\\%s',options.caso,typewei,nameso,fname_vtk),'w'); 
%Write head informations
fprintf(fid,'# vtk DataFile Version 2.0\n');
fprintf(fid,'Pressure Field Data \r\n');
fprintf(fid,'ASCII \r\n');
%Information about grid type
fprintf(fid,'DATASET UNSTRUCTURED_GRID \r\n\r\n');


%Write the POINT informations:
nfnode = 0;
for i=1:nfrac
    nfnode = nfnode + size(fracture.coord{i},1);
end
%Head (POINT)
fprintf(fid,'POINTS %i float \r\n\r\n',nfnode);
%Distribution (POINT)
sumnodebyfrac = 0; nelfrac = 0; K = []; pfrac = []; tfrac = []; ref = 0;
for i=1:nfrac
    data1 = [fracture.coord{i}(:,1:3)']; 
    fprintf(fid,'%26.16E %26.16E %26.16E \r\n',data1);
    k = convhull(fracture.coord{i}(:,1),fracture.coord{i}(:,2),fracture.coord{i}(:,3)) + ref;
    K = [K;k];
    sumnodebyfrac = sumnodebyfrac + 3*size(k,1);
    nelfrac = nelfrac + size(k,1); ref = ref + size(fracture.coord{i},1);
    pfrac = [pfrac; p(size(element.vertices,1)+i)*ones(size(k,1),1)];
%     tfrac = [tfrac; double(fracture.type(i))*ones(size(k,1),1)];
    tfrac = [tfrac; fracture.permeability{fracture.type(i)}(1,1)*ones(size(k,1),1)];
end

%Jump a line
fprintf(fid,'\r\n');

%Write the CELL informations:
%Head (CELL)
fprintf(fid,'CELLS %i %i \r\n\r\n',nelfrac,sumnodebyfrac+nelfrac);

nodebyelem = 3;
%"numdata" is the first column of data to be printed in the CELL sect.
numdata = nodebyelem*ones(nelfrac,1);

%Considering which the initial data (counter) begining from zero, the 
%node number is diminish of one (-1) in the node definition (command 
%below). The decision command below is used in order to define if the 
%elemnts ploted is a triangle ("nodebyelem" == 3) or quadrangle 
%("nodebyelem" == 4)
data2 = [numdata';(K(1:nelfrac,1)-1)';...
                  (K(1:nelfrac,2)-1)';...
                  (K(1:nelfrac,3)-1)'];
fprintf(fid,'%i %i %i %i \r\n',data2);
%Definition of element type (5 is a code to triangle)
type_frac(1:nelfrac) = 5;

%Jump a line
fprintf(fid,'\r\n');

%Cell type information
fprintf(fid,'CELL_TYPES %i \r\n\r\n',nelfrac);
fprintf(fid,'%i \r\n',type_frac);
%Jump a line
fprintf(fid,'\r\n');

%Write the CELL_DATA informations:
%Head (CELL_DATA)
fprintf(fid,'CELL_DATA %i \r\n',nelfrac);
%Write data related to PRESSURE
fprintf(fid,'SCALARS Pressure float 1 \r\n');
fprintf(fid,'LOOKUP_TABLE default \r\n\r\n');
fprintf(fid,'%26.16E \r\n',full(pfrac));

%Jump a line
fprintf(fid,'\r\n');

%Write data related to FRAC TYPE
fprintf(fid,'SCALARS Permeability float 1 \r\n');
fprintf(fid,'LOOKUP_TABLE default \r\n\r\n');
fprintf(fid,'%26.16E \r\n',full(tfrac));

%Jump a line
fprintf(fid,'\r\n');
    
end

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

fclose(fid);

end  %End of function

