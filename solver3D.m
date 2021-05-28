function p = solver3D( wells, Keq, GV )
%
global sist element options

s = determinasource; element.sourceterm = s;

sist.Btpfa = sist.Btpfa + s;

if strcmp(options.limitadora,'sim')
    
    p = slopelim( wells, Keq, GV );
%     p = slopetestinho( wells, Keq, GV );
    
else

    M = sist.Mtpfa + sist.Mcdt;
    B = sist.Btpfa + sist.Bcdt;

    [ M, B ] = incluiwells( M, B, wells );

    p = M\B;

end

end

function s = determinasource
%
global options

if strcmp(options.limitadora,'sim')
    typewei = [options.tipopeso '_' 'Limited'];
else
    typewei = options.tipopeso;
end

[~,nameso,~] = fileparts(options.malha);   
if exist(sprintf('3D - Caso %u\\%s\\%s\\source.mat',options.caso,typewei,nameso),'file') == 2
    load(sprintf('3D - Caso %u\\%s\\%s\\source',options.caso,typewei,nameso));
else
    [ s ] = sourceterm3D;
    if exist(sprintf('3D - Caso %u\\%s\\%s',options.caso,typewei,nameso),'dir')==0
        mkdir(sprintf('3D - Caso %u\\%s\\%s',options.caso,typewei,nameso));
    end
    save(sprintf('3D - Caso %u\\%s\\%s\\source',options.caso,typewei,nameso),'s');
end

end

function [ M, B ] = incluiwells( M, B, wells )
%

    for iw=1:size(wells,1)      
        if (wells(iw,1)~=0)&&(wells(iw,5)>400)&&(wells(iw,5)<600) % Caso seja prescrita alguma pressão.
            M(wells(iw,2),:) = 0*M(wells(iw,2),:);
            M(wells(iw,2),wells(iw,2)) = 1;
            B(wells(iw,2)) = wells(iw,6);
        end
    end

end






