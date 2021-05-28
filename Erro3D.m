function [ eL2, gradL2 ] = Erro3D( p, G, NG )
%
global options element

a = zeros(size(element.vertices,1),1); p = p(1:size(element.vertices,1));

for i=1:size(element.vertices,1)
    
    reg = element.region(i);
    a(i,1) = options.solanalit{reg}(element.centroid(i,1),element.centroid(i,2),element.centroid(i,3));
    gradif(i,1) = norm(G(i,:) - NG(i,:));
    aG(i,1) = norm(G(i,:));
    
    if options.caso == 10 && i == 1, a(i) = options.solanalit{2}(0,0,0); end 
    
end

eL2 = sqrt(sum((abs(a-p).^2).*element.volume)./sum((a.^2).*element.volume));

%--------------------------------------------------------------------------
% Sandve (2012)
% M = max(a)-min(a);
% eL2 = sqrt(er1)/(M*sum(elemvolume));
%--------------------------------------------------------------------------

gradL2 = sqrt(sum((gradif.^2).*element.volume)./sum((aG.^2).*element.volume));

eOU = sqrt(sum(((max(0,p-max(a)).^2).*element.volume))+...
           sum(((max(0,min(a)-p).^2).*element.volume)));

if strcmp(options.limitadora,'sim')
    typewei = [options.tipopeso '_' 'Limited'];
else
    typewei = options.tipopeso;
end
       
[~,nameso,~] = fileparts(options.malha);
if exist(sprintf('3D - Caso %u\\%s\\%s',options.caso,typewei,nameso),'dir')==0
    mkdir(sprintf('3D - Caso %u\\%s\\%s',options.caso,typewei,nameso));
end
save(sprintf('3D - Caso %u\\%s\\%s\\eL2',options.caso,typewei,nameso),'eL2');
save(sprintf('3D - Caso %u\\%s\\%s\\eOU',options.caso,typewei,nameso),'eOU');
save(sprintf('3D - Caso %u\\%s\\%s\\gradL2',options.caso,typewei,nameso),'gradL2');

end

