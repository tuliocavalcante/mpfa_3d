clear
clc

nel = 10; npo = nel + 1; d = 1/nel;

coord = zeros(npo^3,3);

t = 1;
for i=1:npo
    for j=1:npo
        for k=1:npo
            
            coord(t,:) = [d*(k-1) d*(j-1) d*(i-1)];
            t = t + 1;
            
            if k<npo && j<npo
                coord(t,:) = [d*(k-1)+0.5*d d*(j-1)+0.5*d d*(i-1)];
                t = t + 1;
            end
            
            if k<npo && i<npo
                coord(t,:) = [d*(k-1)+0.5*d d*(j-1) d*(i-1)+0.5*d];
                t = t + 1;
            end
            
            if i<npo && j<npo
                coord(t,:) = [d*(k-1) d*(j-1)+0.5*d d*(i-1)+0.5*d];
                t = t + 1;
            end
                
            if k<npo && j<npo && i<npo
                coord(t,:) = [d*(k-1)+0.5*d d*(j-1)+0.5*d d*(i-1)+0.5*d];
                t = t + 1;
            end
            
        end
    end
end

sprintf('Coordenadas geradas! São %u pontos!',size(coord,1))

DT = delaunayn(coord);
n1 = size(DT,1);
tetras = zeros(size(DT,1),9);
tetras(:,6:9) = DT;
tetras(:,1) = 1:n1; tetras(:,2) = 4; tetras(:,3) = 2; 
tetras(:,4) = 1; tetras(:,5) = 1; 

sprintf('Tetraedros gerados! São %u células!',size(tetras,1))

t = 1;
for i=1:size(tetras,1)
    
    allfaces(t,:) = sort([tetras(i,6) tetras(i,7) tetras(i,8)]); t = t + 1;
    allfaces(t,:) = sort([tetras(i,7) tetras(i,8) tetras(i,9)]); t = t + 1;
    allfaces(t,:) = sort([tetras(i,8) tetras(i,9) tetras(i,6)]); t = t + 1;
    allfaces(t,:) = sort([tetras(i,9) tetras(i,6) tetras(i,7)]); t = t + 1;
    
end

b = 1;
for i=1:size(allfaces,1)
    
    linhas = find(allfaces(:,1)==allfaces(i,1));
    
    if isempty(linhas)==0 && allfaces(i,1)~=0
    
        facei = allfaces(i,:);

        [~, index] = ismember(allfaces(linhas,:),facei,'rows');

        teste = find(index);

        if size(teste,1)==1
           bound(b,:) = allfaces(linhas(teste),:);
           allfaces(linhas(teste),:) = zeros(size(allfaces(linhas(teste),:)));
           b = b + 1;
        else
           allfaces(linhas(teste),:) = zeros(size(allfaces(linhas(teste),:)));
        end
        
    end
    
    andamento = i/size(allfaces,1)
    
end

faces = zeros(size(bound,1),8); n2 = size(bound,1);
faces(:,1) = 1:n2; faces(:,2) = 2; faces(:,3) = 2; 
faces(:,4) = 101; faces(:,5) = 1; faces(:,6:8) = bound; 
tetras(:,1) = n2+1:n2+n1;
pontos = zeros(size(coord,1),4);
pontos(:,1) = 1:size(coord,1);
pontos(:,2:4) = coord;

fid = fopen(sprintf('C:\\Users\\Túlio\\OneDrive\\Documentos\\doutorado\\mpfa3d\\Malhas\\cubodeyang_%u.msh',4),'w');
fprintf(fid,'$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n');
fprintf(fid,sprintf('%u\n',size(coord,1)));
fprintf(fid,sprintf('%u %f %f %f\n',pontos'));
fprintf(fid,'$EndNodes\n$Elements\n');
fprintf(fid,sprintf('%u\n',size(tetras,1)+size(bound,1)));
fprintf(fid,sprintf('%u %u %u %u %u %u %u %u\n',faces'));
fprintf(fid,sprintf('%u %u %u %u %u %u %u %u %u\n',tetras'));
fprintf(fid,'$EndElements\n');
fclose(fid);


