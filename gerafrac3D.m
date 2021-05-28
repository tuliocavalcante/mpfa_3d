function [ fracture ] = gerafrac3D(options)
%

if strcmp(options.fractarch,'novo')==1

    addpath('ADFNE15');

    [num, text, raw] = xlsread('ADFNE1.5_Parameters.xlsx',1,'C2:C4');
    dim = num(1); nfractype = num(2); nome = char(text);

    coordfrac = [];
    typefrac = [];

    for i=1:nfractype

        typeparam1 = xlsread('ADFNE1.5_Parameters.xlsx',i+1,'C2:C6');
        region = xlsread('ADFNE1.5_Parameters.xlsx',i+1,'C7:H7');
        typeparam2 = xlsread('ADFNE1.5_Parameters.xlsx',i+1,'C8:C10');

        nfrac(i) = typeparam1(1);
        meangleori(i) = typeparam1(2);
        varmeangleori(i) = typeparam1(3);
        minlen(i) = typeparam1(4);
        maxlen(i) = typeparam1(5);
        meangleprof(i) = typeparam2(1);
        varmeangleprof(i) = typeparam2(2);
        numvertfrac(i) = typeparam2(3);

        if numvertfrac(i)==0
            set = Field(DFN('dim',dim,'n',nfrac(i),'dir',meangleori(i),'ddir',varmeangleori(i),'minl',minlen(i),... 
            'mu',0.1,'maxl',maxlen(i),'bbx',region,'dip',typeparam2(1),'ddip',typeparam2(2)),'Poly');
        else
            set = Field(DFN('dim',dim,'n',nfrac(i),'dir',meangleori(i),'ddir',varmeangleori(i),'minl',minlen(i),...
            'mu',0.1,'maxl',maxlen(i),'bbx',region,'dip',typeparam2(1),'ddip',typeparam2(2),'shape','e','q',typeparam2(3)),'Poly');  
        end

        coordfrac = [coordfrac; set];
        typefrac = [typefrac; i*ones(size(set))];

    end

    fid = fopen(sprintf('ADFNE15\\%s.txt',nome),'w');
    fprintf(fid,'Número de Fraturas:\n');
    fprintf(fid,'%u\n',size(coordfrac,1));
    fprintf(fid,'Classificação:\n');
    fprintf(fid,'%u\n',typefrac);
    fprintf(fid,'Fraturas:\n');
    for i=1:size(coordfrac,1)
        fprintf(fid,'%u\n',size(coordfrac{i},1));
        fprintf(fid,'%f %f %f\n',coordfrac{i}');
    end
    fclose(fid);
    
    fracture.coord = coordfrac;
    fracture.type = typefrac;

    rmpath('ADFNE15');
    
    fid = fopen('ADFNE15\\OtherOptions.txt','r');
    celnumtype = textscan(fid,'%u',1,'headerlines',1); numtype = celnumtype{1};
    nomeperm = textscan(fid,'%u',1,'headerlines',1);
    for i=1:numtype
        celnumperm = textscan(fid,'%u',1,'headerlines',1);
        celperm = textscan(fid,'%f %f %f',3,'headerlines',1); fracture.permeability{i,1} = [celperm{1} celperm{2} celperm{3}];
    end
    celaperture = textscan(fid,'%f',numtype,'headerlines',3); fracture.aperture = celaperture{1};
    fclose(fid);
    
elseif strcmp(options.fractarch,'nenhum')==1
    
    coordfrac = {};
    fracture.coord = coordfrac;
    
else
    
    fid = fopen(sprintf('ADFNE15\\%s.txt',options.fractarch),'r');
    celnumfrac = textscan(fid,'%u',1,'headerlines',1); numfrac = celnumfrac{1};
    celtypefrac = textscan(fid,'%u',numfrac,'headerlines',2); typefrac = celtypefrac{1};
    celnome = textscan(fid,'%s',1,'headerlines',1);
    
    for i=1:numfrac
        celnumvert = textscan(fid,'%u',1,'headerlines',1); numvert = celnumvert{1};
        celvert = textscan(fid,'%f %f %f',numvert,'headerlines',1); coordfrac{i,1} = [celvert{1} celvert{2} celvert{3}];
    end
    
    fracture.coord = coordfrac;
    fracture.type = typefrac;

    fclose(fid);
    
    fid = fopen('ADFNE15\\OtherOptions.txt','r');
    celnumtype = textscan(fid,'%u',1,'headerlines',1); numtype = celnumtype{1};
    nomeperm = textscan(fid,'%u',1,'headerlines',1);
    for i=1:numtype
        celnumperm = textscan(fid,'%u',1,'headerlines',1);
        celperm = textscan(fid,'%f %f %f',3,'headerlines',1); fracture.permeability{i,1} = [celperm{1} celperm{2} celperm{3}];
    end
    celaperture = textscan(fid,'%f',numtype,'headerlines',3); fracture.aperture = celaperture{1};
    fclose(fid);
    
end

end

