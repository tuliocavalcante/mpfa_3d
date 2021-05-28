clear
clc

casoi = 7; casof = 7;
malhai = 1; malhaf = 6;
interpi = 3; interpf = 3;
nomesh = 'benchtetra_';
% nomesh = 'cubodeyang_';

ErroL2 = zeros(malhaf-malhai+1,casof-casoi+1,12);
GradL2 = zeros(malhaf-malhai+1,casof-casoi+1,12);
N = zeros(malhaf-malhai+1,casof-casoi+1);
Re = zeros(malhaf-malhai+1,casof-casoi+1);
Rg = zeros(malhaf-malhai+1,casof-casoi+1);

for caso=casoi:casof
    
    for tw = [3 9]
        
        if tw==1, typewei = 'LPEW1_Limited'; end
        if tw==2, typewei = 'LPEW2_Limited'; end
        if tw==3, typewei = 'LPEW3_Limited'; end
        if tw==4, typewei = 'LSW_Limited'; end
        if tw==5, typewei = 'IDW_Limited'; end
        if tw==6, typewei = 'YG2019_Limited'; end
        if tw==7, typewei = 'LPEW1'; end
        if tw==8, typewei = 'LPEW2'; end
        if tw==9, typewei = 'LPEW3'; end
        if tw==10, typewei = 'LSW'; end
        if tw==11, typewei = 'IDW'; end
        if tw==12, typewei = 'YG2019'; end
        
        for malha=malhai:malhaf  
            
            nome = sprintf('%s%u.msh',nomesh,malha);
            [~,nameso,~] = fileparts(nome);

            load(sprintf('3D - Caso %u\\%s\\%s\\eL2',caso,typewei,nameso));
            ErroL2(malha,caso,tw) = eL2;

            load(sprintf('3D - Caso %u\\%s\\%s\\gradL2',caso,typewei,nameso));
            GradL2(malha,caso,tw) = gradL2;

            load(sprintf('3D - Caso %u\\%s\\%s\\source',caso,typewei,nameso));
            N(malha,caso,tw) = size(s,1);

        end
        
    end 
    
end

for caso=1:casof
    
    for tw = [3 9]
        
        if tw==1, typewei = 'LPEW1_Limited'; end
        if tw==2, typewei = 'LPEW2_Limited'; end
        if tw==3, typewei = 'LPEW3_Limited'; end
        if tw==4, typewei = 'LSW_Limited'; end
        if tw==5, typewei = 'IDW_Limited'; end
        if tw==6, typewei = 'YG2019_Limited'; end
        if tw==7, typewei = 'LPEW1'; end
        if tw==8, typewei = 'LPEW2'; end
        if tw==9, typewei = 'LPEW3'; end
        if tw==10, typewei = 'LSW'; end
        if tw==11, typewei = 'IDW'; end
        if tw==12, typewei = 'YG2019'; end
        
        for malha=2:malhaf

            if ErroL2(malha-1,caso,tw)>1e-12 

                Re(malha,caso,tw) = -3*(log(ErroL2(malha,caso,tw)/ErroL2(malha-1,caso,tw))/...
                                log(N(malha,caso,tw)/N(malha-1,caso,tw)));

                Rg(malha,caso,tw) = -3*(log(GradL2(malha,caso,tw)/GradL2(malha-1,caso,tw))/...
                                log(N(malha,caso,tw)/N(malha-1,caso,tw)));

            end

        end
        
    end
    
end

for caso=casoi:casof
           
    for tw = [3 9]
        
        if tw==1, typewei = 'LPEW1_Limited'; end
        if tw==2, typewei = 'LPEW2_Limited'; end
        if tw==3, typewei = 'LPEW3_Limited'; end
        if tw==4, typewei = 'LSW_Limited'; end
        if tw==5, typewei = 'IDW_Limited'; end
        if tw==6, typewei = 'YG2019_Limited'; end
        if tw==7, typewei = 'LPEW1'; end
        if tw==8, typewei = 'LPEW2'; end
        if tw==9, typewei = 'LPEW3'; end
        if tw==10, typewei = 'LSW'; end
        if tw==11, typewei = 'IDW'; end
        if tw==12, typewei = 'YG2019'; end
        
        fid = fopen(sprintf('3D - Caso %u\\%s\\TabeladeResultados.txt',caso,typewei),'w');
                
        for malha=malhai:malhaf
            
            fprintf(fid,'Caso %u\n',caso);
            fprintf(fid,'%u %1.3f %1.3f %1.3f %1.3f\n',N(malha,caso,tw),ErroL2(malha,caso,tw),Re(malha,caso,tw),GradL2(malha,caso,tw),Rg(malha,caso,tw));
        
        end
        
        fclose(fid);
    
    end
    
    minerr = 1e50; 
    mingrd = 1e50;
    for i = [3 9]
        minerr = min(minerr,ErroL2(end,caso,i));
        mingrd = min(mingrd,GradL2(end,caso,i));
    end

    n = N(:,caso,tw);
    refx = [min(n) 10*min(n)];
    refy = [minerr minerr*exp((-2/3)*log(10))];

    rgfx = [min(n) 10*min(n)];
    rgfy = [mingrd mingrd*exp((-1/3)*log(10))];

%     figure(1)
%     loglog(n,ErroL2(:,caso,1),'-g',n,ErroL2(:,caso,2),'-k',n,ErroL2(:,caso,3),'-m',n,ErroL2(:,caso,4),'-r',n,ErroL2(:,caso,5),'-b',...
%            n,ErroL2(:,caso,6),'-y',n,ErroL2(:,caso,7),'--^g',n,ErroL2(:,caso,8),'--*k',n,ErroL2(:,caso,9),'--dm',...
%            n,ErroL2(:,caso,10),'--sr',n,ErroL2(:,caso,11),'--ob',n,ErroL2(:,caso,12),'--+y',refx,refy,'-k');
%     legend('LPEW1 Limited','LPEW2 Limited','LPEW3 Limited','LSW Limited','IDW Limited','YG2019 Limited','LPEW1','LPEW2','LPEW3','LSW','IDW','YG2019');
%     xlabel('\it\fontname{Times} n');
%     ylabel('\it\fontname{Times} l^{2}_{u}');
%     title(sprintf('Caso %u - Escalar',caso));
%     
%     saveas(gcf,sprintf('3D - Caso %u\\ConvergenceEL2.png',caso));
%     
%     figure(2)
%     loglog(n,GradL2(:,caso,1),'-g',n,GradL2(:,caso,2),'-k',n,GradL2(:,caso,3),'-m',n,GradL2(:,caso,4),'-r',n,GradL2(:,caso,5),'-b',...
%            n,GradL2(:,caso,6),'-y',n,GradL2(:,caso,7),'--^g',n,GradL2(:,caso,8),'--*k',n,GradL2(:,caso,9),'--dm',...
%            n,GradL2(:,caso,10),'--sr',n,GradL2(:,caso,11),'--ob',n,GradL2(:,caso,12),'--+y',rgfx,rgfy,'-k');
%     legend('LPEW1 Limited','LPEW2 Limited','LPEW3 Limited','LSW Limited','IDW Limited','YG2019 Limited','LPEW1','LPEW2','LPEW3','LSW','IDW','YG2019');
%     xlabel('\it\fontname{Times} n');
%     ylabel('\it\fontname{Times} l^{2}_{u}');
%     title(sprintf('Caso %u - Gradiente',caso));
%     
%     saveas(gcf,sprintf('3D - Caso %u\\ConvergenceGL2.png',caso));
    
end


% figure(1)
% loglog(n,ErroL2(:,caso,1),'-g',n,ErroL2(:,caso,2),'-k',n,ErroL2(:,caso,3),'-m',n,ErroL2(:,caso,4),'-r',n,ErroL2(:,caso,5),'-b',...
%        n,ErroL2(:,caso,6),'^g',n,ErroL2(:,caso,7),'*k',n,ErroL2(:,caso,8),'dm',n,ErroL2(:,caso,9),'sr',n,ErroL2(:,caso,10),'ob',refx,refy,'-k');
% legend('LPEW1 Limited','LPEW2 Limited','LPEW3 Limited','IDW Limited','LSW Limited','LPEW1','LPEW2','LPEW3','IDW','LSW');
% xlabel('\it\fontname{Times} n');
% ylabel('\it\fontname{Times} l^{2}_{u}');
% title(sprintf('Caso %u - Escalar',caso));
% 
% saveas(gcf,sprintf('3D - Caso %u\\ConvergenceEL2.png',caso));
% 
% figure(2)
% loglog(n,GradL2(:,caso,1),'-g',n,GradL2(:,caso,2),'-k',n,GradL2(:,caso,3),'-m',n,GradL2(:,caso,4),'-r',n,GradL2(:,caso,5),'-b',...
%        n,GradL2(:,caso,6),'^g',n,GradL2(:,caso,7),'*k',n,GradL2(:,caso,8),'dm',n,GradL2(:,caso,9),'sr',n,GradL2(:,caso,10),'ob',refx,refy,'-k');
% legend('LPEW1 Limited','LPEW2 Limited','LPEW3 Limited','IDW Limited','LSW Limited','LPEW1','LPEW2','LPEW3','IDW','LSW');
% xlabel('\it\fontname{Times} n');
% ylabel('\it\fontname{Times} l^{2}_{u}');
% title(sprintf('Caso %u - Gradiente',caso));
% 
% saveas(gcf,sprintf('3D - Caso %u\\ConvergenceGL2.png',caso));



figure(1)
loglog(n,ErroL2(:,caso,3),'-*r',n,ErroL2(:,caso,9),'-ob',refx,refy,'-k','LineWidth',3);
h1 = legend('MPFA-DNL','MPFA-D');
xlabel('\it\fontname{Times} n');
ylabel('\it\fontname{Times} l^{2}_{u}');
% title('Scalar Variable')
set(h1,'FontSize',20);
set(gca,'FontSize',20)

saveas(gcf,sprintf('3D - Caso %u\\ConvergenceEL2.png',caso));

figure(2)
loglog(n,GradL2(:,caso,3),'-*r',n,GradL2(:,caso,9),'-ob',rgfx,rgfy,'-k','LineWidth',3);
h2 = legend('MPFA-DNL','MPFA-D');
xlabel('\it\fontname{Times} n');
ylabel('\it\fontname{Times} l^{2}_{\nablau}');
% title('Variable Gradient')
ylim([0.05 1.5])
set(h2,'FontSize',20);
set(gca,'FontSize',20)

saveas(gcf,sprintf('3D - Caso %u\\ConvergenceGL2.png',caso));
