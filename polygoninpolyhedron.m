function [ element, fracture, face ] = polygoninpolyhedron ( element, fracture, face, vertex, options )
% Esse algoritmo está limitado para tetraedros, desde a determinação dos
% vértices I, J e K, e o uso deles para determinar as 4 normais. Essas
% regras estão limitadas a tetraedros. Consequentemente os 6 sistemas de
% equações que surgem das combinações (2 a 2) dos 4 planos que contêm as 
% faces dos tetraedros também estão dentro dessa limitação. Para
% generalizar esse algoritmo para poliedros, precisaríamos de uma função
% que retornasse as normais de maneira geral, fossem quantas fossem, além
% de guardar um vértice de cada plano (de cada normal). Depois 
% precisaríamos fazer a combinação dos planos 2 a 2 (fossem quantos 
% fossem). Depois de determinados os pontos de intersecção entre os planos,
% o resto do código pode ser aproveitado, porque as funções que ordenam os
% vértices, organizam os polígonos e calculam a intersecção entre eles
% podem lidar com quaisquer números de vértices.

fracmatcon = zeros(size(fracture.coord,1),size(element.vertices,1));
fracfraccon = zeros(size(fracture.coord,1),size(fracture.coord,1));
fracconfra = zeros(size(fracture.coord,1),size(fracture.coord,1),3);

montelem = [face.inner.montelem; face.bound.montelem];
E1 = -2*(montelem(element.faces(:,1),:) ~= [1:size(element.vertices,1)]') + 1;
E2 = -2*(montelem(element.faces(:,2),:) ~= [1:size(element.vertices,1)]') + 1;
E3 = -2*(montelem(element.faces(:,3),:) ~= [1:size(element.vertices,1)]') + 1;
E4 = -2*(montelem(element.faces(:,4),:) ~= [1:size(element.vertices,1)]') + 1;

N = [face.inner.normal; face.bound.normal];
N1 = N(element.faces(:,1),:).*[E1 E1 E1];
N2 = N(element.faces(:,2),:).*[E2 E2 E2];
N3 = N(element.faces(:,3),:).*[E3 E3 E3]; 
N4 = N(element.faces(:,4),:).*[E4 E4 E4];

F = [face.inner.vertices; face.bound.vertices];
F1 = F(element.faces(:,1),:); F2 = F(element.faces(:,2),:);
F3 = F(element.faces(:,3),:); F4 = F(element.faces(:,4),:);


% I = coord(elem(:,1),:); J = coord(elem(:,2),:); K = coord(elem(:,3),:); Q = coord(elem(:,4),:);
% JI = I - J; JK = K - J; JQ = Q - J; IQ = Q - I; IK = K - I; CJ = J - centelem; CI = I - centelem;
% N1 = 0.5*cross(JI,JK).*[sign(diag(cross(JI,JK)*CJ')) sign(diag(cross(JI,JK)*CJ')) sign(diag(cross(JI,JK)*CJ'))];
% N2 = 0.5*cross(JK,JQ).*[sign(diag(cross(JK,JQ)*CJ')) sign(diag(cross(JK,JQ)*CJ')) sign(diag(cross(JK,JQ)*CJ'))]; 
% N3 = 0.5*cross(JQ,JI).*[sign(diag(cross(JQ,JI)*CJ')) sign(diag(cross(JQ,JI)*CJ')) sign(diag(cross(JQ,JI)*CJ'))]; 
% N4 = 0.5*cross(IQ,IK).*[sign(diag(cross(IQ,IK)*CI')) sign(diag(cross(IQ,IK)*CI')) sign(diag(cross(IQ,IK)*CI'))];

% Os planos estão definidos como: ----------------------------------------%
% pi_1: N1(:,1)*x + N1(:,2)*y + N1(:,3)*z = ( N1(:,1)*J(:,1) + N1(:,2)*J(:,2) + N1(:,3)*J(:,3) ) 
% pi_2: N2(:,1)*x + N2(:,2)*y + N2(:,3)*z = ( N2(:,1)*J(:,1) + N2(:,2)*J(:,2) + N2(:,3)*J(:,3) )
% pi_3: N3(:,1)*x + N3(:,2)*y + N3(:,3)*z = ( N3(:,1)*J(:,1) + N3(:,2)*J(:,2) + N3(:,3)*J(:,3) )
% pi_4: N4(:,1)*x + N4(:,2)*y + N4(:,3)*z = ( N4(:,1)*I(:,1) + N4(:,2)*I(:,2) + N4(:,3)*I(:,3) ) 
% pi_f: Nf(:,1)*x + Nf(:,2)*y + Nf(:,3)*z = ( Nf(:,1)*Jf(:,1) + Nf(:,2)*Jf(:,2) + Nf(:,3)*Jf(:,3) )
%-------------------------------------------------------------------------%

for frati=1:size(fracture.coord,1)
    
    If(frati,:) = fracture.coord{frati}(1,:); Jf(frati,:) = fracture.coord{frati}(2,:); 
    Kf(frati,:) = fracture.coord{frati}(3,:);

    JfIf(frati,:) = If(frati,:) - Jf(frati,:); JfKf(frati,:) = Kf(frati,:) - Jf(frati,:); 
    Nf(frati,:) = 0.5*cross(JfIf(frati,:),JfKf(frati,:));
    
    for tetraj=1:size(element.vertices,1)
        
        I(tetraj,:) = vertex.coord(setdiff(element.vertices(tetraj,:),F2(tetraj,:)),:); 
        J(tetraj,:) = vertex.coord(setdiff(element.vertices(tetraj,:),F4(tetraj,:)),:); 
        K(tetraj,:) = vertex.coord(setdiff(element.vertices(tetraj,:),F3(tetraj,:)),:); 
        Q(tetraj,:) = vertex.coord(setdiff(element.vertices(tetraj,:),F1(tetraj,:)),:);
                       
        M1 = [N1(tetraj,:); N2(tetraj,:); Nf(frati,:)];
        B1 = [N1(tetraj,:)*J(tetraj,:)'; N2(tetraj,:)*J(tetraj,:)'; Nf(frati,:)*Jf(frati,:)'];
        
        M2 = [N1(tetraj,:); N3(tetraj,:); Nf(frati,:)];
        B2 = [N1(tetraj,:)*J(tetraj,:)'; N3(tetraj,:)*J(tetraj,:)'; Nf(frati,:)*Jf(frati,:)'];
        
        M3 = [N2(tetraj,:); N3(tetraj,:); Nf(frati,:)];
        B3 = [N2(tetraj,:)*J(tetraj,:)'; N3(tetraj,:)*J(tetraj,:)'; Nf(frati,:)*Jf(frati,:)'];
        
        M4 = [N2(tetraj,:); N4(tetraj,:); Nf(frati,:)];
        B4 = [N2(tetraj,:)*K(tetraj,:)'; N4(tetraj,:)*K(tetraj,:)'; Nf(frati,:)*Jf(frati,:)'];
        
        M5 = [N1(tetraj,:); N4(tetraj,:); Nf(frati,:)];
        B5 = [N1(tetraj,:)*I(tetraj,:)'; N4(tetraj,:)*I(tetraj,:)'; Nf(frati,:)*Jf(frati,:)'];
        
        M6 = [N4(tetraj,:); N3(tetraj,:); Nf(frati,:)];
        B6 = [N4(tetraj,:)*I(tetraj,:)'; N3(tetraj,:)*I(tetraj,:)'; Nf(frati,:)*Jf(frati,:)'];
        
        P(1,:) = (M1\B1)';
        P(2,:) = (M2\B2)';
        P(3,:) = (M3\B3)';
        P(4,:) = (M4\B4)';
        P(5,:) = (M5\B5)';
        P(6,:) = (M6\B6)';
        
        for k=1:6
            JP = P(k,:) - J(tetraj,:); 
            IP = P(k,:) - I(tetraj,:); 
            ver1 = acos(round(JP*N1(tetraj,:)'/(norm(JP)*norm(N1(tetraj,:))),6))>=(pi/2);  
            ver2 = acos(round(JP*N2(tetraj,:)'/(norm(JP)*norm(N2(tetraj,:))),6))>=(pi/2);  
            ver3 = acos(round(JP*N3(tetraj,:)'/(norm(JP)*norm(N3(tetraj,:))),6))>=(pi/2);  
            ver4 = acos(round(IP*N4(tetraj,:)'/(norm(IP)*norm(N4(tetraj,:))),6))>=(pi/2);            
            verif(1,k) = ver1*ver2*ver3*ver4; 
        end
        posk = find(verif>0);
        
        if isempty(posk)==0
        
            polyf = unique(fracture.coord{frati},'rows');
            if isempty(posk)==0
                polyt = P(posk,:);
            end
                                    
            [ polyf, polyt ] = ordintersecnodes( polyf, polyt );
            
            [ tinf, fint ] = insideintnodes( polyt, polyf );

            [ newnodes ] = intersecnewnodes( polyt, polyf );

            intpolynodes = unique([tinf; fint; newnodes],'rows');
            
            if isempty(intpolynodes)==0
                
                [ intpolynodes, ~ ] = ordintersecnodes( intpolynodes, intpolynodes );

                [ area ] = areapoly3D( intpolynodes );
                
                [ face, element, fracture ] = portionvolumes( vertex.coord(element.vertices(tetraj,:),:), ...
                                              N1, N2, N3, N4, Nf, intpolynodes, area, tetraj, frati, face, ...
                                              element, fracture, vertex, options );
               

%                 figinternodes = [ intpolynodes; intpolynodes(1,:) ];
%                 polyf1 = zeros(size(polyf,1)+1,3);
%                 polyf1(1:size(polyf,1),:) = polyf; polyf1(size(polyf,1)+1,:) = polyf(1,:);
%                 polyt1 = zeros(size(polyt,1)+1,3);
%                 polyt1(1:size(polyt,1),:) = polyt; polyt1(size(polyt,1)+1,:) = polyt(1,:);
%                 tretraface1 = vertex.coord(element.vertices(tetraj,[1 2 3 1]),:); 
%                 tretraface2 = vertex.coord(element.vertices(tetraj,[1 2 4 1]),:);
%                 tretraface3 = vertex.coord(element.vertices(tetraj,[1 3 4 1]),:); 
%                 tretraface4 = vertex.coord(element.vertices(tetraj,[2 4 3 2]),:);
%                 plot3(polyf1(:,1),polyf1(:,2),polyf1(:,3),'--r',...
%                       polyt1(:,1),polyt1(:,2),polyt1(:,3),'-k',...
%                       tretraface1(:,1),tretraface1(:,2),tretraface1(:,3),'-b',...
%                       tretraface2(:,1),tretraface2(:,2),tretraface2(:,3),'-b',...
%                       tretraface3(:,1),tretraface3(:,2),tretraface3(:,3),'-b',...
%                       tretraface4(:,1),tretraface4(:,2),tretraface4(:,3),'-b',...
%                       figinternodes(:,1),figinternodes(:,2),figinternodes(:,3),'-om')
%                                     
%                   xlabel('x');
%                   ylabel('y');
%                   zlabel('z');
%                   
%                   grid on
%                   
%                   hold on

            else

                area = 0;

%                 polyf1 = zeros(size(polyf,1)+1,3);
%                 polyf1(1:size(polyf,1),:) = polyf; polyf1(size(polyf,1)+1,:) = polyf(1,:);
%                 tretraface1 = vertex.coord(element.vertices(j,[1 2 3 1]),:); 
%                 tretraface2 = vertex.coord(element.vertices(j,[1 2 4 1]),:);
%                 tretraface3 = vertex.coord(element.vertices(j,[1 3 4 1]),:); 
%                 tretraface4 = vertex.coord(element.vertices(j,[2 4 3 2]),:);
%                 plot3(polyf1(:,1),polyf1(:,2),polyf1(:,3),'-*r',...
%                       tretraface1(:,1),tretraface1(:,2),tretraface1(:,3),'-w',...
%                       tretraface2(:,1),tretraface2(:,2),tretraface2(:,3),'-w',...
%                       tretraface3(:,1),tretraface3(:,2),tretraface3(:,3),'-w',...
%                       tretraface4(:,1),tretraface4(:,2),tretraface4(:,3),'-w')
%                   
%                 xlabel('x');
%                 ylabel('y');
%                 zlabel('z');
%                 
%                 grid on
%                 
%                 hold on

            end
            
        else

                area = 0;
                
        end
          
        fracmatcon(frati,tetraj) = area;
                
    end
    
    for fratj=1:size(fracture.coord,1) % Verifica qual fratj cruza com frati
        
        if frati ~= fratj
        
            veripos = 0; verineg = 0;
            for y=1:size(fracture.coord{fratj},1)
                Ity = fracture.coord{fratj}(y,:); JfIty = Ity - Jf(frati,:);
                veri = round(Nf(frati,:)*JfIty',6);    
                veripos = veripos + (veri>0);
                verineg = verineg + (veri<0);
            end
            verifc = veripos*verineg;

            if verifc > 0
                
                Ity = fracture.coord{fratj}(1,:); JfIty = Ity - Jf(frati,:);
                veri = round(Nf(frati,:)*JfIty'/(norm(Nf(frati,:))*norm(JfIty)),6); 
                t = 1; seg = zeros(2,2);   
                
                yseq = [1:size(fracture.coord{fratj},1) 1];

                for y=2:size(yseq,2)
                    Ity = fracture.coord{fratj}(yseq(y),:); JfIty = Ity - Jf(frati,:);
                    veri = veri*round(Nf(frati,:)*JfIty'/(norm(Nf(frati,:))*norm(JfIty)),6);    
                    if veri<0
                        seg(t,:) = [yseq(y-1) yseq(y)]; t = t + 1;
                        veri = round(Nf(frati,:)*JfIty'/(norm(Nf(frati,:))*norm(JfIty)),6);
                    else
                        veri = round(Nf(frati,:)*JfIty'/(norm(Nf(frati,:))*norm(JfIty)),6);
                    end
                end
                
                vertdif1 = fracture.coord{fratj}(seg(1,1),:);
                vdif1 = fracture.coord{fratj}(seg(1,2),:) - fracture.coord{fratj}(seg(1,1),:);
                t1 = (Nf(frati,:)*Jf(frati,:)' - Nf(frati,:)*vertdif1')/(Nf(frati,:)*vdif1');
                P1 = vertdif1 + t1*vdif1;

                vertdif2 = fracture.coord{fratj}(seg(2,1),:);              
                vdif2 = fracture.coord{fratj}(seg(2,2),:) - fracture.coord{fratj}(seg(2,1),:);                
                t2 = (Nf(frati,:)*Jf(frati,:)' - Nf(frati,:)*vertdif2')/(Nf(frati,:)*vdif2');
                P2 = vertdif2 + t2*vdif2; 

                P = [P1; P2];

                polyfr = fracture.coord{frati};

                [ ~ , pinfr ] = insideintnodes( polyfr, P );               

%                 polyf1 = fracture.coord{fratj}; 
%                 polyf2 = fracture.coord{frati}; 
%                 plot3(polyf1(:,1),polyf1(:,2),polyf1(:,3),'-*k',...
%                       polyf2(:,1),polyf2(:,2),polyf2(:,3),'-*b',...
%                       P(:,1),P(:,2),P(:,3),'-*r')
%     
%                 xlabel('x');
%                 ylabel('y');
%                 zlabel('z');
%     
%                 grid on
%     
%                 hold on

                if size(pinfr,1)==2
                    fracconfra(frati,fratj,:) = pinfr(1,:);
                    fracconfra(fratj,frati,:) = pinfr(2,:);
                elseif size(pinfr,1)==1
                    fracconfra(frati,fratj,:) = pinfr(1,:);
                end

            elseif sum(ismember(fracture.coord{fratj}, fracture.coord{frati}, 'rows'))==2

                P = fracture.coord{fratj}(ismember(fracture.coord{fratj}, fracture.coord{frati}, 'rows') == 1,:);           
                fracconfra(frati,fratj,:) = P(1,:);  fracconfra(fratj,frati,:) = P(2,:);
                

            end
                
        end
    
    end
    
end

for frati=1:size(fracture.coord,1)
    for fratj=1:size(fracture.coord,1)
        
        centi = mean(fracture.coord{frati},1); 
        veci1 = centi - reshape(fracconfra(frati,fratj,:),[1 3]); veci1 = veci1/(norm(veci1)+eps*(norm(veci1)==0));
        veci2 = reshape(fracconfra(fratj,frati,:) - fracconfra(frati,fratj,:),[1 3]); veci2 = veci2/(norm(veci2)+eps*(norm(veci2)==0));
        theta = acos(dot(veci1,veci2)); 
        
        veci1 = centi - reshape(fracconfra(frati,fratj,:),[1 3]);
        veci2 = reshape(fracconfra(fratj,frati,:) - fracconfra(frati,fratj,:),[1 3]);
        h = norm(veci1)*sin(theta);
        
        fracfraccon(frati,fratj) = (1/h)*norm(veci2);
        
%         if fracfraccon(frati,fratj)~=0
%             clear P
%             P(1,:) = fracconfra(frati,fratj,:); P(2,:) = fracconfra(fratj,frati,:);
%             polyf1 = [fracture.coord{fratj};fracture.coord{fratj}(1,:)]; 
%             polyf2 = [fracture.coord{frati};fracture.coord{frati}(1,:)]; 
%             plot3(polyf1(:,1),polyf1(:,2),polyf1(:,3),'-k',...
%                   polyf2(:,1),polyf2(:,2),polyf2(:,3),'-b',...
%                   P(:,1),P(:,2),P(:,3),'-om')
% 
%             xlabel('x');
%             ylabel('y');
%             zlabel('z');
% 
%             grid on
% 
%             hold on
%         end
                                                                
    end
end

% for frati=1:size(fracture.coord,1)
%     
%     [ a ] = areapoly3D( fracture.coord{frati} );
%     
%     areaconf(frati,1) = a;
%     
% end
% 
% teste = [areaconf sum(fracmatcon,2)]

fracture.connecwithfract = fracfraccon;
fracture.connecwithcell = fracmatcon;

end

%-------------------------------------------------------------------------%
function [ area ] = areapoly3D( nodes )
%
area = 0;

if isempty(nodes)==0 && size(nodes,1)>2

    for i=2:size(nodes,1)
        v(i-1,:) = nodes(i,:) - nodes(1,:);    
    end

    for i=1:size(v,1)-1
        area = area + 0.5*norm(cross(v(i,:),v(i+1,:)));
    end
    
end

end

