function p = slopelim( wells, Keq, GV )
%
global vertex sist element

p = sist.Mtpfa\sist.Btpfa; pref = p; iterlim = 1000; iterGVlim = 20;
crit = 1e18; CPAR = (max(p)-min(p))*1e-3; iter = 1; res = 0; iterres = 0;
a = zeros(size(vertex.coord,1),1); % ENQUANTO NÃO TEM FACES DE NEUMANN
R = []; Q = []; G = []; verimod = 1; verispec = 1;
AAcrit = (max(p)-min(p))*1e-1; AAstart = 0; mAA = 0; f_old = p; g_old = p;

while crit > CPAR && iter < iterlim
    
    [ minpvec, maxpvec ] = determinmax( p ); critGV = 0; iterGV = 1;
    
    if iter==1, minref = min(minpvec); maxref = max(maxpvec); end
    
    while critGV == 0 && iterGV < iterGVlim
        
        if verimod==1    
            [ GV, verimod ] = modifysitem( p, Keq, GV, minpvec, maxpvec );
        end

        M = sist.Mtpfa + sist.Mcdt;
        B = sist.Btpfa + sist.Bcdt;

        [ M, B ] = incluiwells( M, B, wells );

        p = p + (B - M*p)./diag(M);

        r = p - pref;
        
        iterGV = iterGV + 1;
        
        if min(p) > (minref - 1e-9) && max(p) < (maxref + 1e-9)
            critGV = 1; verimod = 0;
            p(p < minref & p > (minref - 1e-9)) = minref;
            p(p > maxref & p < (maxref + 1e-9)) = maxref;
        elseif iterGV < iterGVlim
            p = pref; verimod = 1;
        end
        
    end
    
    % Verifico se o sistema da iteração k contribui para a convergência --%
    if verispec==1 || verimod==1
        W = eye(size(M))-inv(diag(diag(M)))*M;
        specradius = abs(eigs(W, 1, 'lm'));
        if specradius>1
            verispec = 1; verimod = 1; convjac = 0;
            while convjac==0
                deltaro = 0.0840*specradius + 0.0359;
                if (1/deltaro)*((specradius - 0.99)/specradius) > 1
                    beta = 0.75^specradius; 
                else
                    beta = sqrt(1 - (1/deltaro)*((specradius - 0.99)/specradius)); 
                end
                sist.Mcdt = beta*sist.Mcdt; sist.Bcdt = beta*sist.Bcdt; GV = beta*GV;
                M = sist.Mtpfa + sist.Mcdt; B = sist.Btpfa + sist.Bcdt;
                W = eye(size(M))-inv(diag(diag(M)))*M;
                specradius = abs(eigs(W, 1, 'lm'));
                convjac = specradius<1;
            end
            p = pref; verispec = 0;
            p = p + (B - M*p)./diag(M);
        else
            verispec = 0; 
        end
    end
    %---------------------------------------------------------------------%
    
    %---------------------------------------------------------------------%
    crit = sqrt(sum((abs(r).^2).*element.volume)./sum((pref.^2).*element.volume));
    
%     if crit < CPAR*1e2
%         [ p, f_old, g_old, mAA, G, R, Q ] = AndersonAcelerator( mAA, p, r, iter, ...
%                                             AAstart, G, f_old, g_old, 1, 4, 1e10, ...
%                                             R, Q, minpvec, maxpvec );   
%     end
    
    if AAstart ~= 1 && crit < AAcrit
    	AAstart = 1;
    end
    %---------------------------------------------------------------------%
    
    %---------------------------------------------------------------------%
    res(iter) = crit;
    iterres(iter) = iter;
    
    if iter == 1
        CPAR = crit*1e-3;
    end
    
    pref = p; 
    iter = iter + 1;
    %---------------------------------------------------------------------%
    
end

savegraph( iterres, res )
    
end

%-------------------------------------------------------------------------%
function [ M, B ] = vertexcontrib( M, B, Keq, GV, i )
%
global vertex options face

if vertex.flag(i)>200
            
    [ M, B ] = contribuvertex3D( M, B, Keq, GV, vertex.weights, vertex.flowcontribs, i );

else

    reg4 = vertex.flag(i)-100; 
    pV = options.solanalit{reg4}(vertex.coord(i,1),vertex.coord(i,2),vertex.coord(i,3));
    for j=1:3
        faces = find(face.inner.vertices(:,j)==i);
        if isempty(faces)==0              
            leftofaces = face.inner.montelem(faces);
            rigtofaces = face.inner.juselem(faces);
            for t=1:size(leftofaces,1)
                B(leftofaces(t)) = B(leftofaces(t)) + Keq(faces(t))*GV(faces(t),j)*pV;
                B(rigtofaces(t)) = B(rigtofaces(t)) - Keq(faces(t))*GV(faces(t),j)*pV;
            end
        end
    end

end

end

%-------------------------------------------------------------------------%
function [ M, B ] = incluiwells( M, B, wells )
%

    for iw=1:size(wells,1)      
        if (wells(iw,1)~=0)&&(wells(iw,5)>400)&&(wells(iw,5)<600) % Caso seja prescrita alguma pressão.
            M(wells(iw,1),:) = 0*M(wells(iw,1),:);
            M(wells(iw,1),wells(iw,1)) = 1;
            B(wells(iw,1)) = wells(iw,6);
        end
    end

end

%-------------------------------------------------------------------------%
function [ GV, verimod ] = modifysitem( p, Keq, GV, minpvec, maxpvec )
%
global sist element face

iter = 1; res = 0; iterres = 0; crit = 0;
    
while crit < 1 && iter < 10
    
    GVnew = GV;

    alpha = determalpha( p, minpvec, maxpvec );
    
    elemmuda = find(alpha<1);
    
    if isempty(elemmuda)==0 
        
       verimod = 1;
        
       facecorr = unique(sort([element.faces(elemmuda,1); element.faces(elemmuda,2); ...
                               element.faces(elemmuda,3); element.faces(elemmuda,4)]));  
    
        for i = facecorr'
            if i>size(face.inner.vertices,1)
                jusalpha = 1;
                montalpha = alpha(face.bound.montelem(i-size(face.inner.vertices,1)));
            else
                jusalpha = alpha(face.inner.juselem(i));
                montalpha = alpha(face.inner.montelem(i));
            end
            GVnew(i,:) = min(jusalpha,montalpha)*GV(i,:);
        end

        vertcorr = unique(sort([element.vertices(elemmuda,1); element.vertices(elemmuda,2); ...
                                element.vertices(elemmuda,3); element.vertices(elemmuda,4)]));

        for j = vertcorr'       
            [ sist.Mcdt, sist.Bcdt ] = vertexcontrib( sist.Mcdt, sist.Bcdt, Keq, -GV, j );
            [ sist.Mcdt, sist.Bcdt ] = vertexcontrib( sist.Mcdt, sist.Bcdt, Keq, GVnew, j );
        end

    else
        
        verimod = 0;
    
    end
    
    GV = GVnew;
    
    crit = round(min(alpha),4); 
    res(iter) = crit; 
    iterres(iter) = iter;
    iter = iter + 1;
    
end

end

%--------------------------------------------------------------------------
function savegraph( iterres, res )
%
global options

    if strcmp(options.limitadora,'sim')
        typewei = [options.tipopeso '_' 'Limited'];
    else
        typewei = options.tipopeso;
    end

    [~,nameso,~] = fileparts(options.malha);
    
    figure(1); loglog( iterres, res );
    xlabel('Iterations'); ylabel('Residue'); grid on;    
    
    if exist(sprintf('3D - Caso %u\\%s\\%s',options.caso,typewei,nameso),'dir')==0
        mkdir(sprintf('3D - Caso %u\\%s\\%s',options.caso,typewei,nameso));
    end    
    saveas(gcf,sprintf('3D - Caso %u\\%s\\%s\\Resid_vs_Iter.png',options.caso,typewei,nameso))

end

%--------------------------------------------------------------------------
function [ minpvec, maxpvec ] = determinmax( p )
%
global vertex element

maxpvec = []; minpvec = []; pno = zeros(size(p,1),2);
% DETERMINAÇÃO DAS PRESSÕES MÁXIMA E MÍNIMA DO ESTÊNCIL ESTENDIDO ---------
for i=1:size(p,1)
    pvec = []; ivec = [];
    for j=1:size(element.vertices,2)
        pno = recuperapressnos3D( p, pno, element.vertices(i,j) );
        inicesurn = vertex.elsurvertpointer(element.vertices(i,j));
        if element.vertices(i,j)==size(vertex.elsurvertpointer,2)
            fimesurn = size(vertex.elsurvertex,2);
        else
            fimesurn = vertex.elsurvertpointer(element.vertices(i,j)+1)-1;
        end
        pvec = [pvec; p(vertex.elsurvertex(inicesurn:fimesurn))];
        ivec = [ivec; vertex.elsurvertex(inicesurn:fimesurn)'];
    end

    vecpno = pno(element.vertices(i,1:size(element.vertices,2)),1);
    pospno = find(vertex.flag(element.vertices(i,1:size(element.vertices,2)))<200);
    if isempty(pospno)==0
        vecaccno = vecpno(pospno);
        pvec = [pvec; vecaccno];
    end

    maxpvec(i,1) = max(pvec) + element.sourceterm(i)/element.volume(i); 
    minpvec(i,1) = min(pvec) + element.sourceterm(i)/element.volume(i);
    
end
%--------------------------------------------------------------------------

end

%-------------------------------------------------------------------------%
function [ pno ] = recuperapressnos3D( p, pno, i )
% SE FOR UM VÉRTICE NUMA FACE DE NEUMANN DEVE DAR ERRADO!!! --------------
global options vertex

if pno(i,2)==0
    if vertex.flag(i)<200
        for j=1:4
            reg4 = vertex.flag(i)-100; 
        end
        I = vertex.coord(i,:)';
        pno(i,1) = options.solanalit{reg4}(I(1),I(2),I(3)); pno(i,2) = 1;
    else
        inicesurn = vertex.elsurvertpointer(i);
        if i==size(vertex.elsurvertpointer,2)
            fimesurn = size(vertex.elsurvertex,2);
        else
            fimesurn = vertex.elsurvertpointer(i+1)-1;
        end
        pre = p(vertex.elsurvertex(inicesurn:fimesurn));
        wei = vertex.weights(inicesurn:fimesurn);
        pno(i,1) = wei*pre; pno(i,2) = 1;
    end
end

end



