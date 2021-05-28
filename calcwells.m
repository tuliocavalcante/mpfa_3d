function wells = calcwells( element, vertex )
%

dirnodes = find(vertex.flag<200);

for i=1:size(dirnodes,1)
    no = dirnodes(i);
    col_1 = find(element.vertices(:,1)==no);
    col_2 = find(element.vertices(:,2)==no);
    col_3 = find(element.vertices(:,3)==no);
    col_4 = find(element.vertices(:,4)==no);
    if vertex.flag(i)==101
        inj = unique([col_1;col_2;col_3;col_4]);
    else
        prod = unique([col_1;col_2;col_3;col_4]);
    end
end

wells = zeros(size(inj,1)+size(prod,1),6);

for iw=1:size(inj,1)
    wells(iw,1) = 1;
    wells(iw,2) = inj(iw);
    wells(iw,5) = 400 + iw;
    wells(iw,6) = 1; % AJEITAR ISSO!!!
end

for iw=1:size(prod,1)
    wells(iw + size(inj,1),1) = 2;
    wells(iw + size(inj,1),2) = prod(iw);
    wells(iw + size(inj,1),5) = 400 + iw + size(inj,1);
    wells(iw + size(inj,1),6) = 0; % AJEITAR ISSO!!!
end


end

