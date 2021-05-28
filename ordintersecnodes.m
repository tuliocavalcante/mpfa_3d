function [ polyf, polyt ] = ordintersecnodes( polyf, polyt )
%

if isempty(polyt)==0
    centpolyt = mean(polyt);
    V(1,:) = polyt(1,:) - centpolyt; thetat(1)=0;
    for i=2:size(polyt,1)
        V(i,:) = polyt(i,:) - centpolyt;
        c = cross(V(1,:),V(i,:));
        if c(3)<0
            thetat(i) = 2*pi - acos((V(1,:)*V(i,:)')/(norm(V(1,:))*norm(V(i,:))));
        elseif c(3)>=0
            thetat(i) = acos((V(1,:)*V(i,:)')/(norm(V(1,:))*norm(V(i,:))));
        end   
    end
    [~,idxt] = sort(thetat);
    polyt = polyt(idxt,:);
end

if isempty(polyf)==0
    centpolyf = mean(polyf);
    U(1,:) = polyf(1,:) - centpolyf; thetaf(1)=0;
    for i=2:size(polyf,1)
        U(i,:) = polyf(i,:) - centpolyf;
        c = cross(U(1,:),U(i,:));
        if c(3)<0
            thetaf(i) = 2*pi - acos((U(1,:)*U(i,:)')/(norm(U(1,:))*norm(U(i,:))));
        elseif c(3)>=0
            thetaf(i) = acos((U(1,:)*U(i,:)')/(norm(U(1,:))*norm(U(i,:))));
        end   
    end
    [~,idxf] = sort(thetaf);
    polyf = polyf(idxf,:);
end

end

