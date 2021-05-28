function [ newnodes ] = intersecnewnodes( polyt, polyf )

u = 1;
t = [1:size(polyt,1) 1];
newnodes = [];
for i=1:size(polyt,1)
    
    f = [1:size(polyf,1) 1];
    for j=1:size(polyf,1)
        
       M = [polyt(t(i+1),1:2)'-polyt(t(i),1:2)' polyf(f(j),1:2)'-polyf(f(j+1),1:2)'];
       B = [polyf(f(j),1:2)'-polyt(t(i),1:2)'];
       
       tf = M\B;
       
       if sum(isinf(tf))==0
           intt = round(polyt(t(i),:) + (polyt(t(i+1),:)-polyt(t(i),:))*tf(1),12);
%            intf = round(polyf(f(j),:) + (polyf(f(j+1),:)-polyf(f(j),:))*tf(2),12);
%            if sum(intt==intf)<3
%                ddd=81;
%            end
           [ ~, fint ] = insideintnodes( polyt, intt );
           [ tinf, ~ ] = insideintnodes( intt, polyf );
           if isempty(fint)==0 && isempty(tinf)==0
               newnodes(u,:) = fint;
               u = u + 1;
           end
       end
        
    end
    
end

end

