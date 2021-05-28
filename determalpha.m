function alpha = determalpha( p, minpvec, maxpvec )
%
global sist
% CÁLCULO DO FATOR QUE RESTRINGIRÁ A CONTRIBUIÇÃO DO TERMO DE
% DIFUSÃO CRUZADO --------------------------------------------------------
alphamax = ones(size(p)); alphamin = ones(size(p));
dm = 0.01*( maxpvec - minpvec ); 
vecpmax = [ maxpvec-dm maxpvec maxpvec+dm ]; 
vecpmin = [ minpvec-dm minpvec minpvec+dm ];

num1 = diag(sist.Mtpfa).*p + sist.Btpfa - sist.Mtpfa*p;
den1 = - diag(sist.Mcdt).*p - sist.Bcdt + sist.Mcdt*p + 1e-16;

nummax = [num1 num1 num1] - [diag(sist.Mtpfa) diag(sist.Mtpfa) diag(sist.Mtpfa)].*vecpmax;
nummin = [num1 num1 num1] - [diag(sist.Mtpfa) diag(sist.Mtpfa) diag(sist.Mtpfa)].*vecpmin;

denmax = [diag(sist.Mcdt) diag(sist.Mcdt) diag(sist.Mcdt)].*vecpmax + [den1 den1 den1];
denmin = [diag(sist.Mcdt) diag(sist.Mcdt) diag(sist.Mcdt)].*vecpmin + [den1 den1 den1];

vecfactor1 = nummax./denmax; vecfactor2 = nummin./denmin;

facmax = vecfactor1(:,2); facmin = vecfactor2(:,2);
dpmax_dfacmax = (vecpmax(:,3)-vecpmax(:,1))./(vecfactor1(:,3)-vecfactor1(:,1));
dpmin_dfacmin = (vecpmin(:,3)-vecpmin(:,1))./(vecfactor2(:,3)-vecfactor2(:,1));

posmax = find((facmax>=0).*(facmax<=1).*(dpmax_dfacmax>0)~=0);
posmin = find((facmin>=0).*(facmin<=1).*(dpmin_dfacmin<0)~=0);

if isempty(posmax)==0, alphamax(posmax) = facmax(posmax); end
if isempty(posmin)==0, alphamin(posmin) = facmin(posmin); end

alpha = min(alphamin,alphamax);   

end
