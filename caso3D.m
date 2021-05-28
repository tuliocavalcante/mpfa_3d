function [ f, K, wells ] = caso3D( caso, element, vertex )
%

wells = zeros(1,8);

% LINEAR
if caso==1
    syms x y z
    f = {@(x,y,z) x};
    K = {@(x,y,z) [1 0 0;0 1 0;0 0 1]};

% QUADRÁTICO
elseif caso==2
    syms x y z
    f = {@(x,y,z) y^2};
    K = {@(x,y,z) [2 1 0;1 2 1;0 1 2]};

% PLANO
elseif caso==3
    syms x y z
    f = {@(x,y,z) 2*x+3*y+z};
    K = {@(x,y,z) [100 1 0;1 100 1;0 1 100]};

% CÚBICO
elseif caso==4
    syms x y z
    f = {@(x,y,z) z^3};
    K = {@(x,y,z) [1 0 0;0 1 0;0 0 1]};
    
% ALTAMENTE HETEROGÊNEO
elseif caso==5
    syms x y z
    f = {@(x,y,z) 2*x+2*y};
    K = {@(x,y,z) [x+1 0 0;0 y+1 0;0 0 z+1]};
   
% TESTE 1 DO BENCHMARK DE EYMARD ET AL. (2011) 
elseif caso==6
    syms x y z
    f = {@(x,y,z) 1+sin(pi*x)*sin(pi*(y+(1/2)))*sin(pi*(z+(1/3)))};
    K = {@(x,y,z) [1 0.5 0;0.5 1 0.5;0 0.5 1]};

% TESTE 2 DO BENCHMARK DE EYMARD ET AL. (2011)
elseif caso==7
    syms x y z
    f = {@(x,y,z) x^3*y^2*z+x*sin(2*pi*x*z)*sin(2*pi*x*y)*sin(2*pi*z)};
    K = {@(x,y,z) [y^2+z^2+1 -x*y -x*z;-x*y x^2+z^2+1 -y*z;-x*z -y*z x^2+y^2+1]};
    
% TESTE 3 DO BENCHMARK DE EYMARD ET AL. (2011)
elseif caso==8
    syms x y z
    f = {@(x,y,z) sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z)};
    K = {@(x,y,z) [1 0 0;0 1 0;0 0 10^3]};

% TESTE 5 DO BENCHMARK DE EYMARD ET AL. (2011)
elseif caso==9
    syms x y z
    f = {@(x,y,z) 0.1*sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z),...
        @(x,y,z) 10*sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z),...
        @(x,y,z) 100*sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z),...
        @(x,y,z) 0.01*sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z)};   
    K = {@(x,y,z) [1 0 0; 0 10 0; 0 0 0.01],...
        @(x,y,z) [1 0 0; 0 0.1 0; 0 0 100],...
        @(x,y,z) [1 0 0; 0 0.01 0; 0 0 10],...
        @(x,y,z) [1 0 0; 0 100 0; 0 0 0.1]};
    
% CASCA ALTAMENTE ANISOTRÓPICA - VIOLA DMP
elseif caso==10
    syms x y z
    Rx = [1 0 0; 0 cos(pi/3) -sin(pi/3); 0 sin(pi/3) cos(pi/3)];
    Ry = [cos(pi/4) 0 sin(pi/4); 0 1 0; -sin(pi/4) 0 cos(pi/4)];
    Rz = [cos(pi/6) -sin(pi/6) 0; sin(pi/6) cos(pi/6) 0; 0 0 1];
    A = [100 0 0; 0 10 0; 0 0 1];
    f = {@(x,y,z) 0, @(x,y,z) 2};
    K = {@(x,y,z) Rz'*Ry'*Rx'*A*Rx*Ry*Rz};
                    
% OBLIQUE FLOW
elseif caso==11
    syms x y z
    f = {@(x,y,z) 0.5, @(x,y,z) 1, @(x,y,z) 0};
    theta = 40*pi/180;
    A = [1 0 0; 0 0.001 0; 0 0 1];
    R = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
    K = {@(x,y,z) R*A*R',...
         @(x,y,z) R*A*R',... 
         @(x,y,z) R*A*R'};
     
% OBLIQUE DRAIN
elseif caso==12
    syms x y z
    f = {@(x,y,z) -x-0.2*y, @(x,y,z) -x-0.2*y, @(x,y,z) -x-0.2*y};
    R = [cos(atan(0.2)) -sin(atan(0.2)) 0; sin(atan(0.2)) cos(atan(0.2)) 0; 0 0 1];
    A = [1 0 0; 0 0.1 0; 0 0 1];
    B = [100 0 0; 0 10 0; 0 0 1];
    K = {@(x,y,z) R*A*R',...
         @(x,y,z) R*B*R',... 
         @(x,y,z) R*A*R'};
     
% OBLIQUE BARRIER - SOLUÇÃO CONTÍNUA
elseif caso==13
    syms x y z
    f = {@(x,y,z) -(y-0.2*(x-0.5)-0.475), @(x,y,z) -(y-0.2*(x-0.5)-0.475)/0.01, @(x,y,z) -(y-0.2*(x-0.5)+0.02)};
    K = {@(x,y,z) [1 0 0; 0 1 0; 0 0 1],...
         @(x,y,z) [0.01 0 0; 0 0.01 0; 0 0 1],... 
         @(x,y,z) [1 0 0; 0 1 0; 0 0 1]};
     
% OBLIQUE BARRIER - SOLUÇÃO DESCONTÍNUA
elseif caso==14
    syms x y z
    f = {@(x,y,z) -(y-0.2*(x-0.5)-0.475), @(x,y,z) -(y-0.2*(x-0.5)-0.475)/0.01, @(x,y,z) -(y-0.2*(x-0.5)-0.475-0.05)-0.05/0.01};
    K = {@(x,y,z) [1 0 0; 0 1 0; 0 0 1],...
         @(x,y,z) [0.01 0 0; 0 0.01 0; 0 0 1],... 
         @(x,y,z) [1 0 0; 0 1 0; 0 0 1]};
     
% SOLUÇÃO PLANAR POR PARTES
elseif caso==15
    syms x y z
    f = {@(x,y,z) x-2*z+1.5, @(x,y,z) x+0.5};
    K = {@(x,y,z) [5 0 2; 0 1 0; 2 0 1],...
         @(x,y,z) [1 0 0; 0 1 0; 0 0 1]};
     
% TESTE 6 DE YANG AND GAO (2019) 
elseif caso==16
    syms x y z
    f = {@(x,y,z) 1+sin(pi*x)*sin(pi*(y+(1/2)))*sin(pi*(z+(1/3)))};
    K = {@(x,y,z) [1.5 0.5 0;0.5 1.5 0.5;0 0.5 1]};
    
% TESTE SOUZA-QUEIROZ VERSÃO 3D - FALTA AJUSTAR!!!
elseif caso==17
    syms x y z
    Rx = [1 0 0; 0 cos(pi/3) -sin(pi/3); 0 sin(pi/3) cos(pi/3)];
    Ry = [cos(pi/4) 0 sin(pi/4); 0 1 0; -sin(pi/4) 0 cos(pi/4)];
    Rz = [cos(pi/6) -sin(pi/6) 0; sin(pi/6) cos(pi/6) 0; 0 0 1];
    A = [100 0 0; 0 10 0; 0 0 1]; csi = 1000;
    f = {@(x,y,z) 0, @(x,y,z) 2};
    K = {@(x,y,z) Rz'*Ry'*Rx'*A*Rx*Ry*Rz, ...
         @(x,y,z) [csi*(x+csi^-1)^2+(y+csi^-1)^2+(z+csi^-1)^2 -(1-csi)*(x+csi^-1)*(y+csi^-1) -(1-csi)*(x+csi^-1)*(z+csi^-1);...
                   -(1-csi)*(x+csi^-1)*(y+csi^-1) (x+csi^-1)^2+csi*(y+csi^-1)^2+(z+csi^-1)^2 -(1-csi)*(y+csi^-1)*(z+csi^-1);...
                   -(1-csi)*(x+csi^-1)*(z+csi^-1) -(1-csi)*(y+csi^-1)*(z+csi^-1) (x+csi^-1)^2+(y+csi^-1)^2+csi*(z+csi^-1)^2]};
    
% TESTE 3 DO BENCHMARK DE EYMARD ET AL. (2011) ALTERADO POR ZHANG E ALKOBASI (2019)
elseif caso==18
    syms x y z
    f = {@(x,y,z) sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z) + 1};
    K = {@(x,y,z) [1 0 0;0 1 0;0 0 10^3]};
    
% YANG AND GAO (2019) - LP VERIFICATION
elseif caso==19
    syms x y z
    f = {@(x,y,z) 3*x + 4*y + 6*z};
    K = {@(x,y,z) [1 0 0;0 1 0;0 0 1]};

% FIVE SPOT 3D
elseif caso==20
    syms x y z
    f = {@(x,y,z) x+y+z, @(x,y,z) x+y+z};
    K = {@(x,y,z) [1 0 0;0 1 0;0 0 1]};
    wells = calcwells( element, vertex );
    
end   

end

