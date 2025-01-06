clear;
close all;

el=elasticity("../data","p01");

function r=R3(theta)
    r=[cosd(theta),sind(theta),0;-sind(theta),cosd(theta),0;0,0,1];    
end
function r=R1(theta)
    r=[1,0,0;0,cosd(theta),sind(theta);0,-sind(theta),cosd(theta)];    
end
function r=R2(theta)
    r=[cosd(theta),0,sind(theta);0,1,0;-sind(theta),0,cosd(theta)];    
end

bound=[];

for a=1:size(el.coord,1)
    if abs(el.coord(a,1))<1e-3
        bound=[bound;a,-el.coord(a,:)+el.coord(a,:)*R3(40)'+[-1,0.5,0.3]];
    elseif abs(el.coord(a,1)-3)<1e-3
        bound=[bound;a,0,0,0];
    end
end

writematrix(bound,'../data/p01/p01.dat','Delimiter','tab');

el=elasticity("../data","p01");
clc;

el.lambda=0;
el.mu=0.3;
[u,~,~]=el.solve("J");
figure;
el.plot3();
