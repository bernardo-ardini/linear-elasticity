clear;
close all;

quale="mesh0";
stile="none";

el=elliptic("../data",quale);
subplot(1,2,1);
el.solve("IC");
el.plot(stile);

par=parabolic("../data",quale);

par.dt=0.02;
par.T=10;
par.zeta=@(t) (t/5)*(t<=5)+(t>5);
par.dzeta=@(t) 1/5*(t<=5);
par.u0=zeros(size(par.coord,1),1);
w=par.solve("IC");

subplot(1,2,2);
par.plot(100,stile);
par.animate(stile);
%par.movie(stile,"movie");
