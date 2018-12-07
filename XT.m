clear
nx=101;ny=101;[xp,tp]=MakeRectMesh(nx,ny);
[xu,tu,tuplot,bnd]=Triangle3to6(xp,tp);
%% locate boundaries
minx=min(xu);maxx=max(xu);smallno=1e-8;% any smaller makes trouble
S=find(xu(:,2)<smallno+minx(2));
N=find(xu(:,2)>maxx(2)-smallno);
W=find(xu(:,1)<smallno+minx(1));
E=find(xu(:,1)>maxx(1)-smallno);
%% for Cooscillate remove duplicate points from N and S
%[~,order]=sort(xu(N,1));N=N(order);N=N(2:end-1);
%[~,order]=sort(xu(S,1));S=S(order);S=S(2:end-1);
%% for SV remove duplicate points from W and E (means changing ReadXT.f90)
[~,order]=sort(xu(W,2));W=W(order);W=W(2:end-1);
[~,order]=sort(xu(E,2));E=E(order);E=E(2:end-1);
[~,order]=sort(xu(N,1));N=N(order);
[~,order]=sort(xu(S,1));S=S(order);
save XT.mat xp tp xu tu tuplot N W S E
