clear;
## usage: TitleFrame (str,array)
##
##
function TitleFrame (str,array)
  maxa=max(array);mina=min(array);
  titlestring=strcat(str,{', range: '},num2str(mina,"%3.1f"),{' to '},num2str(maxa,"%3.1f"));
  title(titlestring)
endfunction

load XT
nx=length(xp(:,1));nu=length(xu(:,1))
load Soln
u=Soln(1:nu);v=Soln(nu+1:2*nu);p=Soln(2*nu+1:2*nu+nx);
%omega=M\(Gx*v'-Gy*u');psi=K\(M*omega);psi=psi-mean(psi);

set(gcf,"position",[50,700,1200,400]);
%% Map u,omega,p

clf;axes("position",[0.03,0.05,0.3,0.8]);
mapfem(xu,tuplot,u);colormap(jet);
TitleFrame ('u',u)

axes("position",[0.35,0.05,0.3,0.8]);%
mapfem(xu,tuplot,v);
TitleFrame ('v',v)

axes("position",[0.68,0.05,0.3,0.8]);
mapfem(xp,tp,p)
TitleFrame ('Pressure',p)
