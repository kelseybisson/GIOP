
cd /Users/bissonk/Documents/MATLAB/GIOP/GIOP_July2014
load  /Users/bissonk/Documents/MATLAB/GIOP/GIOP_July2014/longhurst.mat

l = reshape(ids, [64800 1]); 
bad= find(isnan(l)==1);
l(bad)=[];
flat = inptsm(:,3);
flon = inptsm(:,4);

lat = 89.5:-1:-89.5;
lon = -179.5:1:179.5;

[long,latg] = meshgrid(lon,lat);
lidx = [reshape(latg,[64800 1]),reshape(long,[64800 1])];
lidx(bad,:)=[];

[idx]=knnsearch([lidx],[flat flon]);
check = [flat flon lidx(idx,:) l(idx)];
longidx = l(idx);

% From Rrs, correct for raman scattering of water 
a  = [0.003 0.004 0.011 0.015 0.017 0.018];
B1 = [0.014 0.015 0.010 0.010 0.010 0.010];
B2 = [ -0.022 -0.023 -0.051 -0.070 -0.080 -0.081];

% The lee et al 2013 algorithm only 
% (https://agupubs.onlinelibrary.wiley.com/doi/10.1002/jgrc.20308)
wl = [412 443 488 531 555 667];

% Only use WV - 412nm, 443nm, 488nm, 531nm, 551nm, 667nm
RF = a.*(Rrs443./Rrs555) + B1.*(Rrs555).^B2;

Rrs_t = [Rrs412 Rrs443 Rrs488 Rrs531 Rrs555 Rrs667];
Rrs_r = Rrs_t ./ (1 + RF);
gopt.aw  = get_aw(wl);
gopt.bbw = get_bbw(wl);


eta1 = 0.5:0.1:2;

for i = 1:length(Rrs_r)
    try
        for j = 1:length(eta1)
            gopt.eta = eta1(j);
[x,mapg,maph,madg,mbbp,mrrs] = giop(wl,Rrs_r(i,:),cm(i),gopt);
bbp_giop(i,:,j) = mbbp;
        end
    catch
        continue 
    end
end

etabbp = squeeze(bbp_giop(:,6,:)); % this is bbp at 667 using diff etas.
% need to make it at 700 using diff etas. this is GIOP with diff etas

eta700= etabbp.*(667 ./ 700).^-eta1;

for i = 1:length(eta700)
    
   bias(i,1) = max(eta700(i,:))-min(eta700(i,:));
   bias(i,2) = (max(eta700(i,:))-min(eta700(i,:))) ./ mean(eta700(i,:));
end
   
   bias(bias<0 | bias > 4)=NaN;
   
h6=histogram(bias(:,2),12,'FaceColor','k','EdgeColor','k','FaceAlpha',0.5);
h6.Normalization = 'probability';
h6.BinWidth = 0.05;
ylabel(' Probability','interpreter','Latex')
xlabel('Fraction b$_{bp}$ varies, f($\eta$)','interpreter','latex')
ax = gca; ax.LineWidth=1; 
set(gca,'TickLabelInterpreter','latex','FontSize',28)
title('Effect of assumed $\eta$ on b$_{bp}$','interpreter','latex')
  axis([0 0.5 0 0.55]) 
   
   
   
   
   
   
   
   
