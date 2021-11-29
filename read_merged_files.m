
files = '/data1/bisson/lidar/MELD/merged_data/meld*old.mat'
[e,r]=unix(['ls -1 ',files]);

% Get mat files into workspace within cell structure
if (~e)
   r(end)=[];  % kill new line at end
   cellLs=eval([abs([123,39]),strrep(r,char(10),char([39,59,39])), ...
                abs([39,125])]); 
else
  warning(r)
  cellLs={};
end


for j = 1:length(cellLs)
    
load(cellLs{j})

str1=cellLs{j};
yr = str2num(str1(43:46));
mo = str2num(str1(47:48));

wl = [412 443 488 531 547 667];

gopt.aw  = get_aw(wl);
gopt.bbw = get_bbw(wl);
cm = get_oc(Rrs_r(:,2),Rrs_r(:,3),-1,Rrs_r(:,5),'oc3m');

for i = 1:length(Rrs_r)
    try
[x,mapg,maph,madg,mbbp,mrrs] = giop(wl,Rrs_r(i,:),cm(i),gopt);
    catch
        continue 
    end
bbp_giop(i,:) = mbbp;
end

bbp_giop(bbp_giop<0)=NaN;
bbp_giop(imag(bbp_giop) ~= 0) = 0;

eta = 2 * (1 - 1.2 * exp(-0.9 * Rrs_r(:,2)./Rrs_r(:,5)));

bbp = bbp_giop(:,2).*(532/443).^-eta;

ocb(j,1) = nanmedian(bbp);
lidarb(j,1) = nanmedian(lidar1deg(:,3)).*0.5; % phase function diff
lid8(j,1) =  datenum(datetime(yr, mo,1));
clear bbp_giop eta
end

save('/data1/bisson/lidar/MELD/ts_bbp.mat','ocb','lidarb','lid8')


