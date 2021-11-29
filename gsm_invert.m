function IOPs = gsm_invert(rrs,aw,bbw,bbpstar,aphstar,admstar);
% 2013-02-19 changed initial guess to be close to PnB values.
% old values commented out.  Don't think it makes a bit of difference.

options = optimset('fminsearch');
options = optimset(options,'TolX',1e-9,'TolFun',1e-9, 'MaxFunEvals', 2000, 'MaxIter', 2000);
IOPs=repmat(NaN,[size(rrs,1),3]);

% use init for PnB
IOPSinit=[.02,.01,.0029]; % used for global.  changing this to see
%if there will be less max iterations reached.  
%IOPSinit=[1.5 .08 .0033];
    
for i = 1:size(rrs,1)
   rrs_obs = rrs(i,:);
   [iops,cost(i),exitFlag] = fminsearch('gsm_cost',IOPSinit, ...
					options,rrs_obs,aw,bbw,bbpstar,aphstar,admstar); 

   if (exitFlag==1) % if converged then use value as IOP and inital
                    % guess
     IOPs(i,:)=iops;  % 
%     IOPSinit=iops;
   end
   
   %% stephane uses [.002 .01 .0029]
   %NOTE THAT THE CHL RETRIEVALS ARE ESPECIALLY SENSITIVE TO ITS INITIAL
   %GUESS. IF WE USE 0.002 AS SUGGESTED< MANY NEGATIVES AND OTHER
   %OUTRAGEOUS STUFF OCCURS. 

end
% removed 2012-12-17 correction for using synthetic data- no longer needed.
%IOPs(:,2)=IOPs(:,2) .* 0.754188;    %%%%% CORRECTION FOR BIAS IN OPTIMIZED MODEL				    %%%%% SEE MARITORENA et al., (2002)
return