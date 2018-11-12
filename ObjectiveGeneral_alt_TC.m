function varargout = ObjectiveGeneral_alt_TC(x,y,idx,params)
d=params.order;
params.factors{d}=reshape(x,params.numc,params.nr);

%T=vec(cpdgen(params.factors));
T=vec(conj(params.L)*(conj(params.factors{d}')));

err=T-y;
%norm(err(idx))
if nargout > 0
    
   %F = -sum(log(y(idx).*min(max(f(T(idx)),1e-7),1-1e-7) - (y(idx)-1)/2));
   F = (norm(err(idx)))^2/2;
   %F = (norm(err(idx)));
  
                        
   if isnan(F)
       err
       x
       
       params.factors{1}
       params.factors{2}
       params.factors{3}
       params.factors{4}
       y(idx(cxc))
       T(idx(cxc))
       %log(y(idx).*f(T(idx)) - (y(idx)-1)/2)
       fedf
   end
   varargout{1} = F;
end

if nargout > 1
   G = zeros(size(T));
%    v = err(idx);
%    if length(find(abs(v)<1e-8))>0
%        cxc=find(abs(v)<1e-8);
%        T(idx(cxc))
%        f(T(idx(cxc)))
%        y(idx(cxc))
%        'dfv********@#%*%*%*#@*%*@#%*#%*#*%*#%*#@*%*#@*%@#*%*#%*_#%%_%*_@%#*%*#%_#*_*%#@*#_#%_@#'
%    end
%    w = -fprime(T(idx));
   G(idx) = err(idx);
   G=(reshape(G,params.numr,params.numc));
   doinv=0;%rand>0.5;
   L=params.L;
   if doinv
       varargout{2}=vec((conj(G')*(params.Linv)'));
   else
    varargout{2}=(vec(conj(G')*(params.L)));
   end
   %MM = conj(vec(G'*params.L));
  
   %varargout{2} = conj(vec((params.Linv*G)'));
   %varargout{2}=conj(vec(G'*params.L));
   
end
