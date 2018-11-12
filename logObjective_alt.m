function varargout = logObjective_alt(x,y,idx,f,fprime,params)
%% Computes the negative log-likelihood function and its gradient
%
% Usage:   [F,G] = logObjectiveGeneral(x,y,idx,f,fprime)
% Inputs:  x - vectorized version of the matrix
%          y - observations
%          idx - vector of indices corresponding to indices of observations
%          f - link function (f(x[i]) gives the probability that y[i]=1)
%          fprime - gradient of f with respect to x
% Outputs: F - the negative log-likelihood of x
%          G - the gradient of F (optional)
%
d=params.order;
params.factors{d}=reshape(x,params.numc,params.nr);

T=vec(cpdgen(params.factors));

if nargout > 0
   F = -sum(log(y(idx).*min(max(f(T(idx)),1e-9),1-1e-9) - (y(idx)-1)/2));
                        
   if isnan(F)
       cxc=find(isnan(f(T(idx))));
       params.factors{1}
       params.factors{2}
       params.factors{3}
       y(idx(cxc))
       T(idx(cxc))
       %log(y(idx).*f(T(idx)) - (y(idx)-1)/2)
       fedf
   end
   varargout{1} = F;
end

if nargout > 1
   G = zeros(size(T));
   v = (min(max(f(T(idx)),1e-9),1-1e-9)+(y(idx)-1)/2);
%    if length(find(abs(v)<1e-8))>0
%        cxc=find(abs(v)<1e-8);
%        T(idx(cxc))
%        f(T(idx(cxc)))
%        y(idx(cxc))
%        'dfv********@#%*%*%*#f(@*%*@#%*#%*#*%*#%*#@*%*#@*%@#*%*#%*_#%%_%*_@%#*%*#%_#*_*%#@*#_#%_@#'
%    end
   w = -fprime(T(idx));
   G(idx) = w./v;
   G=reshape(G,params.numr,params.numc);
%    max(abs(G(:)))
    %G=min(G,0.1*params.alpha);
    % G=max(G,-0.1*params.alpha);
   L=params.L;
   Gvec=vec(G'*params.L);
   %max(abs(Gvec(:)))
 %  Gvec=Gvec/max(abs(Gvec(:)))*params.alpha;
 %      Gvec=min(Gvec,0.1*params.alpha);
 %      Gvec=max(Gvec,-0.1*params.alpha);
    %Gvec=sign(Gvec);
   % Gvec=Gvec+randn(size(Gvec))/(mean(abs(Gvec(:))*10));
   varargout{2} = Gvec;
end
