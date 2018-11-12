function varargout = NucNorm_project(x, B,params,projTol,projData)
projData=[];

X = reshape(x,params.numc,params.nr);
d=params.order;


    
 if norm(X(:))>B
        X = X*B/norm(X(:));
 end
    
params.factors{d}=X;
 T=vec(cpdgen(params.factors));
scale_fac=max(abs(T(:)))/params.alpha;

% if scale_fac>1
%     X=X/scale_fac;
% end

varargout{1} = X(:);
if (nargout > 1)
   projData    = []; % Dummy solution
   varargout{2} = projData;
end
end