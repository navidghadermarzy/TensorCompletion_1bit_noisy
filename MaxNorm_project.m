function varargout = MaxNorm_1bit_project_active(x, B,params,projTol,projData)
projData=[];
%%%% Force the rows of L and R to have norm at most B.
%
% (LOut, ROut) = arg min_(U,V) || [L;R]-[U;V] ||_F^2 s.t. || [U;V] ||_mr<B
%
% Where || A ||_mr is the maximum Euclidean norm of a row of A.
%

X = reshape(x,params.numc,params.nr);
d=params.order;
passiveranks=setdiff(1:params.nr,params.activeranks);

for k = 1:params.numc
    Bold=norm(X(k,passiveranks));
    
    if norm(X(k,:)) > B
        X(k,params.activeranks) = X(k,params.activeranks) *(sqrt(B^2-Bold^2))/norm(X(k,params.activeranks));
    end
end
% % 

% params.factors{d}=X;
% T=vec(cpdgen(params.factors));
% scale_fac=max(abs(T(:)))/params.alpha;
% if scale_fac>1
%    X=X/scale_fac;
% end
%   


varargout{1} = X(:);
if (nargout > 1)
   projData    = []; % Dummy solution
   varargout{2} = projData;
end
end