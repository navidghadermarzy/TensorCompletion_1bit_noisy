function T=cp_to_ten(params,x)
if nargin>1
    params.factors{params.current_factor}=reshape(x,params.dimensions(params.current_factor),params.rank);
end
if(params.order==1)
    T=params.factors{1};
else
    T=cpdgen(params.factors);
end
    
