function p = MaxNorm_Primal(x, weights,params)

X = reshape(x,params.numc,params.nr);
v = zeros(params.numc,1);

for k = 1:params.numc
    v(k) = norm(X(k,:));
end

p = max(v);

end