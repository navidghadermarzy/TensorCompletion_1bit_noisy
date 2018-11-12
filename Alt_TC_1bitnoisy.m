function [T_constrained,params_opt] = Alt_TC_1bitnoisy(TE, E, rmin,rmax, f, fprime, Indc,indcc,iter,alpha,constraint_bound,paramsin,alt,pcc,doing1bit,doingmax)

%main_noisyor1bit_TC: outputs T_constrained from the measurements TE when
%constrained to the bound  "constraint_bound"

%Inputs:
%TE:       The mesurements tensor on the observed indices(E) can be 1-bit or noisy
%E:        Binary tensor indicaiting observed indices.
%rmin,rmax:Estimation 
%p_cc:     Fraction of observed sampled used for cross validation
%doingmax: Set to 1 for max-qnorm constrained TC. 0 for factor-fron constrained TC
%doing1bit:Set to 1 for 1-bit measurements. 0 for noisy(or clean partial) measurements. 

%f, fprime:f and fprime are used for defining obsevarions in 1-bit TC 
%check https://arxiv.org/abs/1804.00108 : Learning tensors from partial binary measurements

%iter:     Number of iteration for each alternation subproblem
%alternations: Number of alternation
%infbound .  : Infinity norm of the original tensor (this is just useful
%when used for recoverung the original tensor in 1bit TC) don't use for
%recovering the sign of the tensor 

doingnuclear=1-doingmax;
rank_cur=rmin;
Indcin=Indc;
iterinit=iter;
TE_init=TE;

dim=length(size(TE));

activeinds=1:rmin; %%In the process usually some of the columns of the factors go to zero in these cases, its makes th ealgorithm to remove those columns in all the factors


params.nr = rank_cur;
params.alpha=alpha;
params.rank_cur=rank_cur;
params.dimensions=size(TE);
ind=setdiff(1:prod(params.dimensions),Indc);
params.Ind = ind;



params.order=length(params.dimensions);
params.mode=1;
params.ls=1;
params.logical=0;

addpath(pwd);


if length(indcc)>0
    
    params.Indc=Indc;
    
    num_obs=length(Indc);
   
    num_cc=floor(num_obs*pcc);
    
    indcc=Indc(1:num_cc);
    
    Indc=Indc(num_cc+1:end);
    
    TE(indcc)=0;
    
    ind=[ind indcc];
end


params.Ind=ind;
params.Indc=Indc;

if paramsin.init==-1
    
    for i=1:params.order
        params.factors{i}=((randn(params.dimensions(i),params.nr)));
        X=params.factors{i};
        B=constraint_bound;
        for k = 1:params.dimensions(i)
            if norm(X(k,:))>B
                X(k,:) = X(k,:)*B/norm(X(k,:));
            end
        end
        params.factors{i}=X;
    end
    
else
    
    
    params.factors=paramsin.factors;
    
    for i=1:params.order
        params.factors{i}(:,activeinds)=randn(params.dimensions(i),length(activeinds));
    end
    params_opt=paramsin;
end



T_constrained=cp_to_ten(params);


if length(Indcin)==0
    errorindcc=-1000000;
else
    
end

f_ccopt=1000000;
params.activeranks=activeinds;

%% Begin Alternations
for i = 1:alt
     
            if  i>1 && i<6 %balancing the factors
    
                mult=prod(vv)^(1/params.order);
                for j=1:params.order
                    params.factors{j}=params.factors{j}*(mult/vv(j));
                    vv(j)=mult;
                end
            end
    
    for ii=params.nr:-1:activeinds(1)
        for j=1:params.order
            normvec(j)=norm(params.factors{j}(:,ii));
        end
        if prod(normvec)<1e-8
            activeinds=activeinds-(activeinds>ii);
            
            activeinds(ii-activeinds(1)+1)=[];
            params.activeranks=activeinds;
            for j=1:params.order
                params.factors{j}(:,ii)=[];
            end
            params.nr=params.nr-1;
            params.rank_cur=params.rank_cur-1;
            rank_cur=rank_cur-1;
        end
    end
   
    
    
    for j=1:params.order
        
        iter=iterinit;
        if j==1
            iter=iterinit;
        end
        params_cur=params;
        params_cur.dimensions=size(TE);
        params_cur.current_factor=j;
        row_dims = setdiff(1:params.order,j);
        y_cur=vec((matricize(TE,row_dims,params_cur.current_factor)));
        y_rot = vec(y_cur);
        
        params_cur.Ind = find(y_cur==0);
        
        
        swapfac=params_cur.factors{j};
        for shiftfac=j+1:params.order
            params_cur.factors{shiftfac-1}=params_cur.factors{shiftfac};
        end
        params_cur.factors{params.order}=swapfac;
        
        swapdim=params_cur.dimensions(j);
        for shiftfac=j+1:params.order
            params_cur.dimensions(shiftfac-1)=params_cur.dimensions(shiftfac);
        end
        params_cur.dimensions(params.order)=swapdim;
        
        
        %%
        params_cur.current_factor=params.order;
        reduced_sizes=params_cur.dimensions;
        reduced_sizes(params_cur.current_factor)=[];
        params_cur.numc=params_cur.dimensions(params.order);
        params_cur.numr=prod(reduced_sizes);
        B=zeros(prod(reduced_sizes),params_cur.rank_cur);
        for i2=params.activeranks
            for j2=1:length(params_cur.factors)
                if j2<params_cur.current_factor
                    cell_cols.factors{j2}=((params_cur.factors{j2}(:,i2)));
                elseif j2>params_cur.current_factor
                    cell_cols.factors{j2-1}=((params_cur.factors{j2}(:,i2)));
                end
            end
            
            cell_cols.order=params_cur.order-1;
            cell_cols.rank_cur=params_cur.rank_cur;
            
            B(:,i2)=(vec(cp_to_ten(cell_cols)));
            
        end
        params_cur.L = B;
        params_cur.Linv=pinv(B);
        tic
        params_cur.Indc=setdiff(1:prod(params_cur.dimensions),params_cur.Ind);
        
        
        Indcinit=params_cur.Indc;
        if mod(i,2)==3 && i<alt-3
            rp=randperm(length(Indcinit));
            Indcur=Indcinit(rp(1:ceil(length(Indcinit)/10)));
        else
            rp=randperm(length(Indcinit));
            Indcur=Indcinit(rp(1:ceil(length(Indcinit))));
        end
            
        params_cur.Indc=Indcur;%  union(Indcur,diff_inds);
        
        R1=conj(params_cur.factors{params.order});
        
        R1=R1(:,1:rank_cur);
        
        xinit = R1(:);
        
        if doingmax
        %%funproj for using spgsolver
        funProj_tensor = @(x,projTol,projData) MaxNorm_project(x,constraint_bound,params_cur,projTol,projData);
        %funProj_tensor = @(x) MaxNorm_project(x,constraint_bound,params_cur);
        else
            
        funProj_tensor = @(x,projTol,projData) NucNorm_project(x,constraint_bound,params_cur,projTol,projData);
        end
        if doing1bit
        funObj_tensor  = @(x) logObjective_alt(x,y_rot,params_cur.Indc,f,fprime,params_cur);
        else
        funObj_tensor  = @(x) ObjectiveGeneral_alt_TC(x,y_rot,params_cur.Indc,params_cur);
        end
        
        %warning off
        optionsPQN.optTol=1e-7;
        optionsPQN.progTol=1e-7;
        optionsPQN.corrections=30;
        optionsPQN.iterations=iter;
        optionsPQN.maxIter=iter;
        optionsPQN.suffDec=1e-7;
        optionsPQN.SPGoptTol=1e-7;
        optionsPQN.verbosity=0;
        optionsPQN.verbose=0;
        
        [xLS,info] = spgSolver(funObj_tensor, funProj_tensor, R1, optionsPQN);
        
        %[xLS,infp]=minConf_PQN(funObj_tensor, funProj_tensor, R1(:), optionsPQN);
        
        
        params_cur.factors{params.order}=reshape(xLS,params_cur.dimensions(params.order),params.rank_cur);
        params.factors{j}=((params_cur.factors{params.order}));
        
        T_constrained=cp_to_ten(params);
        T_constrained=min(T_constrained,params.alpha);
        T_constrained=max(T_constrained,-1*params.alpha);
        if doingmax
        vv(j)=MaxNorm_Primal(xLS,[],params_cur);
        else
        vv(j)=norm(xLS(:));
        end
        
        
 
        
        f_c=logObjectiveGeneral(T_constrained(:),TE(:),Indc,f,fprime);
        f_cc=logObjectiveGeneral(T_constrained(:),TE_init(:),indcc,f,fprime);
        if f_cc < f_ccopt
            params_opt=params;
            f_ccopt=f_cc;
        end
        
        
        
        f_cold=f_c;
        
    end
end

params_opt=params;
params_opt.init=0;
params_opt.res=prod(vv)^(1/params_opt.order);

T_constrained=cp_to_ten(params_opt);
end
