%% demo.m
% 
% This file contains the code for a simple demo of the one bit matrix 
% completion software package that supplements the paper 
% http://arxiv.org/abs/1209.3672. 
%
% Most recent change - 5/16/2014
%
% Copyright 2013, M. Davenport, Y. Plan, E. van den Berg, M. Wootters
%
% This file is part of 1BMC Toolbox version 1.2.
%
%    The 1BMC Toolbox is free software: you can redistribute it and/or 
%    modify it under the terms of the GNU General Public License as 
%    published by the Free Software Foundation, either version 3 of the 
%    License, or (at your option) any later version.
%
%    The 1BMC Toolbox is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with the 1BMC Toolbox. If not, see <http://www.gnu.org/licenses/>.


%% Initialize
clear all;
close all;

d1 = 20; % Number of rows in matrix
d2 = 20; % Number of columns in matrix
d3=20;
r = 3; % Rank
pct = 25; % Percent of entries to observe
sigma = 0.0; % Noise level

%strm = RandStream('mt19937ar','Seed',5);

%% Generate a random low-rank matrix using uniform matrices
n = [100 150 200];
r1 = 4;
r2 = 5;
r3 = 6;
r = [r1, r2, r3];

% sparsity: p = nnz / (n^3)
p = 0.01;

opts = struct( 'maxiter', 100, 'tol', 1e-9 );

% set the seed of the mersenne twister to 11 for reproducible results
rng( 11 );

% create the data to complete:
% create random rank-r tensor ...
A = makeRandTensor( n, r );
% create the sampling set ...
subs = makeOmegaSet( n, round(p*prod(n)) );
% get the values of A at the sampling points ...
vals = getValsAtIndex(A, subs);

%% Define observation model (probit/Gaussian noise)
f      = @(x) gausscdf(x,0,sigma);
fprime = @(x) gausspdf(x,0,sigma);

% Logistic model
%f       = @(x) (1 ./ (1 + exp(-x)));
%fprime  = @(x) (exp(x) ./ (1 + exp(x)).^2);  

%% Obtain signs of noisy measurements
Y = sign(f(M)-rand(strm,d1,d2));
y = Y(:);
    
%% Observe 'pct' percent of the entries in Y
idx = find(rand(strm,d1*d2,1) <= pct/100);
m = length(idx);

% %% Set up optimization problem
% options = struct();
% options.iterations = 10000; 
% options.stepMax    = 10000;
% options.stepMin    = 1e-4;
% options.optTol     = 1e-3;
% options.stepMax    = 1e9;
%         
% funObj  = @(x) logObjectiveGeneral(x,y,idx,f,fprime);
% 
% %% Define alpha to be the correct maximum using an oracle
% alpha   = 1;
% radius  = alpha * sqrt(d1*d2*r);
% 
% %% Define constraints
% % Use nuclear-norm constraint only
% funProj = @(x,projTol,projData) projNucnorm(x,d1,d2,radius,projTol,projData);
% 
% % Use nuclear-norm plus infinity-norm constraints
% %funProj = @(x,projTol,projData) projectKappaTau(x,d1,d2,radius,alpha,projTol,projData);
% 
% %% Recover estimate Mhat of M
% [Mhat,info] = spgSolver(funObj, funProj, zeros(d1*d2,1), options);
% Mhat = reshape(Mhat,d1,d2);
% 
% [U,S,V] = svd(Mhat);
% Mhat_debias = U(:,1:r)*S(1:r,1:r)*V(:,1:r)'; % Project onto actual rank if known
% 
% %% Compute relative error
% norm(Mhat_debias-M,'fro')/norm(M,'fro')
