function T = cpdgen(U,S)
%CPDGEN Generate full tensor given a polyadic decomposition.
%   T = cpdgen(U) computes the tensor T as the sum of R rank-one tensors
%   defined by the columns of the factor matrices U{n}.
%
%   See also btdgen, lmlragen, cpdres.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
if nargin < 2, S=ones(size(U{1},2),1); end
U{1}=U{1}*diag(S);
T = reshape(U{1}*kr(U(end:-1:2)).',cellfun('size',U(:).',1));
