function B = gendiag(A,dims)
% generalized "diag" function
% dims is a row vector of dimensions that are to be merged

sz = size(A);

if max(dims) > length(sz)
    sz = cat(2,sz,ones(1,max(dims)-length(sz)));
end

assert(nnz(diff(sz(dims)))==0,'specified dimensions must be equal')

alldims = 1:length(sz);
otherdims = setdiff(alldims,dims);

B = permute(A,cat(2,otherdims,dims));
B = reshape(B,cat(2,prod(sz(otherdims)),sz(dims)));
M = logical(eye(sz(dims(1))));
B = B(:,M);

B = reshape(B,cat(2,sz(otherdims),sz(dims(1))));


end

