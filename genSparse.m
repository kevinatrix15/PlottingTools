function A = genSparse(n,nz)
    density = nz /(n*n);
    A = sprandsym(n,density); % generate a random n x n matrix
    A = A + n * speye(n);
    [i,j,s] = find(A);
    nnZ = nnz(A);
    fid = fopen('C:\Users\Kevin\OneDrive\3DSIM\A.mtx','Wt');
    fprintf(fid, '%%%%MatrixMarket matrix coordinate real general\n');
    fprintf(fid, '%d %d %d\n', n, n, nnZ);    
    for k = 1:nnZ   
        fprintf(fid, '%d %d %f\n', i(k), j(k), s(k));
    end
    fclose(fid);
end
