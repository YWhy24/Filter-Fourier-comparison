function filtered = laplace1d(matrix)

    N = length(matrix);
    filtered = zeros(1,N);
    for i = 2:N-1
        filtered(i) = matrix(i+1)+matrix(i-1)-2*matrix(i);
    end

end

