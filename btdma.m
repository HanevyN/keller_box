function X = btdma(A,B,C,R)
% Block TriDiagonal Matrix Algorithm 
% following Cebici & Bradshaw 
% Physical and Computational Aspects of Convective Heat Transfer pg 394

    X = zeros(size(R));
    N = size(A,3); 

    Gamma = zeros(size(A));
    Delta = Gamma;
    w = zeros(size(R));

    Delta(:,:,1) = B(:,:,1);
    w(:,:,1) = R(:,:,1);

    for i = 2:N
        Gamma(:,:,i) = A(:,:,i)/Delta(:,:,i-1);
        Delta(:,:,i) = B(:,:,i) - Gamma(:,:,i)*C(:,:,i-1);
        w(:,:,i) = R(:,:,i) - Gamma(:,:,i)*w(:,:,i-1);
    end

    X(:,:,N) = Delta(:,:,N)\w(:,:,N);

    for i = N-1:-1:1
        X(:,:,i) = Delta(:,:,i)\(w(:,:,i) - C(:,:,i)*X(:,:,i+1)); 
    end
    
end