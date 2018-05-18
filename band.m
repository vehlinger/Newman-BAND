function [C,E] = band(J,A,B,C,D,E,G,X,Y,N,NJ)

NP1 = N+1;

if J == 1
    for I = 1:N
        D(I,2*N+1) = G(I);
        for L = 1:N
            LPN = L+N;
            D(I,LPN) = X(I,L);
        end
    end
    [det,D] = matinv(N,2*N+1,B,D);
    if det == 0
        fprintf('Determinant = 0 at j = %d\n',J)
        return
    else
        for K = 1:N
            E(K,NP1,1) = D(K,2*N+1);
            for L = 1:N
                E(K,L,1) = -D(K,L);
                LPN = L+N;
                X(K,L) = -D(K,LPN);
            end
        end
        return
    end
end

if J == 2 
    for I = 1:N
        for K = 1:N
            for L = 1:N
                D(I,K) = D(I,K) + A(I,L)*X(L,K);
            end
        end
    end
end

if J == NJ
    for I = 1:N
        for L = 1:N
            G(I) = G(I)-Y(I,L)*E(L,NP1,J-2);
            for M = 1:N
                A(I,L) = A(I,L)+Y(I,M)*E(M,L,J-2);
            end
        end
    end
end

for I = 1:N
    D(I,NP1) = -G(I);
    for L = 1:N
        D(I,NP1) = D(I,NP1)+A(I,L)*E(L,NP1,J-1);
        for K = 1:N
            B(I,K) = B(I,K)+A(I,L)*E(L,K,J-1);
        end
    end
end

[det,D] = matinv(N,NP1,B,D);
if det == 0
    fprintf('Determinant = 0 at j = %d\n',J)
    return
else
    for K = 1:N
        for M = 1:NP1
            E(K,M,J) = -D(K,M);
        end
    end
end

if J < NJ
    return
else
   for K = 1:N
        C(K,J) = E(K,NP1,J);
    end
    for JJ = 2:NJ
        M = NJ-JJ+1;
        for K = 1:N
            C(K,M) = E(K,NP1,M);
            for L = 1:N
                C(K,M) = C(K,M)+E(K,L,M)*C(L,M+1);
            end
        end
    end
    for L = 1:N
        for K = 1:N
            C(K,1) = C(K,1)+X(K,L)*C(L,3);
        end
    end
end

end