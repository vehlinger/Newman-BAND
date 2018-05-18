function [DET,D] = matinv(N,M,B,D)

DET = 1.0;
ID = zeros(1,N);

for NN = 1:N
    BMAX = 1.1;
    for I = 1:N
        if ID(I) == 0
            BNEXT = 0.0;
            BTRY = 0.0;
            for J = 1:N
                if ID(J) == 0 
                    if abs(B(I,J)) > BNEXT
                        BNEXT = abs(B(I,J));
                        if BNEXT > BTRY
                            BNEXT = BTRY;
                            BTRY = abs(B(I,J));
                            JC = J;
                        end
                    end
                end
            end        
            if BNEXT < BMAX*BTRY 
                BMAX = BNEXT/BTRY;
                IROW = I;
                JCOL = JC;
            end
        end
    end
    if ID(JC) == 0
        ID(JCOL) = 1;
    else
        DET = 0;
        return
    end
    if JCOL ~= IROW
        for J=1:N
            SAVE = B(IROW,J);
            B(IROW,J) = B(JCOL,J);
            B(JCOL,J) = SAVE;
        end
        for K=1:M
            SAVE = D(IROW,K);
            D(IROW,K) = D(JCOL,K);
            D(JCOL,K) = SAVE;
        end
    end
    F = 1/B(JCOL,JCOL);
    for J = 1:N
        B(JCOL,J) = B(JCOL,J)*F;
    end
    for K = 1:M
        D(JCOL,K) = D(JCOL,K)*F;
    end         
    for I=1:N
        if (I~=JCOL)
            F = B(I,JCOL);
            for J = 1:N
                B(I,J)=B(I,J)-F*B(JCOL,J);
            end
            for K=1:M
                D(I,K) = D(I,K)-F*D(JCOL,K);
            end
        end
    end
end
return

end         

                            
            
