function [C,E,G] = autoband(j,n,nj,C,dC,E)

G = zeros(1,n);     % matrix of governing equations
A = zeros(n,n);     % matrix of dG/dC at j-1 
B = zeros(n,n);     % matrix of dG/dC at j
D = zeros(n,2*n+1); % matrix of dG/dC at j+1
X = zeros(n,n);     % matrix of dG/dC at j-2 
Y = zeros(n,n);     % matrix of dG/dC at j+2

% initialize G (k = 1, dC = 0)
eq = eqn(j,1,1,0,C,nj);
for i = 1:n
    G(i) = eq(i);
end

% generate A,B,D matrices
for k = 1:n
    eq = eqn(j,j,k,dC(k),C,nj);
    for i = 1:n
        B(i,k) = -(eq(i)-G(i))/dC(k);
    end
    if j == 1
        eq = eqn(j,j+2,k,dC(k),C,nj);
        for i = 1:n
            X(i,k) = -(eq(i)-G(i))/dC(k);
        end
    else
        eq = eqn(j,j-1,k,dC(k),C,nj);
        for i = 1:n
            A(i,k) = -(eq(i)-G(i))/dC(k);
        end
    end
    if j == nj
        eq = eqn(j,j-2,k,dC(k),C,nj);
        for i = 1:n
            Y(i,k) = -(eq(i)-G(i))/dC(k);
        end
    else
        eq = eqn(j,j+1,k,dC(k),C,nj);
        for i = 1:n
            D(i,k) = -(eq(i)-G(i))/dC(k);
        end
    end
end

% solve for C using band function
[C,E] = band(j,A,B,C,D,E,G,X,Y,n,nj);
end