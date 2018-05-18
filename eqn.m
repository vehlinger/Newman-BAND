function [eq,C] = eqn(j,jp,k,dC,C,nj)
    
    % j = current mesh point
    % jp = differentiation point
    % k = unknown variable being varied
    % dC = amount variable k is varied by
    % C = matrix of unknowns
    
    % physical constants
    D12 = 8.48e-2;     % cm2/s
    D13 = 13.72e-2;    % cm2/s
    D23 = 19.91e-2;    % cm2/s
    ctot = 36.395;     % mol/cm3
    L = 238;           % modeling length (cm)
    h = L/nj;          % mesh spacing
    
    % store old C value
    C_save = C(k,jp);
    % increment unknown k
    C(k,jp) = C(k,jp)+dC;
    
    if j == 1
    % boundary condition 1: Dirichlet
        eq(1) = (C(1,j)*C(5,j)-C(2,j)*C(4,j))/(ctot*D12)...
            -(C(3,j)*C(4,j))/(ctot*D13)-(C(1,j+1)-C(1,j))/h;
        eq(2) = (C(2,j)*C(4,j)-C(1,j)*C(5,j))/(ctot*D12)...
            -(C(3,j)*C(5,j))/(ctot*D23)-(C(2,j+1)-C(2,j))/h;
        eq(3) = 1-C(1,j)-C(2,j)-C(3,j);
        eq(4) = 0.319-C(1,j);             % x1(0) = 0.319
        eq(5) = 0.528-C(2,j);             % x2(0) = 0.528
    elseif j == nj
    % boundary condition 2: Dirichlet
        eq(1) = C(1,j);                   % x1(L) = 0
        eq(2) = C(2,j);                   % x2(L) = 0
        eq(3) = 1-C(1,j)-C(2,j)-C(3,j);   % sum(xi) = 1
        eq(4) = C(4,j)-C(4,j-1);
        eq(5) = C(5,j)-C(5,j-1);
    else
        eq(1) = (C(1,j)*C(5,j)-C(2,j)*C(4,j))/(ctot*D12)...
            -(C(3,j)*C(4,j))/(ctot*D13)-(C(1,j+1)-C(1,j))/h;
        eq(2) = (C(2,j)*C(4,j)-C(1,j)*C(5,j))/(ctot*D12)...
            -(C(3,j)*C(5,j))/(ctot*D23)-(C(2,j+1)-C(2,j))/h;
        eq(3) = 1-C(1,j)-C(2,j)-C(3,j);
        eq(4) = C(4,j)-C(4,j-1);
        eq(5) = C(5,j)-C(5,j-1);
    end
    C(k,jp) = C_save;
end