function dW = dBL(rho0, rho1, Ng, dx)

%rho0, rho1 are normalized
L = Ng*dx;

drho = (rho1-rho0)*dx;

A = spdiags([-ones(Ng,1)/dx ones(Ng,1)/dx],[-1 0],Ng,Ng);
A(1,Ng) = -1/dx;
A = [A; -A];
b = ones(2*Ng,1);

lb = -100*ones(Ng,1); ub = 100*ones(Ng,1);

% options = optimoptions('linprog','Display','iter');

[j1, dW1, exitflag,output] = linprog(drho,A,b,[],[],lb,ub);
if( exitflag ~= 1 )
    output
end
[j2, dW2, exitflag,output] = linprog(-drho,A,b,[],[],lb,ub);
if( exitflag ~= 1 )
    output
end
[dW,I] = min([dW1,dW2]); dW = -dW;
j = [j1,j2]; j = j(:,I);

end