function dW = dBL_pm(xp0, rho1, N, Ng, L)

dx = L/Ng;
xg = (1:Ng)'*dx - 0.5*dx;

[xp, I] = sort([xp0;xg],'ascend');
drho = zeros(N+Ng,1);
drho(I<=N) = -1/N; drho(I>N) = rho1*dx;

A = spdiags([-ones(N+Ng,1) ones(N+Ng,1)],[-1 0],N+Ng,N+Ng);
A(1,N+Ng) = -1;
A = [A; -A];
b = zeros(N+Ng,1);
b(2:N+Ng) = xp(2:N+Ng)-xp(1:N+Ng-1);
b(1) = xp(1)+L - xp(N+Ng);
b = min(abs(b),L-abs(b)); b = repmat(b,[2,1]);

lb = -100*ones(N+Ng,1); ub = 100*ones(N+Ng,1);

options = optimoptions('linprog','MaxIter',1e3);

[j1, dW1, exitflag,output] = linprog(drho,A,b,[],[],lb,ub,[],options);
if( exitflag ~= 1 )
    output
end
[j2, dW2, exitflag,output] = linprog(-drho,A,b,[],[],lb,ub,[],options);
if( exitflag ~= 1 )
    output
end
[dW,minI] = min([dW1,dW2]); dW = -dW;
j = [j1,j2]; j = j(:,minI);

end