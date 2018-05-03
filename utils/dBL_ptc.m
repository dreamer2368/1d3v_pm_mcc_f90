function dW = dBL_ptc(xp0, xp1, N, L)

coder.extrinsic('spdiags');

[xp, I] = sort([xp0;xp1],'ascend');
drho = -1/N*ones(2*N,1);
drho(I>N) = 1/N;

A = spdiags([-ones(2*N,1) ones(2*N,1)],[-1 0],2*N,2*N);
A(1,2*N) = -1;
A = [A; -A];
b = zeros(2*N,1);
b(2:2*N) = xp(2:2*N)-xp(1:2*N-1);
b(1) = xp(1)+L - xp(2*N);
b = min(abs(b),L-abs(b)); b = repmat(b,[2,1]);

lb = -100*ones(2*N,1); ub = 100*ones(2*N,1);

options = optimoptions('linprog','MaxIter',1e7,'Algorithm','dual-simplex');

[j1, dW1, exitflag,output] = linprog(drho,A,b,[],[],lb,ub,[],options);
if( exitflag ~= 1 )
    output
end
[j2, dW2, exitflag,output] = linprog(-drho,A,b,[],[],lb,ub,[],options);
if( exitflag ~= 1 )
    output
end
[dW,I] = min([dW1,dW2]); dW = -dW;
j = [j1,j2]; j = j(:,I);

end
