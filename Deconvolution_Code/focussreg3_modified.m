function [x, J, Reg] = focussreg3_modified(x, b, G, K)
x = x(:);
b = b(:);
[m n] = size(G);
lambd=2.0e-3;
%K=0.001;
itermax = 1000;
for index_xiter=1:30
    W = diag(abs(x).^(3/2));
    x =W * G' * ((G * W* G' + lambd * eye(length(b))) \ b);
    for i = 1:n
        if(x(i) <= K)
            x(i) = 0.0;
        end
    end
    [U,S,V] = svd(G*sqrt(W));
    E=diag(U'*b*b'*U);
    vlambd=@(l)sum(E.*(l./(diag(S).^2+l^2)).^2)/(1/m*sum(l./(diag(S).^2+l^2))^2);
    lambd = fminbnd(vlambd, 0, 10);
    lambdaer(index_xiter)=lambd;
%     x2(index_xiter)=norm(x);
%     x1(index_xiter)=sum(abs(x));
%     x0(index_xiter)=sum(sign(x));
%     J(index_xiter)= norm(G*x - b);
    Reg=lambd;
    J= 1/2*norm(G*x - b)^2+Reg*(sum(abs(x).^(1/2)));
end
count = 0;
fprintf('lambda:  %d,  ', lambdaer(index_xiter));
while (round(lambdaer(index_xiter)*10^3)/10^3)~=(round(lambdaer(index_xiter-1)*10^3)/10^3)
    index_xiter = index_xiter+1;
    W = diag(abs(x).^(3/2));
    x =W * G' * ((G * W* G' + lambd * eye(length(b))) \ b);
    for i = 1:n
        if(x(i) < K)
            x(i) = 0.0;
        end
    end
    [U,S,V] = svd(G*sqrt(W));
    E=diag(U'*b*b'*U);
    vlambd=@(l)sum(E.*(l./(diag(S).^2+l^2)).^2)/(1/m*sum(l./(diag(S).^2+l^2))^2);
    lambd = fminbnd(vlambd, 0, 1);
    fprintf('%d,  ',lambd);
    lambdaer(index_xiter)=lambd;
%     x2(index_xiter)=norm(x);
%     x1(index_xiter)=sum(abs(x));
%     x0(index_xiter)=sum(sign(x));
%     J(index_xiter)= norm(G*x - b);
    Reg=lambd;
    J= 1/2*norm(G*x - b)^2+Reg*(sum(abs(x).^(1/2)));
    count = count+1;
    if(count>itermax)
        break;
    end
end
fprintf('\n');
end