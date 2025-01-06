function [x,rv,k]=PCG(A,b,x0,tol,kmax,L)
    rv=zeros(kmax+1,1);

    r=b-A*x0;
    rv(1)=norm(r);
    p=L'\(L\r);
    k=0;
    x=x0;
    c=r'*p;

    while norm(r)>tol*norm(b) && k<kmax
        z=A*p;
        alpha=c/(p'*z);
        x=x+alpha*p;
        r=r-alpha*z;
        g=L'\(L\r);
        d=r'*g;
        beta=d/c;
        c=d;
        p=g+beta*p;
        k=k+1;
        rv(k+1)=norm(r);
    end

    if k==kmax
        disp("PCG: ho raggiunto il massimo numero di iterazioni")
    end

    rv=rv(1:k+1);
end