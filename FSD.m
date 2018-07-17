function [BFSD] = FSD(y, A, O, N)

% y is the y = H*B + noiseV*Noise vector with noise;
% A is the matrix of G,R,F.
% N is the length of the matrix C.


BFSD = sign(randn(N,1));

s = zeros(N+1);
% Initialize s;
L = N+1;
s(L,1) = 1;

j1 = 1;

% Fisrt, get the s(1) using the formula (17)
while(j1<L)
    mid = 0;
    for i=L:-1:L-j1+1
        mid = mid + A(:,L-j1)' * A(:,i) * s(i,1);
    end
    s(L-j1,1) = sign(mid);
    j1 = j1+1;
end


% Based on s(1), we can get s(n) using the formula (19) and (17)
n = 2;

while(n<=L)
    s(L-n+1,n) = -s(L-n+1,1);
    
    for i=L-n+2:L
        s(i,n) = s(i,1);
    end
    
    j = n;
    
    while(j<L)
        
        mid = 0;
        
        for i=L:-1:L-j+1
            mid = mid + A(:,L-j)' * A(:,i) * s(i,n);
        end
        s(L-j,n) = sign(mid);
        
        j = j+1;
    end
    
    n = n+1;
end

best = 0;
% calculate which s is the best, first put the best one as the 1st vector
be1 = norm(y - O * s(1:N,1),'fro');;
be2 = norm(y + O * s(1:N,1),'fro');;

if(be1 < be2)
    BFSD = s(1:N,1);
    best = be1;
else
    BFSD = -s(1:N,1);
    best = be2;
end

% Calculate which is the best one and put the best one in s
for q = 2:L
    
    Cal1 = norm(y - O * s(1:N,q),'fro');
    Cal2 = norm(y + O * s(1:N,q),'fro');
    
    if(best > min(Cal1,Cal2))
        
        best = min(Cal1,Cal2);
        
        if(Cal1 < Cal2)
            BFSD = s(1:N,q);
        else
            BFSD = -s(1:N,q);
        end
    end
    
end
