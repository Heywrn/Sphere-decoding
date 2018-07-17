function funct(z, R, sym, k, dist, b_hat)
global bEst;
global radius;

if k == 1
    for ii = 1:2
        b_hat(k) = sym(ii);
        d = (z(1)- R(1,:)*b_hat)^2 + dist;
        if (d<=radius)
            bEst = b_hat;
            radius = d;
        end
    end
else
    for ii = 1:2
        b_hat(k) = sym(ii);
        d = abs(z(k) - R(k,[k:end])*b_hat(k:end))^2 + dist;
        if (d <= radius)
            funct(z, R, sym, k-1, d, b_hat);
        end
    end
end