function series = fgn_sim(n, H)
% SERIES = FGN_SIM(N, H)
%
% Fractional Gaussian Noise Simulation
% 
% fgn_sim(n, H) create a time series of fractional Gaussian noise. This
% implementation is based on a function implemented by Diethelm Wuertz in
% R and the original algorithm in Beran (1994).
% 
% Input Parameters:
% n is an integer indicating the desired length of the series
% H is a real number on the interval (0, 1) and represents the desired Hurst
% exponent for the output series.
%
% Output Parameters:
% Series is a real valued time series of length n with Hurst expoent equal
% to H.
% 
%   Example:
%       x = fgn_sim(1000, 0.9);
%       plot(x);
% Author:  Author: Aaron D. Likens (2022)
%
% References:
% Beran, J. (1994). Statistics for long-memory processes. Chapman & Hall.


        % FUNCTION:
        
        % Settings:
        mu = 0; % output mean
        sd = 1; % output standard deviation
        
        % Generate Sequence:
        z = randn(1,2*n);
        zr = z(1:n);
        
        zi = z((n+1):(2*n));
        zic = -zi;
        zi(1) = 0;
        zr(1) = zr(1)*sqrt(2);
        zi(n) = 0;
        zr(n) = zr(n)*sqrt(2);
        zr = [zr(1:n), zr(fliplr(2:(n-1)))];
        zi = [zi(1:n), zic(fliplr(2:(n-1)))];
        z = complex(zr, zi);
        
        k = 0:(n-1);
        gammak = ((abs(k-1).^(2*H))-(2*abs(k).^(2*H))+(abs(k+1).^(2*H)))/2;
        ind = [0:(n - 2), (n - 1), fliplr(1:(n-2))];
        gkFGN0 = ifft(gammak(ind + 1))*length(z); % needs to non-normalized
        gksqrt = real(gkFGN0);
        if (all(gksqrt > 0))
            gksqrt = sqrt(gksqrt);  
            z = z.*gksqrt;
            z = ifft(z)*length(z); 
            z = 0.5*(n-1).^(-0.5)*z;
            z = real(z(1:n));
        else
            gksqrt = 0*gksqrt;
            Error("Re(gk)-vector not positive")
        end
        
        % Standardize:
        % (z-mean(z))/sqrt(var(z))
        series = sd*z + mu;

end
