function [start_end] = fit_line_tls(n, c, pixellist) 
    d = n(1) *  pixellist(:, 1) + n(2) *  pixellist(:, 2) + c;
    Li_proj = pixellist - [n(1) * d, n(2) * d];    
%     
%     [~, min_i] = min(Li_proj(:,1));
%     [~, max_i] = max(Li_proj(:,1));
%     
    tau = pixellist * [-n(2); n(1)];
    [~, tau_min_i] = min(tau);
    [~, tau_max_i] = max(tau);
    
    start_end = [Li_proj(tau_min_i,1), Li_proj(tau_min_i,2);
                 Li_proj(tau_max_i,1), Li_proj(tau_max_i,2)];
end