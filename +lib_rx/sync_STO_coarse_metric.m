function metric = sync_STO_coarse_metric(window, M, L, detection_E_rms_threshold)
    
    % separate into repetitive pattern
    window = reshape(window, M, L);

    % power for normalization
    E = sum(sum(window.*conj(window)));

    % we have to have a minimal power to assume there is a frame
    E_mean = sqrt(E/(M*L));
    if E_mean < detection_E_rms_threshold
        metric = 0;
        return;
    end
    
    % correlation from pattern to pattern
    P = 0;
    for j=1:1:L-1
        P = P + sum(window(:,j).*conj(window(:,j+1)));
    end

    metric = (  (L/(L-1)) * abs(P)/E  )^2;

    % We assume the worst case:
    % We can't detect NaN and inf, instead we get maximum metric.
    if isnan(metric) || isinf(metric)
        metric = 1;
    end
end

