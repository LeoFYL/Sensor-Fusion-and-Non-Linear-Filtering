n = length(xf_k);
    switch type
        case 'EKF'
            [~, F] = f(xf_k, T);
            Pxf = Pf_k * F' * Pp_kplus1;
            
        case {'UKF', 'CKF'}
            [SP, W] = sigmaPoints(xf_k, Pf_k, type);
            SP_kplus1 = zeros(n, size(SP, 2));
            for i = 1:size(SP, 2)
                SP_kplus1(:, i) = f(SP(:, i), T);
            end
            
            Pxf = zeros(n, n);
            for i = 1:size(SP, 2)
                Pxf = Pxf + W(i) * (SP(:, i) - xf_k) * (SP_kplus1(:, i) - xp_kplus1)';
            end
            
        otherwise
            error('Incorrect type of non-linear smoother')
    end

    % Compute smoother gain
    Ks = Pxf / Pp_kplus1;

    % Update the smoothed state and covariance
    xs = xf_k + Ks * (xs_kplus1 - xp_kplus1);
    Ps = Pf_k - Ks * Pp_kplus1 * Ks';