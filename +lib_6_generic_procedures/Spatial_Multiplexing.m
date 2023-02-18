function [y_PxC_ss] = Spatial_Multiplexing(x_PxC, N_SS)

    % sanity check
    if N_SS <= 1
        error('Spatial Multiplexing ony applicable for N_SS>1, is %d.', N_SS);
    end
        
    % 6.3.2
    M_symb = numel(x_PxC);

    % sanity check
    if mod(M_symb,N_SS) ~= 0
        error('M_symb is not a multiple of N_SS.');
    end

    M_stream_symb = M_symb/N_SS;

    % deinterleave into spatial streams
    y_PxC_ss = cell(N_SS,1);
    for i=1:1:N_SS

        x_s = x_PxC(i:N_SS:end);

        % sanity check
        if numel(x_s) ~= M_stream_symb
            error('Number of elements in single spatial stream is not M_stream_symb.');
        end            

        y_PxC_ss(i) = {x_s};
    end
end

