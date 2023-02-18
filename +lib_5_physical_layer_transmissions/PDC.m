%
%   PDC is mapped across all spatial streams 0. The mapping is the same for each stream.
%
%   physical_resource_mapping_PDC_cell is always a cell(1,x)
%
%       physical_resource_mapping_PDC_cell(1,1)
%       physical_resource_mapping_PDC_cell(1,2)
%       physical_resource_mapping_PDC_cell(1,3)
%       ...
%       physical_resource_mapping_PDC_cell(1,x-2) contain subcarrier indices for all ofdm symbols where the PDC is placed
%
%       physical_resource_mapping_PDC_cell(1,x-1) contains the corresponding ofdm symbol indices
%
%       physical_resource_mapping_PDC_cell(1,x) contains linear indices for fast demapping
%
%   These indices are the same for each spatial stream.
%
%   example size fo Figure 4.5-2 b):
%
%        physical_resource_mapping_PDC_cell =
%
%          1×7 cell array
%
%            {56×1 double}    {56×1 double}    {56×1 double}    {42×1 double}    {56×1 double}    {56×1 double}    {1×6 double}    {322×1 double}
%
%
function [physical_resource_mapping_PDC_cell, N_PDC_subc] = PDC(u, numerology, k_b_OCC, N_PACKET_symb, N_TS, N_eff_TX,...
                                                                physical_resource_mapping_STF_cell,...
                                                                physical_resource_mapping_DRS_cell,...
                                                                physical_resource_mapping_PCC_cell)

    % Technical Specification assumes first index is 0, matlab 1
    MATLAB_INDEX_SHIFT = 1;
    
    % copy all relevant variables
    N_b_DFT = numerology.N_b_DFT;
    N_b_OCC = numerology.N_b_OCC;

    [~, mat_STF_DRS_PCC_all_streams] = lib_util.matrix_STF_DRS_PCC_PDC(N_b_DFT, N_PACKET_symb, N_TS, [],...
                                                                        physical_resource_mapping_STF_cell,...
                                                                        physical_resource_mapping_DRS_cell,...
                                                                        physical_resource_mapping_PCC_cell,...
                                                                        []);

    %% 5.2.5
    
    switch u
        case 1
            N_GI_plus_STF_symb = 2;
        case {2,4}
            N_GI_plus_STF_symb = 3;            
        case 8
            N_GI_plus_STF_symb = 4;            
    end
    
    N_DF_symb = N_PACKET_symb - N_GI_plus_STF_symb;
    
    if N_eff_TX <= 2
        N_step = 5;
    elseif N_eff_TX >= 4
        N_step = 10;
    else
        error('N_eff_TX is %d, cannot determine N_step according to 5.2.3.', N_eff_TX);
    end    

    nof_OFDM_symbols_carying_DRS = floor(N_PACKET_symb/N_step);

    % this is an error in the standard
    if N_step == 10 && mod(N_PACKET_symb, 10) ~= 0
        
        % sanity check
        if mod(N_PACKET_symb, 5) ~= 0
            error("N_PACKET_symb not a multiple of 5 or 10.");
        end
        
        nof_OFDM_symbols_carying_DRS = nof_OFDM_symbols_carying_DRS + 1;
    end

    N_DRS_subc = N_eff_TX * N_b_OCC/4 * nof_OFDM_symbols_carying_DRS;
    
    % according to 5.2.4
    N_PCC_subc = 98;
    
    N_PDC_subc = N_DF_symb * N_b_OCC - N_DRS_subc - N_PCC_subc;

    %% extract subcarrier and symbol index for PDC
    
    % we save data in one cell
    %
    %   PDC in mapped into subcarriers which are not used in transmit streams by DRS or PCC in spatial stream 0
    %
    %   number of rows: always 1
    %   first cell:     subcarrier indices for one ofdm symbol
    %   second cell:    subcarrier indices for one ofdm symbol
    %   third cell:     subcarrier indices for one ofdm symbol
    %   ...
    %   last cell:   	indices of all aforementioned ofdm symbols (see ofdm_symbol_indices)
    %    
    physical_resource_mapping_PDC_cell = cell(0);
    ofdm_symbol_indices = [];
    
    % first ofdm symbol at index 0 is STF, so we have to start at index 1
    l = 1;
    
    % counter for found PDC subcarriers
    cnt_PDC_subc = 0;
    
    % loop over each ofdm sybol in Data Field (DF)
    for q = 1:1:N_DF_symb
        
        % find all indices of unused subcarriers in current OFDM-symbol
        [k_i_matlab, ~] = find(~mat_STF_DRS_PCC_all_streams(:,l + MATLAB_INDEX_SHIFT));
        
        % convert to notation of technical specification
        k_i = lib_util.index_conversion_TS_matlab(N_b_DFT, k_i_matlab);
        
        % remove any subcarriers not within k_b_OCC, namely guards and DC
        k_i = intersect(k_i, k_b_OCC);
        
        % no free subcarriers in this symbol, move to next ofdm symbol
        if numel(k_i)== 0
            l = l + 1;
            continue;
        end
        
        % sort from lowest subcarrier to highest
        k_i = sort(k_i);

        % append subcarriers for this ofdm symbol
        physical_resource_mapping_PDC_cell = [physical_resource_mapping_PDC_cell {k_i}];
        cnt_PDC_subc = cnt_PDC_subc + numel(k_i);
                
        % remember current ofdm symbol
        ofdm_symbol_indices = [ofdm_symbol_indices, l];

        % move to next of ofdm symbol
        l = l + 1;
    end   
    
    %% create a vector with linear indices
    linear_indices_matlab = zeros(cnt_PDC_subc, 1);
    idx = 1;
    for i=1:1:numel(ofdm_symbol_indices)
        
        k_i = cell2mat(physical_resource_mapping_PDC_cell(i));
        n_k_i = numel(k_i);
        
        % linear indices
        rows = lib_util.index_conversion_TS_matlab(N_b_DFT,k_i);
        cols = repmat(ofdm_symbol_indices(i) + MATLAB_INDEX_SHIFT, numel(k_i), 1);
        li_matlab = sub2ind([N_b_DFT N_PACKET_symb], rows, cols);
        
        linear_indices_matlab(idx : idx + n_k_i - 1) = li_matlab;
        
        idx = idx + n_k_i;
    end  
    
    % sanity check
    if idx-1 ~= cnt_PDC_subc
        error('Incorrect number of linear indices.');
    end
    
    %% write result
    physical_resource_mapping_PDC_cell = [physical_resource_mapping_PDC_cell {ofdm_symbol_indices} {linear_indices_matlab}];
end

