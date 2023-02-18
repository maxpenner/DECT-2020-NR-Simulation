function [Y, idx, prefactor] = Transmit_diversity_precoding_Y(N_TS)

    % matlab complains if it's not 1i but just i
    q = 1i;
    
    % The TS does not include any power scaling for transmit diversity.
    % Transmit diversity takes one spatial stream and maps it onto multiple antennas.
    % We use a prefactor >= 1 to increase the power per antenna.
    %   
    %   1 stream on 2 antennas -> prefactor = 1 as the matrix in Table 6.3.3.2-1 already incorporates a prefactor of sqrt(2), so (sqrt(2))^2=2 times the power in total
    %   1 stream on 4 antennas -> prefactor = sqrt(2) as the matrix in Table 6.3.3.2-2 already incorporates a prefactor of sqrt(2), so (sqrt(2)*sqrt(2))^2=4 times the power in total
    %   1 stream on 8 antennas -> prefactor = 2 as the matrix in Table 6.3.3.2-3 already incorporates a prefactor of sqrt(2), so (sqrt(2)*2)^2=8 times the power in total
    %
    % Later, either beamforming or the Antenna_mapper_scaler(), if no beamforming is used, will decrease the power.
    %
    
    % !!!!!!!!!!! NOT STANDART COMPLIANT !!!!!!!!!!!
    % !!!!!!!!!!! NOT STANDART COMPLIANT !!!!!!!!!!!
    % !!!!!!!!!!! NOT STANDART COMPLIANT !!!!!!!!!!!
    % !!!!!!!!!!! NOT STANDART COMPLIANT !!!!!!!!!!!
    %dont_use_preafactor_as_in_TS = true;
    dont_use_preafactor_as_in_TS = false;
    % !!!!!!!!!!! NOT STANDART COMPLIANT !!!!!!!!!!!
    % !!!!!!!!!!! NOT STANDART COMPLIANT !!!!!!!!!!!
    % !!!!!!!!!!! NOT STANDART COMPLIANT !!!!!!!!!!!
    % !!!!!!!!!!! NOT STANDART COMPLIANT !!!!!!!!!!!
    
    prefactor = 1.0;

    switch N_TS
        case 2
            Y0 = [  1  0  q  0;...
                    0 -1  0  q;...
                    0  1  0  q;...
                    1  0 -q  0;...
                ];
            idx0 = [1, 2];
            
            if dont_use_preafactor_as_in_TS == true
                Y = Y0;
            else
                Y = prefactor * Y0;
            end
            
            idx = zeros(2,2);
            idx(1,:) = idx0;
            idx(2,:) = idx0;
            
        case 4
            Y0 = [  1  0  q  0;...
                    0 -1  0  q;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  1  0  q;...
                    1  0 -q  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                ];
            idx0 = [1, 2];
            
            Y1 = [  0  0  0  0;...
                    0  0  0  0;...
                    1  0  q  0;...
                    0 -1  0  q;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  1  0  q;...
                    1  0 -q  0;...
                ];
            idx1 = [3, 4];
            
            Y2 = [  1  0  q  0;...
                    0  0  0  0;...
                    0 -1  0  q;...
                    0  0  0  0;...
                    0  1  0  q;...
                    0  0  0  0;...
                    1  0 -q  0;...
                    0  0  0  0;...
                ];
            idx2 = [1, 3];
            
            Y3 = [  0  0  0  0;...
                    1  0  q  0;...
                    0  0  0  0;...
                    0 -1  0  q;...
                    0  0  0  0;...
                    0  1  0  q;...
                    0  0  0  0;...
                    1  0 -q  0;...
                ];
            idx3 = [2, 4];
            
            Y4 = [  1  0  q  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0 -1  0  q;...
                    0  1  0  q;...
                    0  0  0  0;...
                    0  0  0  0;...
                    1  0 -q  0;...
                ];
            idx4 = [1, 4];
            
            Y5 = [  0  0  0  0;...
                    1  0  q  0;...
                    0 -1  0  q;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  1  0  q;...
                    1  0 -q  0;...
                    0  0  0  0;...
                ];
            idx5 = [2, 3];
            
            Y = zeros(8,4,6);
            if dont_use_preafactor_as_in_TS == true
                Y(:,:,1) = Y0;
                Y(:,:,2) = Y1;
                Y(:,:,3) = Y2;
                Y(:,:,4) = Y3;
                Y(:,:,5) = Y4;
                Y(:,:,6) = Y5;
            else
                prefactor = sqrt(2);
                Y(:,:,1) = prefactor * Y0;
                Y(:,:,2) = prefactor * Y1;
                Y(:,:,3) = prefactor * Y2;
                Y(:,:,4) = prefactor * Y3;
                Y(:,:,5) = prefactor * Y4;
                Y(:,:,6) = prefactor * Y5;               
            end
            
            idx = zeros(12,2);
            idx(1,:) = idx0;
            idx(2,:) = idx0;
            idx(3,:) = idx1;
            idx(4,:) = idx1;
            idx(5,:) = idx2;
            idx(6,:) = idx2;
            idx(7,:) = idx3;
            idx(8,:) = idx3;
            idx(9,:) = idx4;
            idx(10,:) = idx4;
            idx(11,:) = idx5;
            idx(12,:) = idx5;
            
        case 8
            Y0 = [  1  0  q  0;...
                    0 -1  0  q;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  1  0  q;...
                    1  0 -q  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                ];
            idx0 = [1, 2];
            
            Y1 = [  0  0  0  0;...
                    0  0  0  0;...
                    1  0  q  0;...
                    0 -1  0  q;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  1  0  q;...
                    1  0 -q  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                ];
            idx1 = [3, 4];
            
            Y2 = [  0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    1  0  q  0;...
                    0 -1  0  q;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  1  0  q;...
                    1  0 -q  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                ];
            idx2 = [5, 6];
            
            Y3 = [  0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    1  0  q  0;...
                    0 -1  0  q;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  1  0  q;...
                    1  0 -q  0;...
                ];
            idx3 = [7, 8];
            
            Y4 = [  1  0  q  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0 -1  0  q;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  1  0  q;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    1  0 -q  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                ];
            idx4 = [1, 5];
            
            Y5 = [  0  0  0  0;...
                    1  0  q  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0 -1  0  q;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  1  0  q;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    1  0 -q  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                ];
            idx5 = [2, 6];
            

            Y6 = [  0  0  0  0;...
                    0  0  0  0;...
                    1  0  q  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0 -1  0  q;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  1  0  q;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    1  0 -q  0;...
                    0  0  0  0;...
                ];
            idx6 = [3, 7];
            
            Y7 = [  0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    1  0  q  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0 -1  0  q;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  1  0  q;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    1  0  -q  0;...
                ];
            idx7 = [4, 8];
            
            Y8 = [  1  0  q  0;...
                    0  0  0  0;...
                    0 -1  0  q;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  1  0  q;...
                    0  0  0  0;...
                    1  0 -q  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                ];
            idx8 = [1, 3];
            
            Y9 = [  0  0  0  0;...
                    1  0  q  0;...
                    0  0  0  0;...
                    0 -1  0  q;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  1  0  q;...
                    0  0  0  0;...
                    1  0 -q  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                ];
            idx9 = [2, 4];
            
            Y10 = [ 0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    1  0  q  0;...
                    0  0  0  0;...
                    0 -1  0  q;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  1  0  q;...
                    0  0  0  0;...
                    1  0 -q  0;...
                    0  0  0  0;...
                ];
            idx10 = [5, 7];
            
            Y11 = [ 0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    1  0  q  0;...
                    0  0  0  0;...
                    0 -1  0  q;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  0  0  0;...
                    0  1  0  q;...
                    0  0  0  0;...
                    1  0 -q  0;...
                ];
            idx11 = [6, 8];
            
            Y = zeros(16,4,12);
            if dont_use_preafactor_as_in_TS == true
                Y(:,:,1) = Y0;
                Y(:,:,2) = Y1;
                Y(:,:,3) = Y2;
                Y(:,:,4) = Y3;
                Y(:,:,5) = Y4;
                Y(:,:,6) = Y5;
                Y(:,:,7) = Y6;
                Y(:,:,8) = Y7;
                Y(:,:,9) = Y8;
                Y(:,:,10) = Y9;
                Y(:,:,11) = Y10;
                Y(:,:,12) = Y11;
            else
                prefactor = 2.0;
                Y(:,:,1) = prefactor * Y0;
                Y(:,:,2) = prefactor * Y1;
                Y(:,:,3) = prefactor * Y2;
                Y(:,:,4) = prefactor * Y3;
                Y(:,:,5) = prefactor * Y4;
                Y(:,:,6) = prefactor * Y5;
                Y(:,:,7) = prefactor * Y6;
                Y(:,:,8) = prefactor * Y7;
                Y(:,:,9) = prefactor * Y8;
                Y(:,:,10) = prefactor * Y9;
                Y(:,:,11) = prefactor * Y10;
                Y(:,:,12) = prefactor * Y11;
            end
            
            idx = zeros(24,2);
            idx(1,:) = idx0;
            idx(2,:) = idx0;
            idx(3,:) = idx1;
            idx(4,:) = idx1;
            idx(5,:) = idx2;
            idx(6,:) = idx2;
            idx(7,:) = idx3;
            idx(8,:) = idx3;
            idx(9,:) = idx4;
            idx(10,:) = idx4;
            idx(11,:) = idx5;
            idx(12,:) = idx5;
            idx(13,:) = idx6;
            idx(14,:) = idx6;
            idx(15,:) = idx7;
            idx(16,:) = idx7;
            idx(17,:) = idx8;
            idx(18,:) = idx8;
            idx(19,:) = idx9;
            idx(20,:) = idx9;
            idx(21,:) = idx10;
            idx(22,:) = idx10;
            idx(23,:) = idx11;
            idx(24,:) = idx11;
    end
end

