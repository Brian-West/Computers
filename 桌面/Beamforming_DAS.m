function [ image_data ] = Beamforming_DAS( image_data_v,image_data_t, N_points_column, d_z_start,d_z_size,d_e_width,d_e_kerf,N_active,c,fs,varargin)
%**************************************************************************
% Delay and sum beamforming algorithm
%   By Mr_Cjh
%       2016.08.3
%**************************************************************************

tic;

[~, N_elements, N_lines] = size(image_data_v);
d_z_point_spacing = d_z_size/N_points_column;
% d_z_point_spacing = d_z_size/(N_points_column-1);

image_data = zeros(N_points_column,N_lines);
weight_hanning = hanning(N_elements);
% weight_hanning = ones(N_elements,1);
for i=1:N_lines    
    data_scanLines = image_data_v(:,:,i); 
    % data_scanLines = [ zeros( round(image_data_t(i)*fs),N_elements ); image_data_v(:,:,i) ];
    for j=1:N_points_column
        for k=1:N_elements  %= i*(d_e_width+d_e_kerf)-(d_e_width+d_e_kerf)+( N_active*d_e_width + (N_active-1)*d_e_kerf )/2 ;
            dp_emit_middle = (i-1)*(d_e_width+d_e_kerf) + ( N_active*d_e_width + (N_active-1)*d_e_kerf )/2 ;    % the position of scanlines
            % dp_emit_middle = (i-1)*(d_e_width+d_e_kerf);% + ( N_active*d_e_width + (N_active-1)*d_e_kerf )/2 ;    % the position of scanlines
            % dp_emit_middle = i*(d_e_width+d_e_kerf) - d_e_kerf/2 ;          % the position of scanlines
            dp_receive = (k-1)*(d_e_width+d_e_kerf) + d_e_width/2;          % position of receiving elements
            ds_emit_receive = abs( dp_emit_middle - dp_receive );           % middle aperture distance between emit and receive
            ds_point = d_z_start+(j-0.5)*d_z_point_spacing;             % calculating point position
%             ds_point = d_z_start+(j-1)*d_z_point_spacing;               % calculating point position            
            ds_receive = sqrt( ds_emit_receive^2 + ds_point^2 );            % receiving distance
            dt_receive = (ds_point+ds_receive)/c;                           % receiving time
            di_point_value_k = round(dt_receive*fs - image_data_t(i)*fs );  % index of point value
            
            if di_point_value_k>0 && di_point_value_k<=length(data_scanLines)               
                image_data(j,i) = image_data(j,i) + ...
                            weight_hanning(k)*data_scanLines(di_point_value_k,k);
            end
                       
        end
        
    end
end

fprintf('DAS: ');  toc
    
end

