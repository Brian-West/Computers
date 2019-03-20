function [  ] = Beamforming_SimulatedData( Phantom_Name )
% function [  ] = Beamforming_SimulatedData( Phantom_Width,Phantom_Depth,Phantom_Depth_Start,Phantom_Thickness )
%**************************************************************************
% use field II to generate simulation data for test
%   By Mr_Cjh
%       2016.08.3
% 
% the position of scanlines is determined by active aperture_emit
% and there are 3 dimension in focus: [x,y,z],focus->active aperture_emit
% the position of scanline moves according to the change of parameters:focus, 
% so the intervals of scanlines(reflected by focus) should be set to the same value
% that's to say, the width of a pixel is  (distance of first and last scanline)/(the number of scanlines - 1)
% the width of the final scanning image is determined by the sum of   width of  every pixel of a row
% 
%**************************************************************************

tic;
    field_init(0)   % initial simulation system, need package of Field II

%% set basic parameters
    c=1540;             			% Speed of sound [m/s]
    f0=5e6;             			% Transducer center frequency [Hz]
    fs=40e6;                        % Sampling frequency [Hz]
    d_e_height=4/1000;      		% d_e_height of element [m]
    d_e_width=0.2798/1000;       	% d_e_width of element [m]
    d_e_kerf=0.025/1000;            % Distance between transducer elements [m]
    
    N_elements=128;     			% Number of elements
    N_active = 64;                  % Number of active elements
    N_receive = N_elements;         % the Number of receive elements which be wanted to set
    N_emit = N_active;              % Number of emit elements
    
    set_sampling(fs);                               % Set sampling frequency for simulation       
%     lambda=c/f0;        			% Wavelength [m]    
%     N_elements_sub = N_elements/2;  % Number of elements in a sub aperture when calculate MV algorithm
%     N_sub_aperture = N_elements-N_elements_sub+1;  % Number of sub aperture when calculate MV algorithm
%     dw_emit = N_active*d_e_width + (N_active-1)*d_e_kerf;
%     dw_receive = N_elements*d_e_width + (N_elements-1)*d_e_kerf;
    
%% load the computer phantom    
    N_elements=254;     			% Number of elements
    N_receive = 128;   % the Number of receive elements which be wanted to set
    N_active = 128;                  % Number of active elements
    d_x_size = 40/1000;   % d_e_width of phantom[m]
    d_z_size = 55/1000;   % d_e_height of phantom[m]
    d_z_start =5/1000;   % start of phantom surface[m]
    
   % for different phantom 
    if strcmp(Phantom_Name,'Point1') 
        x = 0/1000;
        y = 0;
        z = 25/1000;
    elseif strcmp(Phantom_Name,'Point2') 
        x = 0/1000;
        y = 0;
        z = 40/1000;
    elseif strcmp(Phantom_Name,'Point3') 
        x = 0/1000;
        y = 0;
        z = 55/1000;
    elseif strcmp(Phantom_Name,'Points3')
        x = [0/1000; 0/1000; 0/1000];
        y = [0; 0; 0];
        z = [25/1000; 40/1000; 55/1000];
    elseif strcmp(Phantom_Name,'Points6')
        N_elements=128;     			% Number of elements
        N_active = 64;                  % Number of active elements
        N_receive = N_elements;         % the Number of receive elements which be wanted to set
        N_emit = N_active;              % Number of emit elements        
        d_x_size = 10/1000;   % d_e_width of phantom[m]
        d_y_size = 4/1000;    % transverse d_e_width of phantom[m]
        d_z_size = 40/1000;   % d_e_height of phantom[m]
        d_z_start = 5/1000;   % start of phantom surface[m]
        % d_x_size = Phantom_Width;           % d_x_size = 10/1000;   % d_e_width of phantom[m]
        % d_z_size = Phantom_Depth;           % d_z_size = 40/1000;   % d_e_height of phantom[m]
        % d_z_start = Phantom_Depth_Start;    % d_z_start = 5/1000;   % start of phantom surface[m]
        % d_y_size = Phantom_Thickness;       % d_y_size = 4/1000;    % transverse d_e_width of phantom[m]

        x = [ -2; 2;  -0.4; 0.4;  -2; 2 ]/1000;
        y = zeros(6,1);
        z = [15; 15;  25; 25; 35; 35]/1000;    
    end
        
    phantom_positions = [x y z];
    phantom_amplitudes = ones(length(x),1)*7.8e25;
%     phantom_amplitudes(1:6,1) = 7.8e25;
    
%% Define the transducer
   % set aperture
    focus=[0 0 25]/1000;	ex = 1; ey = 1;    
%     aperture_emit = xdc_linear_array (N_active, d_e_width, d_e_height, d_e_kerf, ex, ey, focus);
    aperture_emit = xdc_linear_array (N_elements, d_e_width, d_e_height, d_e_kerf, ex, ey, focus);
%     aperture_receive = xdc_linear_array (N_elements, d_e_width, d_e_height, d_e_kerf, ex, ey, focus);
    aperture_receive = xdc_linear_array (N_receive, d_e_width, d_e_height, d_e_kerf, ex, ey, focus);
    % set the focus for the aperture_receive
    xdc_center_focus(aperture_receive,[0 0 0]);
    xdc_focus(aperture_receive,0,focus); 
  % Set receive time focus to NO delay, otherwise, the received data will be delayed.
% 	xdc_focus_times(aperture_receive, 0, zeros(1,N_elements));
	xdc_focus_times(aperture_receive, 0, zeros(1,N_receive));
    
  % set impulse response
    impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
    impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
    % xdc_impulse (aperture_emit, impulse_response);   
    % xdc_impulse (aperture_receive, impulse_response);
  % set excitation
    excitation=cos(2*pi*f0*(0:1/fs:1/f0));
    xdc_excitation (aperture_emit, excitation);

    
%% get simulation image data received
 N_lines=N_elements-N_active+1;  % Number of scan lines now is ...    
 image_data_t = zeros(N_lines,1);	% store the time for the first sample in scat
 for i=1:N_lines 
%     fprintf(' Scanning Lines: %d \n',i);
    % set the focus for this direction
    x_e_position = -( N_elements*d_e_width + (N_elements-1)*d_e_kerf )/2 ...
                    + (i-1)*(d_e_width+d_e_kerf) + ( N_active*d_e_width + (N_active-1)*d_e_kerf )/2;                
    xdc_center_focus(aperture_emit,[x_e_position 0 0]);
    xdc_focus(aperture_emit,0,[x_e_position 0 focus(3)]); 

    apo=[zeros(1, i-1) hamming(N_active)' zeros(1, N_elements-N_active-i+1)];
    xdc_apodization (aperture_emit, 0, apo);    
    
    % get simulation data
    [ v, image_data_t(i) ]=calc_scat_multi (aperture_emit, aperture_receive, phantom_positions, phantom_amplitudes);
%     if N_receive ~= N_elements
%         N_receive_start = round((N_elements-N_receive)/2);
%         v = v(:,N_receive_start:N_receive_start+N_receive-1);
%     end
    
    image_data_v(1:size(v,1),:,i) = v;
    
 end

 if N_receive ~= N_elements, N_active=2; end    % for calculating while simulation data was cut
 
 % save necessary data for calculate the image data
 save(['Z_SimulationData_',Phantom_Name],'image_data_v','image_data_t',...
            'd_x_size','d_z_size','d_z_start','d_e_width','d_e_kerf','N_active','c','fs','z');
%             'Phantom_Width','Phantom_Depth','Phantom_Depth_Start','d_e_width','d_e_kerf','N_active','c','fs','z');

fprintf('Simulated Data Completed. ');toc;fprintf('\n');

end

