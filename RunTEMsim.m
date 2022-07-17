%%  RunTEMsim % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%                                                                         %
%     Script to simulate micrographs under different conditions:          %
%       - with/without phase plate                                        %
%       - defocus range                                                   %
%       - motion blur (radiation damage)                                  %
%       - correction factor of the microscope                             %
%       - electron dose                                                   %
%       - size of the micrograph (pixel number)                           %
%       - pixel size of the micrograph                                    %
%       - distance between particles                                      %  
%                                                                         %
%     To generate the micrographs the script calls the function TEMsim,   %
%       which is an adaptation to the code found in Vulovic, 2013         %
%                                                                         %
%     The main changes to the original code include:                      %
%       -Addition of PHASE PLATE shift                                    %
%       -Optimization of memory usage for parallel processing of          %
%         micrographs                                                     %
%                                                                         %
%    IMPORTANT                                                            %
%       The simulator is based on the DIPimage Matlab Toolbox.            %
%       A working version of DIPimage must be installed on this computer  %
%       in order for the simulation to function properly.                 %
%       Download links can be found at the following address:             %
%           http://www.diplib.org/                                        %
%                                                                         %
%                                                                         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
addpath('/dip/common/dipimage')
addpath('./src')
dip_initialise

clc;
close all;
clear;

time = tic;

% Default Parameters
% mg = 1;                     % Number of micrographs to generate
% pp = 1;                     % Phase plate (0 = no; 1 = yes)
% df_range= [300 800];        % Defocus range [nm]
% mb_series = [0.5];          % Motion blur series (If multiple MB, enter them as a vector)
% cf_series = [14];           % Correction factor  (If multiple CF, enter them as a vector)
% dose      = [50];           % Electron dose to the specimen [e-/A2]
% pix = 512;                  % Number of pixels 
% pixsize = 0.5;              % Pixel size [A]
% mindist = 50/pixsize;       % Minimum distance between particles divided by pixsize, 
                             % Depends on type of protein (apo ferrtin ~150)
% dir1 =  './Micrographs'; % Select folder where to save micrographs
% NT2C Parameters
mg          = 250;                % Number of micrographs to generate
pp          = 0;                % (Default)Phase plate (0 = no; 1 = yes)
% df_range    = [2500 3200];      % Defocus range [nm]
df_range    = [2800 3800];      % Defocus range [nm]
mb_series   = [0];            % (Default)Motion blur series (If multiple MB, enter them as a vector)
cf_series   = [0];             % (Default)Correction factor  (If multiple CF, enter them as a vector)
dose        = [20];             % Electron dose to the specimen [e-/A2]
pix         = 4096;             % Number of pixels 
pixsize     = 1.34;             % Pixel size [A]
mindist     = 150/pixsize;       % (Default = 50)Minimum distance between particles divided by pixsize, 
                                % Depends on type of protein (apo ferrtin ~150)
dir_noiseless =  './Micrographs/testClean5lzf'; % Select folder where to save micrographs


tot = mg*length(mb_series)*length(cf_series)*length(dose);
count  = 0;
rng('shuffle');

disp(...
   [char("###################### Starting new Simulation ######################") newline...
    char("       Micrographs: " + tot) newline...
    char("       Phase Plate: " + pp) newline...
    char("       Size:        " + pix + "x" + pix) newline...
    char("       Particles:   " + "Randomized") newline...
    char("       Defocus:     [" + df_range(1) + " - " + df_range(2)+"]nm") newline...
    char("       Motion blur: " + mb_series) newline...
    char("       Corr factor: " + cf_series) newline...
    char("       Dose:        " + dose) 'e/A^2' newline...
    '                                 ' char(datetime(now,'ConvertFrom','datenum'))])

% Read PDB File
tic
 disp([newline 'Loading PDB file...'])
 % rpdb = pdbread('./PDBs/2wrj.pdb'); % Select PDB file to read
 rpdb = pdbread('./PDBs/5lzf.pdb'); % Select PDB file to read
toc
disp(' ')               
                
for micro = 1:mg
    for cf = 1 : length(cf_series)
        for d = 1:length(dose)
            for mb = 1 : length (mb_series)    
                try
                    count = count + 1; 
                  
                    ou = tic;
                    
                    % Randomize particles positions for each micrograph
                    disp('Generating random particle positions...')
                    tic
                    circles = Randomposition(pix, mindist);
                    particles = length(circles);
                    toc
                    
                    % Randomize defocus for each micrograph
                    defocus = (df_range(2)-df_range(1))*rand(1,1)+df_range(1);

                    disp([newline newline newline])                  
                    disp( "-----------------------   Progress:   "+count+" / "+tot+" Micrographs   -----------------------")
                    disp(' ')
                    disp(char( "Particles:    "+particles) )
                    disp([char( "Defocus:      "+defocus) 'nm' newline])
                    disp(' ')

                    % Run simulation to generate micrograph
                    out = TEMsim(defocus,mb_series(mb),cf_series(cf),1,rpdb,particles,dose(d),circles,pix,pixsize,pp,mindist);
                    ImageNoiseless = out.noiseless_series;
                    delete ./Raw/Particles/*.raw; % delete particle positions for the generated micrographs (for memory) */
                    
                    % Write micrograph to .MRC file
                    format_time = 'mm_dd_HH';
                    s1 = "D4_1_"+ (micro) +"_p"+particles+"_df"+round(defocus)+"-"+datestr(now, format_time)+".mrc";
                    WriteMRC(double(ImageNoiseless),0.5,[dir_noiseless filesep char(s1)]);
                    disp(' ')
                    disp('Successful Micrograph')
                    
                catch ex
                      
                    % In case the simulator fails to produce a micrograph,
                    % it will skip it and display a warning
                      warning(['Failed Micrograph' newline ex.identifier  ...
                                    newline ex.message      ...
                                    newline ex.stack.name])
                      disp(' ')   
                end
            end
        end
    end
end

disp(' ')
disp(...
   [char("###################### End of Simulation ######################") newline newline...
    '                     ' char(datetime(now,'ConvertFrom','datenum')) newline])
toc(time)

