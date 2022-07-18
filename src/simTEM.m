function imStructOut= simTEM(InputVol, params2)
    %simTEM simulates the full image formation given the input interaction
    % potential and acquisition parameters (dose, optics and detector) 
    %
    % SYNOPSIS:
    % [imStructOut,paramsout]= simTEM(InputVol, params2)
    %
    % PARAMETERS:
    %  InputVol: Interaction potenital specimen volume 
    %   params2: Structure containing various input simulation paramters
    %
    % OUTPUT:
    % imStructOut: Structure containg output images (or stacks): noisy, noiseless and exit wave image 
    %              Optional: ctf and mtf (if params.disp.ctf or params.disp.mtfdqe are set)
    
    % (C) Copyright 2013
    %  Quantitative Imaging Group      Leiden University Medical Center
    %  Faculty of Applied Sciences     Department of Molecular Cell Biology
    %  Delft University of Technology  Section Electron Microscopy
    %  Lorentzweg 1                    2300 RC Leiden
    %  2628 CJ Delft
    %  The Netherlands
    %
    %  Milos Vulovic
    %
    voxSz = params2.acquis.pixsize;% the voxel size
    % prealocate memory for the stack (series)
    nTiltAngles = length(params2.acquis.tilt);
    
    if strcmp(params2.seriesout,'tilt')
        Nseries= nTiltAngles;
    elseif strcmp(params2.seriesout,'defocus')
        Nseries= length(params2.acquis.df);
    elseif strcmp(params2.seriesout,'dose')
        Nseries= length(params2.acquis.dose_on_sample);
    else 
        Nseries=1;
    end
    btot_i        = newim(params2.proc.N,params2.proc.N,Nseries);
    
    poten0 = InputVol; clear InputVol;
    thickness = voxSz*size(poten0,3);
    
    switch params2.inter.type
        case 'ms' % if multislice then the exit wave is calculated separately for each tilt angle.  
            psi_exit = newim(params2.proc.N, params2.proc.N, nTiltAngles, 'dcomplex');
            %if strcmp(params2.seriesout,'tilt')
            for ll=1:nTiltAngles
                
                tiltang = params2.acquis.tilt(ll);
                fprintf('Simulate tilt angle %3.0f\n', tiltang*180/pi)
                disp(' ')
                disp('Generating Tilt...')
                tic
                [potenext, n] = tiltingMS(poten0,tiltang, params2);
                clear poten0
                toc
                disp(' ')
                sizepot = size(potenext);
                Nm      = max(sizepot(1), sizepot(2));
                thicknessfull = sizepot(3)*voxSz;
                potenfull = potenext;
                clear potenext
                
                %Fourier domain
                xwm = (voxSz)*Nm;%pixelsize for multislice * size sample
                q_true_pix_m = 1/xwm;
                q_m = rr([Nm Nm])*q_true_pix_m; % frequencies in Fourier domain
        
                % propagator function Fresnel Propagation
                dzprop = thicknessfull/n;
                
                % MULTISLICE
                psi_slice = newim(Nm,Nm,'scomplex')+1;
        
                % Phase grating & Multislice at  same time
                disp('-.-.-.-.-.-.-.-.- STARTING MULTISLICE -.-.-.-.-.-.-.-.-...')
                tic
                for ii = 0:n-1
                    vai = tic;
                        disp(['SLICE ' char(num2str(ii+1)) '/' char(num2str(n))])
                        disp('     Phase grating...')
                        tic
                        % Projected potential in a slice 
                        t = mean(potenfull(:,:,ii*(sizepot(3)/n):(ii+1)*(sizepot(3)/n)-1),[],3); 
                        toc
        
                        disp('     Transmission function...')
                        tic
                        % transmission functions for a single slice as we go through
                        psi_t = exp(1i*params2.inter.sig_transfer*t*dzprop); 
                        toc

                        disp('     Fresnel propagator...')
                        tic
                        % Fresnel propagator
                        P = exp(-1i*pi*params2.inter.lambda*(q_m.^2)*dzprop); % Fresnel propagator
                        toc
        
                        disp('     Wave propagation...')
                        tic
                        % Exit wave from current slice
                        psi_slice = ift(ft(psi_slice*squeeze(psi_t(:,:,0)))*P);
                        toc
                    disp(['FINISHED SLICE ' char(num2str(ii+1)) '/' char(num2str(n))])
                    toc(vai)
                end
                PsiExit = psi_slice;
                psi_exit(:,:,ll-1) = PsiExit;
                disp(' ')
                disp('-.-.-.-.-.-.-.-.- FINISHED MULTISLICE -.-.-.-.-.-.-.-.-')
                toc
        
                if strcmp(params2.spec.source, 'amorph')
                    thicknessfull = params2.spec.thick/cos(tiltang);
                    psi_exit(:,:,ll-1) = psi_exit(:,:,ll-1)*exp(-params2.inter.sig_transfer*params2.spec.potenampl*thicknessfull);
                end
            end
    end
    if strcmp(params2.seriesout,'defocus') || strcmp(params2.seriesout,'dose')
       psi_exit = repmat(psi_exit,[1 1 Nseries]);
    end
    
    %% ---------------------------------- CTF with df, ast, envelopes and optionally phase plate
    
    psi_exit=double(psi_exit);   
    btot_i=double(btot_i);
    for jjj= 1:Nseries
        if strcmp(params2.seriesout,'defocus')
            params2.acquis.df_run=params2.acquis.df(jjj);
        else
            params2.acquis.df_run=params2.acquis.df(1);
        end
        if ~mod(jjj,5)||~mod(jjj,Nseries)
            fprintf(['Calculate the CTF for the ' params2.seriesout sprintf(' series. Image number %3d of %3d\n',  jjj, Nseries)]);
        end
      
      [ctf] = simulateCTF(params2);
      btot = ctf*dip_fouriertransform(dip_image(psi_exit(:,:,jjj)),'forward',[1 1 ]);
      btot_i(:,:,jjj) = double(abs(dip_fouriertransform(btot,'inverse',[1 1])).^2); % intensity in the image without camera influence              
    end
    noiseless_tilt_series = dip_image(btot_i); 
    
    %%  --------------------------------- Camera influence          
    noiseless_tilt_series=double(noiseless_tilt_series);
    imStructOut.noiseless_series = noiseless_tilt_series;