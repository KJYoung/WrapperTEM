function [PartPot, pout] = loadSamples(params2,rpdbsave,saved)
%loadSamples Loads or calculates the interaction potential of the particles
%depending on the required input params2.spec.source : 'pdb', 'map' or 'amorph'
% SYNOPSIS:
% [PartPot,pout] = loadSamples(params2)
%
% PARAMETERS:
%  params2: structure containing various input paramters
%
% OUTPUT:
%  PartPot: Interaction potential of the particle 
%  pout   : possibly changed parameters (number of particles)

% (C) Copyright 2013
%  Quantitative Imaging Group      Leiden University Medical Center
%  Faculty of Applied Sciences     Department of Molecular Cell Biology
%  Delft University of Technology  Section Electron Microscopy
%  Lorentzweg 1                    2300 RC Leiden
%  2628 CJ Delft
%  The Netherlands
%
%  Milos Vulovic

       
pout = params2;
dir0 = params2.proc.rawdir;

if strcmp(params2.spec.source, 'pdb')   % generate the potential maps from pdb
    % first check if the particles are already generated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if params2.spec.imagpot~=3
        list = dir([dir0 filesep 'Particles' filesep params2.spec.pdbin '*MF' sprintf('%3.1f',params2.spec.motblur) '_VoxSize' sprintf('%02.2f',params2.acquis.pixsize*1e10) '*A.raw']);
    else
        list = dir([dir0 filesep 'Particles' filesep params2.spec.pdbin '*MF' sprintf('%3.1f',params2.spec.motblur) '_VoxSize' sprintf('%02.2f',params2.acquis.pixsize*1e10) '*A_Volt' sprintf('%03d',params2.acquis.Voltage/1000) 'kV.raw']);
    end
   params2.NumGenPart = size(list,1);
   if params2.proc.partNum <= params2.NumGenPart %&& ~params2.proc.geom 
         PartPot = 0;
         fprintf('No need for generation of new particles. Load the potential maps from the folder ''Particles''...\n')
   else         
        if params2.proc.geom 
         % load the specified orientation of the particles
         [xt, yt, zt, alphad, betad, gammad] = PartList(params2);
        else
         % random uniform orientation roll-pitch-yaw euler angles
        alpha = 2*pi*rand(1, params2.proc.partNum); 
        beta  = acos(1-2*rand(1, params2.proc.partNum));
        gamma = 2*pi*rand(1, params2.proc.partNum);
        alphad = rad2deg(alpha); betad = rad2deg(beta); gammad = rad2deg(gamma);
        end  
        % calculate the atomic potential from the pdb with given
        % orientation. Save the particles into subfolder 'Particles' (wr=1)
      wr = 1;
      [atompotF, numpart] = AtomPotRot(params2, alphad, betad, gammad, wr,rpdbsave,saved);
      pout.proc.partNum = numpart + params2.NumGenPart;
      PartPot = atompotF;
   end    
else
    error('This type of the input type is not known. Options are ''pdb'' ''map'' or '' amorphous''')
end  



