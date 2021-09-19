function [ptpath]=Decompress(Conditions,bulk,Closed,Thermal,H2O,CO2)
% Magma decompression model to simulate volatile changes during isentropic
% magma ascent. This code uses rhyoliteMELTS v1.2.0 (Gualda and Ghiorso,
% 2015). Please ensure that you have added the relevant alphaMELTS for
% MATLAB folders to your path. This code uses the alphaMELTS for MATLAB
% functions as of Dec 2020.
%   Required Inputs: 
%           'Conditions' (2-by-1 vector containing the initial pressure and 
%           final pressure in MPa)
%           'Bulk' (19-by-1 vector with initial liquid composition)
%           'Closed' (true or false, if true simulated closed system
%           degassing)
%           'Thermal' (interger - 1 or 2, if 1 use isothermal constraints, 
%           if 2, use isenthalpic constraints)
%   Optional Inputs:
%           'H2O' (set H2O content of initial liquid if not present in bulk)
%           'CO2' (set CO2 content of initial liquid if not present in bulk)
%           
%   Outputs:
%           'ptpath' (MELTSdynamic for the ascent path)
    
    % set initial conditions for model
    pressure=Conditions(1)*10; %convert from MPa to bars, initial pressure
    temperature=1400; % set arbitrary T to start calculation.
    
    FPbar=Conditions(2)*10; % convert from MPa to bars, final pressure
    
    if nargin>4
        bulk(15)=H2O;
    else
        H2O=bulk(15);
        CO2=bulk(16);
    end
    
    bulk(16)=0;
    
    bulk(1:14)=bulk(1:14)./sum(bulk(1:14)).*(100-H2O-CO2);
    
    % load melts dynamic
    ptpath = MELTSdynamic(4);
    
    ptpath.engine.set('bulkComposition', bulk);
    ptpath.engine.pressure = pressure;
    ptpath.engine.temperature = temperature;
    
    if Closed~=true
        ptpath.engine.setSystemProperties("Mode", "Fractionate Fluids");
    end
    
    %ptpath.engine.setSystemProperties("Mode", "Isentropic");
    
    ptpath.engine.findLiquidus;
    
    bulk(16)=CO2;
    ptpath.engine.set('bulkComposition', bulk);
    
    %% carry out decompression model
    while ptpath.engine.pressure > FPbar
        ptpath = ptpath.addNodeAfter;
        if ptpath.engine.pressure>=100*FPbar
            ptpath.engine.pressure = ptpath.engine.pressure - 50;
        end
        if ptpath.engine.pressure>=10*FPbar
            ptpath.engine.pressure = ptpath.engine.pressure - 10;
        end
        if ptpath.engine.pressure<10*FPbar
            ptpath.engine.pressure = ptpath.engine.pressure - 0.1;
        end 
        ptpath.engine.calcEquilibriumState(Thermal, 1);       
    end
end