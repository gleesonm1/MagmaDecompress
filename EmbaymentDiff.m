function [H2Od,CO2d,Sd,X]=EmbaymentDiff(PT,L,DR,Volatiles)
% Model to simulate diffusion of volatiles out of an embayment during magma
% ascent. Diffusivities for H2O and CO2 are taken from Zhang and Ni (2010)
% and S diffusivities are taken from Zhang et al. (2010).
%   Required Inputs: 
%           'PT' (n-by-2 matrix containing the pressure and temperature
%           (MPa and K)).
%           'L' (2-by-1 vector containing the length of the embayment and
%           the step size for the diffusion modelling, both in m)
%           'DR' (decompression rate in MPa/s)
%           'Volatiles' (n-by-3 matrix containing the H2O, CO2, and S
%           contents. H2O must be in wt%)
%           
%   Outputs:
%           'H2Od' (H2O profile at the end of calculation)
%           'CO2d' (CO2 profile at the end of calculation)
%           'Sd' (S profile at the end of calculation)
%           'X' (location along profile)
%
%   This code is designed to simulate the diffusive loss of H2O, CO2 or S
%   from an embayment during magma ascent. I have purposefully not included
%   calculation of the degassing systematics in this code so that the user
%   can define their own degassing paths. If you wish to calculate the
%   degassing paths of H2O and CO2 using the MagmaSat module then please
%   see my Decompress code.

    % find stable time step
    DH2O_max=max(Volatiles(1,1).*exp(-8.56-19110./PT(:,2)));
    DCO2_max=max(exp(-13.99-(17367+1944.8*(PT(1,1)/1000))./PT(:,2)+Volatiles(1,1).*(855.2+271.2*(PT(1,1)/1000))./PT(:,2)));
    DS_max=max(exp(-8.21-(27692-651.6*Volatiles(1,1))./PT(:,2)));
    maxD=max([DH2O_max DCO2_max DS_max]);
    
    R=0.49;
    dt=(R*(L(2)^2))/maxD; %time step in seconds
    
    %pressure step in calculations
    dp=DR*dt;
    
    % create pressure array for calculations
    P_diff=linspace(PT(1,1),PT(end,1),(PT(1,1)-PT(end,1))/dp+1);
    
    % re-calculate the volatile content of the melt phase
    Temp=csapi(PT(:,1),PT(:,2),P_diff);
    S_diff=csapi(PT(:,1),Volatiles(:,3),P_diff);
    H2O_diff=csapi(PT(:,1),Volatiles(:,1),P_diff);
    CO2_diff=csapi(PT(:,1),Volatiles(:,2),P_diff);
    
    S_diff(S_diff<0)=0; S_diff(isnan(S_diff))=max(S_diff);
    H2O_diff(H2O_diff<0)=0;
    CO2_diff(CO2_diff<0)=0;

    % create array of x values
    X=linspace(0,L(1),L(1)/L(2)+1);
    
    % set initial concentration profiles
    Sd=zeros(1,length(X))+S_diff(1); Sn=zeros(1,length(X));
    H2Od=zeros(1,length(X))+H2O_diff(1); H2On=zeros(1,length(X));
    CO2d=zeros(1,length(X))+CO2_diff(1); CO2n=zeros(1,length(X));
    
    D_S=zeros(1,length(X));
    D_H2O=zeros(1,length(X));
    D_CO2=zeros(1,length(X));
    
    %% actually do the diffusion models
    for i=1:length(P_diff)
        Sd(1)=S_diff(i);
        H2Od(1)=H2O_diff(i);
        CO2d(1)=CO2_diff(i);
        D_H2O=H2Od.*exp(-8.56-19110./Temp(i));
        D_CO2=exp(-13.99-(17367+1944.8*(P_diff(i)/1000))./Temp(i)+H2Od.*(855.2+271.2*(P_diff(i)/1000))./Temp(i));
        D_S=exp(-8.21-(27692-651.6*H2Od)./Temp(i));
        for j=2:length(X)-1
            Sn(j)=Sd(j)+((dt*D_S(j))/(L(2)^2))*(Sd(j-1)+Sd(j+1)-2*Sd(j));
            H2On(j)=H2Od(j)+((dt*D_H2O(j))/(L(2)^2))*(H2Od(j-1)+H2Od(j+1)-2*H2Od(j));
            CO2n(j)=CO2d(j)+((dt*D_CO2(j))/(L(2)^2))*(CO2d(j-1)+CO2d(j+1)-2*CO2d(j));
        end
        Sn(1)=Sd(1);
        H2On(1)=H2Od(1);
        CO2n(1)=CO2d(1);
        Sn(end)=Sd(end)+((dt*D_S(end))/L(2))*(Sd(end-1)-Sd(end));
        H2On(end)=H2Od(end)+((dt*D_H2O(end))/L(2))*(H2Od(end-1)-H2Od(end));
        CO2n(end)=CO2d(end)+((dt*D_CO2(end))/L(2))*(CO2d(end-1)-CO2d(end));

        Sd=Sn;
        H2Od=H2On;
        CO2d=CO2n;
    end

end