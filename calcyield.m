function [yield, Cgas, Cpar] = calcyield(x,voc,vbs,env)

% Get 'x'
vbs.alpha = x(1:vbs.nbins); %  mass yields, unitless
% vbs.alpha = vbs.alpha / sum(vbs.alpha); % normalize so that total mass = 1
vbs.kaging = 10^(x(vbs.nbins+1)); % aging reaction rate constant, cm3/molecules/s
vbs.mfrag = 10^x(vbs.nbins+2); % probability of fragmentation, unitless

% Calculate probability of fragmentation
for j=1:vbs.nbins
    pfrag(j) = ( (log10(vbs.cstar(end))-log10(vbs.cstar(j))) / (log10(vbs.cstar(end))-log10(vbs.cstar(1))) )^vbs.mfrag; % Changes between 1 (only fragments) to 0 (only functionalizes) as a function of C*, mfrag determines the curvature
    %         pfrag(j) = 1 - exp(vbs.mfrag * (log10(vbs.cstar(j))-log10(vbs.cstarspecies)) / log10(vbs.cstarspecies)); % Based on simpleSOM formulation, Chris Cappa
    %     pfrag(j) = vbs.mfrag;
end

% Base aging kernel, moves things 1 down and 1 up
ak = zeros(vbs.nbins,vbs.nbins);
for j=1:vbs.nbins
    if j==1
        ak(j,j) = 1-pfrag(j)-1;
        ak(j,j+1) = pfrag(j);
    elseif j==vbs.nbins
        ak(j,j-1) = 1-pfrag(j);
        ak(j,j) = pfrag(j)-1;
    else
        ak(j,j-1) = 1-pfrag(j);
        ak(j,j+1) = pfrag(j);
        ak(j,j) = -1;
    end
end

% % Advanced aging kernel, moves things 2 down and 1 up
% ak = zeros(vbs.nbins,vbs.nbins);
% for j=1:vbs.nbins
%     if j==1
%         ak(j,j) = 1-pfrag(j)-1;
%         ak(j,j+1) = pfrag(j);
%     elseif j==2
%         ak(j,j-1) = 1-pfrag(j);
%         ak(j,j+1) = pfrag(j);
%         ak(j,j) = -1;
%     elseif j==vbs.nbins
%         ak(j,j-1) = 1-pfrag(j);
%         ak(j,j) = pfrag(j)-1;
%     else
%         ak(j,j-2) = 1-pfrag(j);
%         ak(j,j+1) = pfrag(j);
%         ak(j,j) = -1;
%     end
% end

% Simulation
for m=1:length(env.COA)
    
    ntime = length(env.time); % number of time intervals
    
    % Initialize time=0
    VOC(1) = voc.initial;
    rVOC = 0;
    Ctot(1,vbs.nbins) = 0;
    Cgas = Ctot;
    Cpar = Ctot;
    yield(m,1) = 0;
    
    for i=2:ntime
        
        % time interval
        dt = env.time(i) - env.time(i-1);
        
        % Copy Cgas, Cpar, and Ctot
        Cpar(i,:) = Cpar(i-1,:);
        Ctot(i,:) = Ctot(i-1,:);
        
        % ad-hoc aging (no frag)
        for j=1:vbs.nbins
            dCgas(j) = Cgas(i-1,j) * (1-exp(-vbs.kaging*env.OHconc*dt));
            for k=1:vbs.nbins
                daging(j,k) = dCgas(j) * ak(j,k);
            end
        end
        for j=1:vbs.nbins
            Cgas(i,j) = Cgas(i-1,j) + sum(daging(:,j));
        end
        
        % VOC reacted
        VOC(i) = VOC(i-1) * exp(-voc.koh*env.OHconc*dt);
        dVOC = (VOC(i-1) - VOC(i)) * (env.Pres/env.Rgas/env.Temp*voc.mw);
        rVOC(i) = rVOC(i-1) + dVOC;
        
        % add organic mass to bins
        Ctot(i,:) = Cgas(i,:) + Cpar(i,:) + dVOC * vbs.alpha;
        
        % g/p partitioning
        for j=1:vbs.nbins
            Cpar(i,j) = (1+vbs.cstar(j)/env.COA(m))^(-1) * Ctot(i,j);
            Cgas(i,j) = Ctot(i,j) - Cpar(i,j);
        end
        
        % calculate SOA yield
        yield(m,i) = sum(Cpar(i,:)) / rVOC(i);
        
    end
    
end

