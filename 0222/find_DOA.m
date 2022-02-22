function [DOA_s,DOA_i] = find_DOA(matX,sample_angle)
    sz = size(matX);
    M = sz(1);
    L = sz(2);
    G=sample_angle/180;
    phase_num = 180*G;

    theta_plot = linspace(-pi/2,pi/2,phase_num);
    a_plot = zeros(M,phase_num);
    for i=1:M
        for j=1:phase_num
            a_plot(i,j) = exp(-1j*pi*(i-1)*sin(theta_plot(j)));
        end
    end     
    
    lambda = 0.9;
    delta = 3*10^(-9);
    range = 44;
    DOA_s = zeros(1,L);
    DOA_i = zeros(1,L);
    for L1=1:L
        Rx = delta*eye(M);
        if L1<=range
            for L2=1:L1
                Rx = Rx+lambda^(L1-L2)*((matX(:,L2))*matX(:,L2)');
            end
        else
            for L2=(L1-range):L1
                Rx = Rx+lambda^(L1-L2)*((matX(:,L2))*matX(:,L2)');
            end
        end

        P_MVDR = zeros(1,phase_num);
        for i=1:phase_num
            P_MVDR(i) = 1/(a_plot(:,i)'*(Rx\a_plot(:,i)));
        end
        beampattern = abs(P_MVDR);

        [pks,loc]=findpeaks(beampattern);

        [P1, loc1] = max(pks);
        pks(loc1) = 0;
        [P2, loc2] = max(pks);
        if (loc(loc1)>=70*G) && (loc(loc1)<=90*G)
            DOA_s(L1) = (loc(loc1)-90*G)/G;
            DOA_i(L1) = (loc(loc2)-90*G)/G;
        elseif (loc(loc2)>=70*G) && (loc(loc2)<=90*G)
            DOA_s(L1) = (loc(loc2)-90*G)/G;
            DOA_i(L1) = (loc(loc1)-90*G)/G;
        else
            flag=0;
            DOA_i(L1) = (loc(loc1)-90*G)/G;
            while(flag<10)
                pks(loc2) = 0;
                [P2, loc2] = max(pks);
                if (loc(loc2)>=70*G) && (loc(loc2)<=90*G)
                    DOA_s(L1) = (loc(loc2)-90*G)/G;
                    flag=10;
                else
                    DOA_s(L1) = (loc(loc2)-90*G)/G;
                    flag=flag+1;
                end
            end
        end
    end
    
    figure()
    subplot(2,1,1)
    plot(DOA_i)
    xlabel('t')
    ylabel('degree')
    ylim([-90 90])
    title('DOA of interfering signal')
    subplot(2,1,2)
    plot(DOA_s)
    xlabel('t')
    ylabel('degree')
    ylim([-90 90])
    title('DOA of source signal')
end

