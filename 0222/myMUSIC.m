function doas = myMUSIC(Rx,sample_points,peak_num)
    sz = size(Rx);
    M = sz(1);
    G=sample_points/180;
    phase_num = 180*G;

    theta_plot = linspace(-pi/2,pi/2,phase_num);
    a_plot = zeros(M,phase_num);
    for i=1:M
        for j=1:phase_num
            a_plot(i,j) = exp(-1j*pi*(i-1)*sin(theta_plot(j)));
        end
    end     
    
    [V,D]=eig(Rx);

    E = zeros(M,M-peak_num);
    for i=1:M-peak_num
        E(:,i) = V(:,i+peak_num);
    end
    P_MUSIC = zeros(1,phase_num);
    for i=1:phase_num
        P_MUSIC(i) = abs(1/(a_plot(:,i)'*E*(E')*a_plot(:,i)));
    end
    doas=zeros(1,peak_num);
    [peak_list,loc_peak] = findpeaks(P_MUSIC);
    for i=1:peak_num
        [a,loc_max] = max(peak_list);
        doas(i) = (loc_peak(loc_max)/G)-90;
        peak_list(loc_max)=0;
    end
end

