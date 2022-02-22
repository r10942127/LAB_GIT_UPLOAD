function [y,w] = LCMV_opt(C,g,x,M,J,Fs,freq,d)
    IT = floor(length(x(1,:))/Fs);
    Rx = zeros(J*M,J*M);
    max_freq = max(freq);
    for i=1:IT*Fs-J+1
        xx=[];
        for j=1:J
            if i-(j-1)>0
                xx = [xx;x(:,i-(j-1))];
            else
                xx = [xx;zeros(M,1)];
            end
        end
        Rx = (Rx*(i-1) + xx*(xx'))/i;
    end
    w = (Rx\C)*((C'*(Rx\C))\g);
    
    y = zeros(1,IT*Fs);
    for it=1:IT
        for j=1:J
            for n=1:Fs
                if (it-1)*Fs+n-(j-1) > 0
                    y((it-1)*Fs+n) = y((it-1)*Fs+n) + x(:,(it-1)*Fs+n-(j-1)).'*conj(w((M*(j-1)+1):(M*j)));
                end
            end
        end
    end
    
    plot_num=180;
    theta_plot = linspace(-pi/2,pi/2,plot_num);
    a_plot = zeros(J*M,plot_num,length(freq));
    for k=1:length(freq)
        for j=1:J
            for m=1:M
                for num=1:plot_num
                    if d==1
                        a_plot((j-1)*M+m,num,k) = exp(-1j*2*pi*(freq(k)/Fs)*((m-1)*sin(theta_plot(num))+(j-1)));
                    elseif d==2
                        a_plot(M*(j-1)+m,num,k) = exp(-1j*pi*(freq(k)/max_freq)*((m-1)*sin(theta_plot(num))+(2*max_freq*(j-1)/Fs)));
                    end
                end
            end 
        end
    end
    B_LCMV = zeros(length(freq),plot_num);
    for k=1:length(freq)
        for i=1:plot_num
            B_LCMV(k,i) = (w')*a_plot(:,i,k);
        end
    end
    
    % plot
    figure()
    subplot(2,1,1)
    F = abs(fft(y));
    plot(linspace(0,Fs,IT*Fs),F)
    title('Frequency response of output y(t)')
    
    subplot(2,1,2)
    for k=1:length(freq)
        B_LCMV_db = 20*log10(abs(abs(B_LCMV(k,:))));
        plot(theta_plot,B_LCMV_db)
        hold on
    end
    title('Beampattern')
    xlim([-pi/2 pi/2])
    xticks([-90*pi/180 -80*pi/180 -70*pi/180 -60*pi/180 -50*pi/180 -40*pi/180 -30*pi/180 -20*pi/180 -10*pi/180 0 10*pi/180 20*pi/180 30*pi/180 40*pi/180 50*pi/180 60*pi/180 70*pi/180 80*pi/180 90*pi/180])
    xticklabels({'$-90^{\circ}$','$-80^{\circ}$','$-70^{\circ}$','$-60^{\circ}$','$-50^{\circ}$','$-40^{\circ}$','$-30^{\circ}$','$-20^{\circ}$','$-10^{\circ}$','0','$10^{\circ}$','$20^{\circ}$','$30^{\circ}$','$40^{\circ}$','$50^{\circ}$','$60^{\circ}$','$70^{\circ}$','$80^{\circ}$','$90^{\circ}$'});
    set(gca,'TickLabelInterpreter', 'latex')
    hold off
    sgtitle('LCMV Optimization') 
end

