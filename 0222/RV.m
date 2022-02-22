function [y,w] = RV(Fs,x,source_freq,interfere_freq,angle,M,d)
    % hyperparameter
    beta = 10;
    sample = 1000;
    max_freq = max(max(source_freq),max(interfere_freq));
    J = input('´X­ÓTap? ');
    a_source = zeros(J*M,1);
    a_interfere = zeros(J*M,1);
    for j=1:J
        for m=1:M
            if d==1
                a_source(M*(j-1)+m,1) = exp(-1j*2*pi*(source_freq(1)/Fs)*((m-1)*sin(angle(1))+(j-1)));
                a_interfere(M*(j-1)+m,1) = exp(-1j*2*pi*(interfere_freq(1)/Fs)*((m-1)*sin(angle(2))+(j-1)));
            elseif d==2
                a_source(M*(j-1)+m,1) = exp(-1j*pi*(source_freq(1)/max_freq)*((m-1)*sin(angle(1))+(2*max_freq*(j-1)/Fs)));
                a_interfere(M*(j-1)+m,1) = exp(-1j*pi*(interfere_freq(1)/max_freq)*((m-1)*sin(angle(2))+(2*max_freq*(j-1)/Fs)));
            end
        end
    end
    C = [real(a_source) imag(a_source) real(a_interfere) imag(a_interfere)];
    g = [1;0;0;0];
    IT = floor(length(x(1,:))/Fs);
    Rx = zeros(J*M,J*M);
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
    interest = [min([source_freq interfere_freq]) max([source_freq interfere_freq])];
    range = interest(2)-interest(1);
    Q0 = zeros(J*M,J*M);
    for s=1:sample
        s_steer = zeros(J*M,1);
        i_steer = zeros(J*M,1);
        for j=1:J
            for m=1:M
                if d==1
                    s_steer(M*(j-1)+m) = exp(-1j*2*pi*((interest(1)+range*(s/sample))/Fs)*((m-1)*sin(angle(1))+(j-1)));
                    i_steer(M*(j-1)+m) = exp(-1j*2*pi*((interest(1)+range*(s/sample))/Fs)*((m-1)*sin(angle(2))+(j-1)));
                elseif d==2
                    s_steer(M*(j-1)+m) = exp(-1j*pi*((interest(1)+range*(s/sample))/max_freq)*((m-1)*sin(angle(1))+(2*max_freq*(j-1)/Fs)));
                    i_steer(M*(j-1)+m) = exp(-1j*pi*((interest(1)+range*(s/sample))/max_freq)*((m-1)*sin(angle(2))+(2*max_freq*(j-1)/Fs)));
                end
            end
        end
        dis = (s_steer-a_source)*(s_steer-a_source)';
        dis = dis + (i_steer-a_interfere)*(i_steer-a_interfere)';
        Q0 = (Q0*(s-1) + dis)/s;
    end
    P = Rx + beta*Q0;
    w = (P\C)*((C'*(P\C))\g);
    
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
    a_plot = zeros(J*M,plot_num,sample);
    for k=1:sample
        for j=1:J
            for m=1:M
                for num=1:plot_num
                    if d==1
                        a_plot((j-1)*M+m,num,k) = exp(-1j*2*pi*((interest(1)+(range*(k/sample)))/Fs)*((m-1)*sin(theta_plot(num))+(j-1)));
                    elseif d==2
                        a_plot(M*(j-1)+m,num,k) = exp(-1j*pi*((interest(1)+(range*(k/sample)))/max_freq)*((m-1)*sin(theta_plot(num))+(2*max_freq*(j-1)/Fs)));
                    end
                end
            end 
        end
    end
    B_RV = zeros(sample,plot_num);
    for k=1:sample
        for i=1:plot_num
            B_RV(k,i) = (w')*a_plot(:,i,k);
        end
    end
    
    % plot
    figure()
    subplot(2,1,1)
    F = abs(fft(y));
    plot(linspace(0,Fs,IT*Fs),F)
    title('Frequency response of output y(t)')
    
    subplot(2,1,2)
    for k=1:sample
        B_RV_db = 20*log10(abs(abs(B_RV(k,:))));
        plot(theta_plot,B_RV_db)
        hold on
    end
    title('Beampattern')
    xlim([-pi/2 pi/2])
    xticks([-90*pi/180 -80*pi/180 -70*pi/180 -60*pi/180 -50*pi/180 -40*pi/180 -30*pi/180 -20*pi/180 -10*pi/180 0 10*pi/180 20*pi/180 30*pi/180 40*pi/180 50*pi/180 60*pi/180 70*pi/180 80*pi/180 90*pi/180])
    xticklabels({'$-90^{\circ}$','$-80^{\circ}$','$-70^{\circ}$','$-60^{\circ}$','$-50^{\circ}$','$-40^{\circ}$','$-30^{\circ}$','$-20^{\circ}$','$-10^{\circ}$','0','$10^{\circ}$','$20^{\circ}$','$30^{\circ}$','$40^{\circ}$','$50^{\circ}$','$60^{\circ}$','$70^{\circ}$','$80^{\circ}$','$90^{\circ}$'});
    set(gca,'TickLabelInterpreter', 'latex')
    hold off
    sgtitle('Response Variation Optimization')
    
    figure()
    xr = theta_plot;
    yr = linspace(interest(1)*2/Fs,interest(2)*2/Fs,sample);
    [xx,yy] = meshgrid(xr,yr);
    brv = 20*log10(abs(B_RV));
    mesh(yy,xx,brv)
    xlim([interest(1)*2/Fs interest(2)*2/Fs])
    ylim([-pi/2 pi/2])
    xlabel('normalized frequency (normalize to sampling frequency)')
    yticks([-90*pi/180 -80*pi/180 -70*pi/180 -60*pi/180 -50*pi/180 -40*pi/180 -30*pi/180 -20*pi/180 -10*pi/180 0 10*pi/180 20*pi/180 30*pi/180 40*pi/180 50*pi/180 60*pi/180 70*pi/180 80*pi/180 90*pi/180])
    yticklabels({'$-90^{\circ}$','$-80^{\circ}$','$-70^{\circ}$','$-60^{\circ}$','$-50^{\circ}$','$-40^{\circ}$','$-30^{\circ}$','$-20^{\circ}$','$-10^{\circ}$','0','$10^{\circ}$','$20^{\circ}$','$30^{\circ}$','$40^{\circ}$','$50^{\circ}$','$60^{\circ}$','$70^{\circ}$','$80^{\circ}$','$90^{\circ}$'});
    set(gca,'TickLabelInterpreter', 'latex')
end

