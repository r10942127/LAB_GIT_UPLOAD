function [y,w] = Frost(Cs,Ci,gs,gi,x,M,J,Fs,freq_s,freq_i,d)
    % hyperparameter 
    mu = 5*10^(-6);
    
    max_freq = max(max(freq_s),max(freq_i));
    AutoFindInterfere = input('interfere signal的DOA是否已知? (Y/N) ','s');
    if AutoFindInterfere=='N'
        C = [Cs Ci];
        g = [gs;gi];
        freq = [freq_s freq_i];
    else
        C = Cs;
        g = gs;
        freq = freq_i;
    end
    IT = floor(length(x(1,:))/Fs);
    P = eye(J*M) - C*((C'*C)\(C'));
    w = C*((C'*C)\g);
    L = C*((C'*C)\g);
    Rx = zeros(J*M,J*M);
    y = zeros(1,IT*Fs);
    for i=1:IT*Fs-J+1
        xx=[];
        for j=1:J
            if i-(j-1)>0
                xx = [xx;x(:,i-(j-1))];
            else
                xx = [xx;zeros(M,1)];
            end
        end
        y(i) = xx.'*conj(w);
        Rx = (Rx*(i-1) + xx*(xx'))/i;
        w_temp = L + P*(w-mu*conj(y(i))*xx);
        w = w_temp;
    end
    plot_num=180;
    theta_plot = linspace(-pi/2,pi/2,plot_num);
    a_plot = zeros(J*M,plot_num,IT);
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
    B_LCMV = zeros(IT,plot_num);
    for k=1:IT
        for i=1:plot_num
            B_LCMV(k,i) = (w')*a_plot(:,i,k);
        end
    end
    
    % plot
    figure()
    subplot(2,1,1)
    plot(abs(y))
    title('Iteration of |y(t)|')
    
    subplot(2,1,2)
    for k=1:IT
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
    sgtitle('Frost') 
end

