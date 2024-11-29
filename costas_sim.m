% Made by Trans-cending (yj Bian from Southeast University)

clear all
close all
clc


%% Parameters
rng(37)
Fs=16.384e6;  % sampling freq
fc=1.024e6;   % carrier freq
Rb=1.024e6;   % bps
phase_off=-pi/6;  % phase offset
nrep = Fs/Rb;
M=4;    % qpsk, M=4
nbit=log2(M);
RB=Rb/log2(M);  % Baud
nsym=8e5;  % total symbols emulated
order = 200;
t_max=nsym/Rb;
ff=[0 fc/(Fs/2) 1.5*fc/(Fs/2) 1];   % LPF design
fa=[1 1 0 0];
c1=5;   % loop filter coefficient
c2=3;
c1_new=2;
c2_new=1;

%% Code sequence generator
Codebook=randi([0,1],[nbit,nsym]);

nsym=nsym+1;
Codebook_diff=zeros([nbit,nsym]);
for i=1:nsym-1
    if Codebook_diff(1,i)==Codebook_diff(2,i)
        Codebook_diff(1,i+1)=xor(Codebook(1,i),Codebook_diff(1,i));
        Codebook_diff(2,i+1)=xor(Codebook(2,i),Codebook_diff(2,i));
    else
        Codebook_diff(1,i+1)=xor(Codebook(2,i),Codebook_diff(1,i));
        Codebook_diff(2,i+1)=xor(Codebook(1,i),Codebook_diff(2,i));
    end
end
Codebook_rep=repmat(Codebook_diff,[nrep,1]);
Codebook_trans=reshape(Codebook_rep,[nbit,nsym*nrep]);



%% Modulation of B type
t_vec=linspace(0,t_max,nsym*nrep);
Codebook_trans(Codebook_trans==0)=-1;
I_flow=Codebook_trans(1,:) * sqrt(2)/2;
Q_flow=Codebook_trans(2,:) * sqrt(2)/2;

trans_sig=I_flow.*cos(2*pi*fc.*t_vec+phase_off)+Q_flow.*sin(2*pi*fc.* ...
    t_vec+phase_off);
% figure
% plot(trans_sig(1:50))
% title('transmit signal qpsk')

%% Add noise
receiv_sig=awgn(trans_sig,5);
receiv_cos=cos(2*pi*fc.*t_vec+phase_off);

h=firpm(order,ff,fa); % LPF
h_fft=fftshift(fft(h));
% figure
% plot(abs(h_fft));

%% Simulation
point_thres=floor(1.5e-3/t_max*nsym*nrep);  % change coefficient after threshold 
phase=zeros([nsym*nrep+1, 1]);
phase_diff=zeros([nsym*nrep, 1]);
lpfi=zeros([nsym*nrep, 1]);
lpfq=zeros([nsym*nrep, 1]);
costa_cos=zeros([nsym*nrep, 1]);
costa_sin=zeros([nsym*nrep, 1]);
e=zeros([nsym*nrep, 1]);
begin_freq=2*pi*fc;
deviation_freq=2*pi*100;        % freq deviation
phase(1)=begin_freq+deviation_freq;
zi=zeros(1,order+1); 
zq=zeros(1,order+1); 
counter=1;  
counter_out = 8;
loop_filter = zeros([floor(nsym*nrep/counter_out)+2,1]);
step = 11*Fs/fc;
orderc=order/2;
flag=0;
k=2;
Codebook_diff_recovery=zeros([nbit, nsym*nrep+1]);
figure('Position', [20, 50, 1500, 600]); 
sgtitle(sprintf('RResults: c1=%.1f, c2=%.1f, c1_{new}=%.1f, c2_{new}=%.1f', c1, c2, c1_new, c2_new))
for i=1:nsym*nrep
    time=t_vec(i);
    phase_diff(i)=-phase(i)*time+(2*pi*fc.*time+phase_off);
    costa_cos(i)=cos(phase(i)*time);
    costa_sin(i)=sin(phase(i)*time);
    mul_i=costa_cos(i)*trans_sig(i);
    mul_q=costa_sin(i)*trans_sig(i);
    zi=[zi(2:order+1), 2*mul_i];
    zq=[zq(2:order+1), 2*mul_q];
    lpfi(i)=fliplr(h)*zi';
    lpfq(i)=fliplr(h)*zq'; 
%     if lpfi(i)>0
%         Codebook_diff_recovery(1,i+1)=1;  
%     else
%         Codebook_diff_recovery(1,i+1)=-1;
%     end
%     if lpfq(i)>0
%         Codebook_diff_recovery(2,i+1)=1;  
%     else
%         Codebook_diff_recovery(2,i+1)=-1;
%     end
    % error function definition
    Codebook_diff_recovery(1,i+1)=sign(lpfi(i));
    Codebook_diff_recovery(2,i+1)=sign(lpfq(i));
    e(i)=sign(lpfq(i))*lpfi(i)-sign(lpfi(i))*lpfq(i);
%     e(i)=-sign(lpfq(i))*lpfi(i)+sign(lpfi(i))*lpfq(i);
    % loop filter
    if i<=point_thres
        disp('larger coefficients')
        if counter==1
            loop_filter_new = c2*e(i);
        elseif counter==2
            loop_filter_new = loop_filter_new + c1*e(i);
        end
    else
        disp('smaller coefficients')
        if counter==1
            loop_filter_new = c2_new*e(i);
        elseif counter==2
            loop_filter_new = loop_filter_new + c1_new*e(i);
        end
    end

    if counter == counter_out
        loop_filter(k)=loop_filter_new;
        counter=1;
        k=k+1;
    end

    phase(i+1)=loop_filter(k-1)+phase(1);
%     phase(i+1)=phase(1);
    counter = counter+1;



    % visualization
    if i>orderc+step+1
        if mod(i,step)==0
            subplot(2,3,1)
            plot(t_vec(1:i), -phase(1:i)/(2*pi)+fc,LineWidth=1.5)
            grid on
            title('freq deviation')
            grid on
            xlabel('t/s')
            ylabel('freq/Hz')
            subplot(2,3,2)
            plot(t_vec(1:i),phase_diff(1:i)*180/pi,LineWidth=1.5)
            grid on
            title('phase deviation')
            xlabel('t/s')
            ylabel('phase/^\circ')
            subplot(2,3,3)
            plot(t_vec(i-step+1-orderc:i-orderc),receiv_sig(i-step+1:i),LineWidth=1.5)
            title('QPSK receive signal')
            xlabel('t/s')
            ylabel('amplitude')
            subplot(2,3,4)
            plot(t_vec(i-step+1-orderc:i-orderc),costa_cos(i-step+1:i),t_vec(i-step+1-orderc:i-orderc), receiv_cos(i-step+1:i),LineWidth=1.5)
            legend('recover', 'expected')
            title('carrier')
            xlabel('t/s')
            ylabel('amplitude')
            subplot(2,3,5)
            plot(t_vec(i-step+1-orderc:i-orderc),sqrt(2)*I_flow(i-step+1-orderc:i-orderc)+1.5,t_vec(i-step+1-orderc:i-orderc),Codebook_diff_recovery(1,i-step+1+1:i+1)-1.5,LineWidth=1.5)
            legend('I transmitted bit', 'I recovered bit')
            xlabel('t/s')
            ylabel('amplitude')
            title('I branch')
            subplot(2,3,6)
            plot(t_vec(i-step+1-orderc:i-orderc),sqrt(2)*Q_flow(i-step+1-orderc:i-orderc)+1.5,t_vec(i-step+1-orderc:i-orderc),Codebook_diff_recovery(2,i-step+1+1:i+1)-1.5,LineWidth=1.5)
            legend('Q transmitted bit', 'Q transmitted bit')
            xlabel('t/s')
            ylabel('amplitude')
            title('Q branch')
            pause(0.01)
        end
        if (((abs(sum(phase_diff(i-step+1:i))*180/pi/step)<5) && ...
            (abs(sum(phase(i-step+1:i)-2*pi*fc)/step)<10*2*pi)) || ...
            (((abs(sum(phase_diff(i-step+1:i))*180/pi/step)<92.5) && ...
            ((abs(sum(phase_diff(i-step+1:i))*180/pi/step)>87.5)))  && ...
            (abs(sum(phase(i-step+1:i)-2*pi*fc)/step)<10*2*pi))) && ~flag
            fprintf('final phase deviation: %f \n', phase_diff(i))
            fprintf('final freq deviation: %f Hz\n', (phase(i)-2*pi*fc)/(2*pi))
            flag=1;
            j=1;
        end
        if flag
            j=j+1;
            if j>step*20    
                break
            end
        end
    end
end








