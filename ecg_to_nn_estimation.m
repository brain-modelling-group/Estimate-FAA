function [out, nn_out] = ecg_to_nn_estimation(ecg,fs,plotopt)

%INPUTS
%ecg = the ECG epoch in samples (as per fs);
%fs = sampling rate of ecg
%plotopt = 'selectepoch' offers an option to select epochs manually (see
%USAGE below) or 'all_epochs' will calculate successive NN intervals without
%inspection (will depend on ECG quality; not built for noisy ECG and may
%yield spurious results)
%

%OUTPUTS
%out = 1 refers to an accepted epoch
%out = 0 refers to a rejected epoch
%nn_out = calculated NN


%USAGE
% This semiautomated script requires button presses to screen valid or invalid ECG epochs
% Each number inputted below on the figure window completes the following
% 11 - add single R peak
% 22 - remove single R peak
% 33 - add multiple R peaks (5 in total)
% 44 - overall assessment this is a valid 5 minute epoch
% 55 - overall assessment this is an invalid 5 minute epoch
% if accepting the epoch or rejecting the epoch without modification press
% either 44 or 55
%
% The screening method will follow through all 5 minute epochs available
% Following the final 5 minute epoch for a given recording, the function will return to command
% window whereby 'out' and 'nn_out' can be saved and used for analysis

% notes:
% sc1 = outlier threshold - if change in nn-interval is greater than
% outlier threshold then ignore nearby r peaks, typically 15
% sc2 = positive negatice difference, threhsold will be different for
% increases in nn-interval than decreases in nn-interval - typically 0.5
% sc3 = % of dodgy R peaks detected to be considered as a bad segent of HRm
% typically 5% or 0.05.

%highpass filter 5 hz
[b,a] = butter(6,(5/(fs/2)), 'high');
dum = conv(double(ecg),[1 1 1])/3;
test_ecg1 = filter(b,a, dum(2:end-1));
exp1 = 100; exp2 = 50;

epoch_size = 5*60;
overlap = epoch_size/2;
block_no = floor(length(test_ecg1)/(overlap*fs))-1;
aa = zeros(block_no, 2);  test_ecg = zeros(size(test_ecg1));
for ii = 1:block_no
    r1 = (ii-1)*overlap*fs+1; r2 = r1+epoch_size*fs-1;
    aa(ii,:) = quantile(test_ecg1(r1:r2), [0.01 0.99]);
    if abs(aa(ii,1))>abs(aa(ii,2)); test_ecg(r1:r2) = -test_ecg1(r1:r2); else; test_ecg(r1:r2) = test_ecg1(r1:r2); end
end

dum1 = abs(hilbert(test_ecg));
dum2 = conv(dum1, ones(1,exp1+1))/(exp1+1); % related to exp1 threshold given below
dum2 = dum2(exp1/2+1:end-exp1/2);
test_ecg = test_ecg./dum2;

N = length(test_ecg);
block_no = floor(length(test_ecg)/(fs*0.5*60)); trr2 = [];
c = -1/(sqrt(2)*erfcinv(3/2)); trr2 = [];
% short 30s segments here
for ii = 1:block_no+1
    if ii>block_no %reach the end of the recording
        r1 = r2; r2 = length(test_ecg); %truncate the final block
        bn1 = floor(length(r1:r2)/4/fs);
    else
        r1 = (ii-1)*fs*0.5*60; r2 = r1+fs*0.5*60;%make a 5 min block isolating the recording at time 5*ii mins
        bn1 = 5*60/4;
    end
    %trx = trr1(find(trr1>r1 & trr1<=r2));%find the peaks that fall within block ii
    tdum = test_ecg(r1+1:r2).^2;  % just take the filtered ECG signal
    %if flt == 5 | flt == 6; tdum = sqrt(tdum); end
    %exp1 = 100; % set threshold ~ 150 bpm - can be generic or maybe based on prelim Pan Tomkins analysis
    if isempty(tdum)==0
        % find peaks using simple thresholding
        th = linspace(quantile(tdum, 0.1), max(tdum), 100);
        erf1 = zeros(1, length(th)-1);
        % run all thresholds - cost function to optimise is the number of
        % detected peaks
        for kk = 1:length(th)-1
            dum = zeros(1, length(tdum));
            dum(tdum>th(kk)) = 1;
            rx = find(diff([0 dum 0]) == 1);
            [erf1(kk),~] = rr_count(rx, tdum, floor(exp1/4)); % this counts peaks above the threshold
        end
        th = th(erf1>(floor(length(tdum)/exp1))/2 & erf1<(floor(length(tdum)/exp1))*2); % this only looks are realistic peak numbers related to previous exp1 threshold
        erf1 = erf1(erf1>(floor(length(tdum)/exp1))/2 & erf1<(floor(length(tdum)/exp1))*2);
        % find the threshold that results in the smallest change in peak number
        % over the longest range of thresholds (the most stable peak number)
        if isempty(th)==0
            q1 = find([erf1 0]<[0 erf1])-1;
            if isempty(q1)==1
                th1 = th(ceil(length(th)/2));
            else
                [q2, q3]= sort(diff(q1), 'descend');
                q4 = find(erf1(q1(q3))>1, 1);
                if isempty(q4)==0
                    th1 = th(q1(q3(q4)));
                else
                    th1 = th(ceil(length(th)/2));
                end
            end
            dum = zeros(1, length(tdum)); dum(tdum>th1) = 1;rx = find(diff([0 dum 0]) == 1);
            [~, trx] = rr_count(rx, tdum, floor(exp1/2));
            trx = trx(tdum(trx)>=median(tdum(trx))/4); % this makes sure detected peaks are of similar amplitude
            trr2 = [trr2 trx+r1];
        end
    end
end

[~, trr2] = rr_count(trr2, test_ecg, floor(exp2/2));

trr1 = trr2; trr2 = [];
block_no = floor(length(test_ecg)/(fs*5*60));
sc1 = 0.75; sc2 = 0.5;
for ii = 1:block_no+1
    if ii>block_no %reach the end of the recording
        r1 = r2; r2 = max(trr1+1); %truncate the final block
        bn1 = floor(length(r1:r2)/4/fs);
    else
        r1 = (ii-1)*fs*5*60; r2 = r1+fs*5*60;%make a 5 min block isolating the recording at time 5*ii mins
        bn1 = 5*60/4;
    end
    trx = trr1(find(trr1>r1 & trr1<=r2)); %find the peaks that fall within block ii
    A = diff(diff(trx));
    th1 = sc1*median(diff(trx)); %sc1*c*median(abs(A-median(A,'omitnan')),'omitnan'); %threshold for change in Heart Rate
    ref = find((A-median(A,'omitnan'))>th1 | (A-median(A,'omitnan'))<-sc2*th1);
    if isempty(ref)==0; trx(ref+1) = NaN; trx(ref+2) = NaN; end
    trr2 = [trr2 trx];
end


if trr2(1)==1; trr2 = trr2(2:end); end
if trr2(end)==N; trr2 = trr2(1:end-1); end
trr = trr2;
nn = trr2;


%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%


nn = nn(isnan(nn)==0);
epoch_size = 5*60;
overlap = epoch_size/2;
block_no = floor(length(test_ecg)/(overlap*fs))-1;
out = zeros(1,block_no);
t = 0:1/fs:(length(test_ecg)-1)/fs;
for ii = 1:block_no
    r1 = (ii-1)*overlap*fs+1; r2 = r1+epoch_size*fs-1;
    %
    ref = find(nn > r1 & nn <= r2);
    if isempty(ref)==0
        acti = 0;
        switch plotopt
            case 'selectepoch'
                ii
                
                try
                    while acti<4
                        figure; subplot(2,1,1); hold on;
                        title('Press "44" to accept epoch or "55" to reject epoch')
                        
                        nn_int = diff(nn(ref));
                        plot(nn(ref(1:end-1))/fs, nn_int)
                        plot([nn(ref(1)) nn(ref(end-1))]/fs, median(nn_int).*[1 1],'k')
                        plot([nn(ref(1)) nn(ref(end-1))]/fs, median(nn_int).*[2 2], 'k--')
                        plot([nn(ref(1)) nn(ref(end-1))]/fs, median(nn_int).*[0.5 0.5], 'k--')
                        axis([nn(ref(1))/fs nn(ref(end-1))/fs median(nn_int)*0.4 median(nn_int)*2.1])
                        subplot(2,1,2); hold on;
                        %rz = find(isnan(nn_int)==0);
                        plot(t(r1:r2), test_ecg(r1:r2))
                        plot(t(nn(ref)), test_ecg(nn(ref)), 'o')
                        axis([r1/fs r2/fs 1.1*min(test_ecg(r1:r2)) 1.1*max(test_ecg(r1:r2))])
                        
                        set(gcf, 'Position', [300 300 1080 420])
                        drawnow;
                        pause
                        k = waitforbuttonpress;
                        acti = str2num(get(gcf,'CurrentCharacter'));
                        if acti == 1
                            aa = ginput(1); bb = round(aa(1)*fs);
                            cc = find(test_ecg(bb-10:bb+10)==max(test_ecg(bb-10:bb+10)));
                            nn = unique(sort([nn bb-11+cc]));
                        end
                        if acti == 2
                            aa = ginput(1); bb = round(aa(1)*fs);
                            rxx = find(abs(nn-bb)==min(abs(nn-bb)));
                            nn(rxx) = 0; nn = nn(nn~=0);
                        end
                        if acti==4; out(ii) = 1; end
                        if acti==5; out(ii) = 0; end
                        if acti==3
                            aa = ginput(5); bb = round(aa(:,1)*fs);
                            for pp = 1:5
                                cc = find(test_ecg(bb(pp)-10:bb(pp)+10)==max(test_ecg(bb(pp)-10:bb(pp)+10)));
                                nn = unique(sort([nn bb(pp)-11+cc]));
                            end
                        end
                        close
                    end
                catch
                end
            case 'all_epochs'
                out(ii) = 1;
        end
    end
end


nn_out = zeros(size(nn));
for jj = 1:length(nn)
    tval = [nn(jj)-1 nn(jj) nn(jj)+1];
    tval1 = [1 2 3]; a = zeros(3);
    for ii = 1:3
        a(ii,1) = tval1(ii).^2;
        a(ii,2)= tval1(ii);
        a(ii,3)=1;
    end
    y = test_ecg(tval);
    b = inv(a)*y';
    nn_out(jj) = -b(2)/(2*b(1))+nn(jj)-2;
end

end



