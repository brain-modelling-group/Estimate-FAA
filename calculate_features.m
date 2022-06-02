function [features_hrv,flist] = calculate_features(data, fs)
% INPUTS
% data = NN intervals from ECG epoch screening
% fs = sampling rate

% OUTPUTS
% features_hrv = computed HRV features (based on metrics below) on epoched NN
% intervals
% flist = feature list 
flist{1} = 'Median NN';
flist{2} = 'SD NN';
flist{3} = 'VLF';
flist{4} = 'LF';
flist{5} = 'HF';
flist{6} = 'SampEn';
flist{7} = 'SD1';
flist{8} = 'SD2';
flist{9} = 'Higuchi';
flist{10} = 'MSE_En1_SE';
flist{11} = 'MSE_En1_SumSE';
flist{12} = 'MSE_En1_MaxSE';
flist{13} = 'MSE_En1_DiffSE';
flist{14} = 'Median Spectral entropy';
flist{15} = 'Max Spectral entropy';
flist{16} = 'PLxmin';
flist{17} = 'PLalpha';
flist{18} = 'PLBSKc';
flist{19} = 'PLBSKm';
flist{20} = 'BAavg';
flist{21} = 'BWHc';
flist{22} = 'BWHm';
flist{23} = 'BN';
flist{24} = 'BN*BAavg';
flist{25} = 'EntropyShan';
flist{26} = 'EntropyRenyi';
flist{27} = 'InstCWT';
flist{28} = 'PermEn';
flist{29} = 'ApproxEn';
flist{30} = 'LZC';
flist{31} = 'skewness';
flist{32} = 'kurtosis';
flist{33} = 'fuzz_en';
flist{34} = 'cond_en';
flist{35} = 'second_moment';
flist{36} = 'third_moment';
flist{37} = 'fourth_moment';
flist{38} = 'max_minus_min_nni';
flist{39} = 'areaunder_nni';
flist{40} = 'NLEO_median';
flist{41} = 'nni_5th_centile';
flist{42} = 'nni_50th_centile';
flist{43} = 'nni_95th_centile';
flist{44} = 'hurst_exponent';
flist{45} = 'VLF normalized to SDNN';
flist{46} = 'LF normalized to SDNN';
flist{47} = 'HF normalized to SDNN';
flist{48} = 'LF / HF Ratio';
flist{49} = 'LF / (LF+VLF) Ratio';
flist{50} = 'LF / (LF+HF)Ratio';

filt_size=3;
nn = data(1:end-1)/fs;
nni = conv(diff(data)/fs, ones(1,filt_size));
nni = nni(3:end-2);
pNa = sum(isnan(nni))/length(nni);
nn = nn(isnan(nni)==0);
nni = nni(isnan(nni)==0);

mnn = mean(nni);
sdnn = std(nni);  % Standard Devition of NN interval in seconds
[pxx,f] = plomb(nni-mean(nni), nn, 1/(2*min(diff(nn)))); % lomb periodogram will estimate the PSD you can extract power from
pxx = 2*pxx*var(nni)/mean(pxx);
band1 = find(f > 0.01 & f <= 0.04); % define band
band2 = find(f > 0.04 & f <= 0.2); % define band
band3 = find(f > 0.2 & f <= 1); % define band
df1 = (0.04-0.01)/max(f);
df2 = (0.2-0.04)/max(f);
df3 = (1-0.2)/max(f);
bp1 = mean(pxx(band1))*df1/(mean(pxx))*sqrt(mean(pxx)/2); % VLF band power
bp2 = mean(pxx(band2))*df2/(mean(pxx))*sqrt(mean(pxx)/2); % LF band power
bp3 = mean(pxx(band3))*df3/(mean(pxx))*sqrt(mean(pxx)/2); % HF power
se_samp = SampEn(2, 0.2*sdnn, nni); % Sample Entropy NN interval
[sd1, sd2] = poincare(nni); % Poincare Analysis NN interval
[fd,~,~] = Higuchi1Dn(nni);
msce = multiscale_entropy(nni);
[pse,te] = pentropy(nni,1:length(nni));
med_se = median(pse);
max_pse = max(pse);
k=PLmeasure(nni,1:length(nni),fs);
[shannonEntr, renyiEntr] = shanRenEntropy(nni,4,3);
[wt,~,~] = cwt_an(nni,fs); inst_spec=trapz(median(abs(wt),1));
max_permEn = PermEn(nni,50);
Ap_En=ApEn(nni,2,1);
fuzz_en = FuzzyEn(nni,2,0.2*sdnn,50);
cond_en = CondEn(nni,2,6);
[lzc, ~] = calc_lz_complexity(nni,'primitive', 1);
nni_sk = skewness(nni);
nni_ku = kurtosis(nni);

m2 = moment(nni,2);
m3 = moment(nni,3);
m4 = moment(nni,4);

maxmin_nni=max(nni)-min(nni);
area_nni=trapz(nni);

median_snleo = median(nlin_energy(nni));
nni_centile=quantile(nni,[0.05 0.5 0.95]);

hurst_nn = estimate_hurst_exponent(nn);

vlf_norm_to_sdnn = bp1*1000./(sdnn*1000)*100;
lf_norm_to_sdnn = bp2*1000./(sdnn*1000)*100;
hf_norm_to_sdnn = bp3*1000./(sdnn*1000)*100;
lf_hf_ratio = bp2*1000./bp3*1000;
lf_lfplusvlf_ratio = bp2*1000./(bp1*1000+bp2*1000);
lf_lfplushf_ratio = bp2*1000./(bp2*1000+bp3*1000);

% feature output
features_hrv = [mnn*1000 sdnn*1000, bp1*1000, bp2*1000, bp3*1000, se_samp, sd1, sd2, fd, msce(1), msce(2), msce(3), msce(4), med_se, max_pse, k(1), k(2), k(3), k(4), k(5), k(6), k(7), k(8), k(9),shannonEntr,renyiEntr,inst_spec,max_permEn,Ap_En,lzc,nni_sk,nni_ku,fuzz_en,cond_en,m2,m3,m4,maxmin_nni,area_nni,median_snleo, nni_centile(1),nni_centile(2),nni_centile(3),hurst_nn,vlf_norm_to_sdnn,lf_norm_to_sdnn,hf_norm_to_sdnn,lf_hf_ratio,lf_lfplusvlf_ratio,lf_lfplushf_ratio]; 


