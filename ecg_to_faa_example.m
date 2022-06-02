% Example based on ECG data
addpath([pwd '\subfunctions']); % dependencies for calculation of HRV features 
load ECG_example.mat % this contains sample ECG epochs for a single subject

% Screen ECG epoch(s)
fs=ECG(1).fs; %set sampling rate
for ii=1:length(ECG)
    [nn_valid, nn_out] = ecg_to_nn_estimation(ECG(ii).epoch,fs,'all_epochs');
    nn_valid_ep(ii)=nn_valid;
    nn_out_ep{ii}=nn_out;
end

% Calculate HRV features based on nn_valid_ep = 1
nn_selected = nn_out_ep(nn_valid_ep==1);
features_hrv=zeros(length(nn_selected),50);
for m=1:length(nn_selected)
    [features_hrv(m,:),flist] = calculate_features(nn_selected{m}, fs);
end
m_features_hrv=median(features_hrv,1);

% Estimate Age based on selected HRV features (see "Bedside tracking of
% functional autonomic age in preterm infants")
load training_Mdl_hrv.mat %LOOCV Model trained from infant ECG dataset
faa_predict=predict(Mdl_HRV,m_features_hrv(:,select_feat));
% based on training set
pma=ECG(1).ca;
PMA_vs_FAA=[pma faa_predict]







