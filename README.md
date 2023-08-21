# SSVEP-TransRCA
A transfer learning based SSVEP classification algorithm


### Usage：

% X_tar: training data of target domain (target * channel * sample * block )\
% X_src: training data of source domain( target * channel * sample * block * source)\
% fs： sampling rate\
% freq: frequency list of stimuli\
% phase: phase list of stimuli\
% num_fbs: number of filter bank

% test_data: test data of target domain (target * channel * sample )\

model = TransRCA_train(X_tar,X_src,fs,freq_list,phase_list,num_fbs);

[results, rho_final, results_6r]  = TransRCA_test_fast(test_data,model,is_ensemble)
