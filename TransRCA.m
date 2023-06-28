function [w_a,w_b,r] = TransRCA(X_t,X_s)
%  Transfer Related Compont Analysis (TransRCA)
% Input:
%   X_t     : EEG data of target domain 
%            (# of  channel_t , sample , trial_t)
%   X_s     : EEG data of source domain 
%            (# of  channel_s , sample , trial_s)

X_t = squeeze(X_t);
X_s = squeeze(X_s);

% sum of variance of target domain
[num_chan_t, num_samp, num_tria_t]  = size(X_t);

UX = reshape(X_t, num_chan_t, num_samp*num_tria_t);
UX =  UX-mean(UX,2);
S11 = UX*UX'/num_tria_t;

% sum of variance of source domain
[num_chan_s, num_samp, num_tria_s]  = size(X_s);
UX = reshape(X_s, num_chan_s, num_samp*num_tria_s);
UX = UX-mean(UX,2);
S22 = UX*UX'/num_tria_s;


% sum of covariance of target domain and source domain(faster)
X_t = X_t - mean(X_t,2);
X_s = X_s - mean(X_s,2);

S12 = pagemtimes(X_t,permute(X_s,[2,1,4,3]));
S12 = sum(S12,[3,4])/(num_tria_t*num_tria_s);
S21=S12';
[w_a,r] = eigs(S11^-1*S12*S22^-1*S21);
[w_b,~] = eigs(S22^-1*S21*S11^-1*S12);

r=diag(r);



% % sum of covariance of target domain and source domain(slower)
% S12 = zeros(num_chan_t,num_chan_s);
% S21 = zeros(num_chan_s,num_chan_t);
% for trial_i = 1:1:num_tria_t
%     x1 = squeeze(X_t(:,:,trial_i));
%     x1 = bsxfun(@minus, x1, mean(x1,2));
%     
%  
%     for trial_j = 1:1:num_tria_s
%         x2 = squeeze(X_s(:,:,trial_j));
%         x2 = bsxfun(@minus, x2, mean(x2,2));
%         S12 = S12 + x1*x2';
%     end % trial_j
%     
% end % trial_i
% S12 = S12/(num_tria_t*num_tria_s);
% S21=S12';
% % S21 = S21/(num_tria_t*num_tria_s);
% S21=S12';
% [w_a,~] = eigs(S11^-1*S12*S22^-1*S21);
% [w_b,~] = eigs(S22^-1*S21*S11^-1*S12);


end