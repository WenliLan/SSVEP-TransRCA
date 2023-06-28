function model = TransRCA_train(X_t,X_s,fs,freq,phase,num_fbs)
% X_tar: training data of of target domain (target * channel * sample * block )
% X_src: training data of source domain( target * channel * sample * block * source)
% fs： sampling rate
% freq: frequency list of stimuli
% phase: phase list of stimuli
% num_fbs: number of filter bank 

% Constract sine-cosine reference
num_harmonic= 5;
[num_targ,num_chan_t,num_samp,num_bloc_t] = size(X_t);
[num_targ,num_chan_s,num_samp,num_bloc_s,num_src] = size(X_s);

t = (0:num_samp-1)/fs;
for targ_i = 1:length(freq)
    Template_Sin(targ_i,:,:) = ref_constract(t, freq(targ_i),phase(targ_i), num_harmonic);
end
%% Concatenation source data to channel domain
% X_src = permute(X_src,[1,2,5,3,4]);
% X_src = reshape(X_src,[num_targ,num_chan*num_sub_src,num_samp,num_bloc]);

%% Concatenation source data to block domain
X_s = reshape(X_s,[num_targ,num_chan_s,num_samp,num_bloc_s*num_src]) ;
[num_targ,num_chan_s,num_samp,num_bloc_s,num_src] = size(X_s);

mean_X_t = zeros(num_targ,num_fbs,num_chan_t,num_samp);         % template of target domain
mean_X_s = zeros(num_targ,num_fbs,num_chan_s,num_samp,num_src); % template of source domain 

%num_targ,num_chan_s,num_fbs,1,num_src
w_tar2ref =  zeros(num_targ,num_chan_t,num_fbs,1,1);  
w_ref2tar =  zeros(num_targ,num_harmonic*2,num_fbs,1,1);  

w_src2ref = zeros(num_targ,num_chan_s,num_fbs,1,num_src);
w_ref2src = zeros(num_targ,num_harmonic*2,num_fbs,1,num_src);

w_tar2src = zeros(num_targ,num_chan_t,num_fbs,1,num_src);
w_src2tar = zeros(num_targ,num_chan_s,num_fbs,1,num_src);

for targ_i = 1:num_targ
    for fb_i = 1:num_fbs
        X_t_tmp = squeeze(X_t(targ_i,:,:,:));
        X_t_tmp = filterbank(X_t_tmp,fs,fb_i);
        mean_X_t(targ_i,fb_i,:,:) = squeeze(mean(X_t_tmp,3));

        % TransRCA between target domian and sin-cosine reference
        [w_tmp1,w_tmp2] = TransRCA(X_t_tmp,squeeze(Template_Sin(targ_i,:,:)));
        w_tar2ref(targ_i,:,fb_i) = w_tmp1(:,1);   
        w_ref2tar(targ_i,:,fb_i) = w_tmp2(:,1);

        for src_i = 1:num_src
            X_s_tmp = squeeze(X_s(targ_i,:,:,:,src_i));
            X_s_tmp = filterbank(X_s_tmp,fs,fb_i);
            
            mean_X_s(targ_i,fb_i,:,:,src_i) = squeeze(mean(X_s_tmp,3));
                    
            % TransRCA between source domian and sin-cosine reference
            [w_tmp1,w_tmp2] = TransRCA(X_s_tmp,squeeze(Template_Sin(targ_i,:,:)));
            w_src2ref(targ_i,:,fb_i,1,src_i) = w_tmp1(:,1);   
            w_ref2src(targ_i,:,fb_i,1,src_i) = w_tmp2(:,1);   

            % TransRCA between target domian and source domian
            [w_tmp1,w_tmp2] = TransRCA(X_t_tmp,X_s_tmp);
            w_tar2src(targ_i,:,fb_i,1,src_i)=w_tmp1(:,1);  %   W_t(targ_i,:,src_sub_i)*mean_X_t
            w_src2tar(targ_i,:,fb_i,1,src_i)=w_tmp2(:,1);

        end %src_i
    end %fb_i
end %targ_i

% save template, spitial filter, filterbank perameter
model = struct('mean_X_t', mean_X_t,'mean_X_s',mean_X_s, 'Template_Sin',Template_Sin,...
    'w_tar2ref',w_tar2ref,'w_ref2tar',w_ref2tar,'w_src2ref',w_src2ref,'w_ref2src',w_ref2src,...
    'w_tar2src', w_tar2src,'w_src2tar',w_src2tar, ...
    'num_fbs',num_fbs,'fs',fs   );
end


