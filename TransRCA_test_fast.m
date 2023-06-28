function [results,rho_final,results_6r]  = TransRCA_test_fast(test_data,model,is_ensemble)
% test_data: targ * chan * samp

[num_targ,num_fbs,num_chan_t,num_samp] = size(model.mean_X_t);
[num_targ,num_fbs,num_chan_s,num_samp,num_src] = size(model.mean_X_s);

% filterbank
fb_coefs = [1:model.num_fbs].^(-1.25)+0.25;

x_tar = permute(model.mean_X_t,[3,4,2,1]); % chan* samp * fb* class

x_src = permute(model.mean_X_s,[3,4,2,1,5]); % chan * samp  * fb * class  * src

for targ_i = 1:num_targ
    test_tmp = squeeze(test_data(targ_i, :, :));

    x_test = zeros(num_chan_t,num_samp,model.num_fbs);

    for fb_i = 1:model.num_fbs
        x_test(:,:,fb_i) = filterbank(test_tmp,model.fs,fb_i);
        for class_i = 1:num_targ
            %% r1
            Y=squeeze(model.Template_Sin(class_i,:,:));
            [~,~,r1_tmp] = canoncorr( x_test(:,:,fb_i)',Y'); 
            r1(fb_i,class_i)=r1_tmp(1);
        end
    end

    %num_targ,num_chan_s,num_fbs,num_bloc,num_src
    w_tar2ref = model.w_tar2ref;
    w_src2ref = model.w_src2ref;
    w_tar2src = model.w_tar2src;
    w_src2tar = model.w_src2tar;

%     w_ref2tar = model.w_ref2tar;


    
    if is_ensemble
        % spitial filtering of test data
        cv1 = pagemtimes(w_tar2ref,x_test); % 40filters 125sample 5fb，1test bloc  1target subject
        cv2 = pagemtimes(w_tar2src,x_test); % 40filters 125sample 5fb，1test bloc  10source subject
      
        % spitial filtering of target domain and source domain template
        cv3 = pagemtimes(w_tar2ref,x_tar); % 40filters 125sample 5fb 40template  1target subject
        cv4 = pagemtimes(w_src2ref,x_src); % 40filters 125sample 5fb 40template  10source subject
        cv5 = pagemtimes(w_tar2src,x_tar); % 40filters 125sample 5fb 40template  10source subject
        cv6 = pagemtimes(w_src2tar,x_src); % 40filters 125sample 5fb 40template  10source subject

    else
        for class_i=1:num_targ
            %spitial filtering of test data
            cv1(class_i,:,:) = pagemtimes(w_tar2ref(class_i,:,:),x_test); % 1filters 125sample 5fb，1test bloc  1target subject
            cv2(class_i,:,:) = pagemtimes(w_tar2src(class_i,:,:),x_test); % 1filters 125sample 5fb，1test bloc  10source subject
            
            % spitial filtering of target domain and source domain template
            cv3(:,:,:,class_i) = pagemtimes(w_tar2ref(class_i,:,:),x_tar(:,:,:,class_i)); % 1filters 125sample 5fb 1template  1target subject
            cv4(:,:,:,class_i) = pagemtimes(w_src2ref(class_i,:,:),x_src(:,:,:,class_i)); % 1filters 125sample 5fb 1template  10source subject
            cv5(:,:,:,class_i) = pagemtimes(w_tar2src(class_i,:,:),x_tar(:,:,:,class_i)); % 1filters 125sample 5fb 1template  10source subject
            cv6(:,:,:,class_i) = pagemtimes(w_src2tar(class_i,:,:),x_src(:,:,:,class_i)); % 1filters 125sample 5fb 1template  10source subject
        end
    end

    r2 = corr2_fast(cv1,cv3);
    r3 = corr2_fast(cv1,cv4);
    r4 = corr2_fast(cv2,cv5);
    r5 = corr2_fast(cv2,cv6);
    
    r = cat(3,r1,r2,mean(abs(r3),3),mean(r4,3),mean(abs(r5),3)); % fusion of all correlation coefficient
   
    r_fb = squeeze(sum(r.*fb_coefs',1)) ;  % filter bank 
    rho = sum(r_fb,2);

    [~,tau] = max(rho);
    results(targ_i,:)=tau;
    rho_final(targ_i,:) = rho;

    rho_all = cat(2,r_fb,rho);
    [~,tau] = max(rho_all);
    results_6r(targ_i,:)=tau;

end % targ_i

end

function r =  corr2_fast(a,b)
a = a - mean(a,[1,2]);  
b = b - mean(b,[1,2]);

cab = sum(a.*b,[1,2]);  
caa = sum(a.*a,[1,2]);
cbb = sum(b.*b,[1,2]);
r = squeeze(cab./sqrt(caa.*cbb));

end



