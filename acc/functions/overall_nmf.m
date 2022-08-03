function xhat1=overall_nmf(x,fs,MAXITER,K,sparsity,TF,nmf_method,reconstruction,supervised,x1,x2)
%K = 3; % number of basis vectors (2 most likely)
%MAXITER = 100; % total number of iterations to run (50 another option)
%sparsity=1 for instance
%supervised = [1 1];
%K = [25 25]; % number of basis vectors
%MAXITER = 500; % total number of iterations to run
if nargin < 9
    supervised=[0,0];
end
%% TF Representation
switch TF
    case 'Q-transform'
        bins=12; %or 48 maybe 
        [X,~,g,fshifts] = cqt(x,'SamplingFrequency' ,fs,'BinsPerOctave',bins);
        %'BinsPerOctave' ? Number of bins per octave
        %1-96 default 12 
        % 'TransformType' The sparse transform is the minimally redundant
        % version of the constant-Q transform. 'Sparse'
        % 'FrequencyLimits' he first element must be greater than or equal to Fs/N, where Fs is the sampling frequency and N is the length of the signal.
        %V=abs(X); 
        V = abs(X(1:(size(X,1)/2+1),:));
        F = size(V,1);
        T = size(V,2);
        
        if supervised(1)
            X1 = cqt(x1,'SamplingFrequency' ,fs,'BinsPerOctave',bins);
            %V1=abs(X1);
            V1 = abs(X1(1:(size(X1,1)/2+1),:));
        end
        if supervised(2)
            X2 = cqt(x2,'SamplingFrequency' ,fs,'BinsPerOctave',bins);
            %V2=abs(X2); 
            V2 = abs(X2(1:(size(X2,1)/2+1),:));
        end
                
    case 'STFT'
        %% STFT
        FFTSIZE = 1024;
        HOPSIZE = 256;
        WINDOWSIZE = 512;
        
        %X=stft(x,'Window',hann(WINDOWSIZE),'OverlapLength',HOPSIZE,'FFTLength',FFTSIZE);
        X = myspectrogram(x,FFTSIZE,fs,hann(WINDOWSIZE),-HOPSIZE);
        %V=abs(X); 
        V = abs(X(1:(FFTSIZE/2+1),:));
        F = size(V,1);
        T = size(V,2);
        
        % imagesc(db(V))
        % set(gca,'YDir','normal')
        % set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
        % set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
        % title('Spectrogram of Mary Had a Little Lamb')
        % ylabel('Frequency')
        % xlabel('Time')
        %sound(x,fs)
        if supervised(1)
            X1 = myspectrogram(x1,FFTSIZE,fs,hann(WINDOWSIZE),-HOPSIZE);
            V1 = abs(X1(1:(FFTSIZE/2+1),:));
        end
        if supervised(2)
            X2 = myspectrogram(x2,FFTSIZE,fs,hann(WINDOWSIZE),-HOPSIZE);
            V2 = abs(X2(1:(FFTSIZE/2+1),:));
        end
    case 'Cochleagram'
        %https://www.mathworks.com/matlabcentral/fileexchange/48622-cochleagram-and-is-nmf2d-for-blind-source-separation
        fRange = [50, 8000];
        gf = gammatoneFast(x,128,fRange);%Construct the cochleagram use Gammatone filterbank
        cg = cochleagram(gf);
        V=cg;
end


%% NMF
Tau=0:7;           % Defines tau shifts
Phi=0:32;          % Defines phi shifts

opt=0;
switch nmf_method
    case 'NMCF' 
        lambda=1; 
        fixedInds = 1:K(1);
        [W, H] = nmcf(V,V1, K, MAXITER,fixedInds,lambda); 
    case 'Optimize'
        %% NMF
        if supervised(1)||supervised(2)
            if supervised(1)
                [W1, ~]=nmf_optimized(V1,K(1),[],[]); 
            else
                W1 = 1+rand(F, K(1));
            end
            if supervised(2)
                [W2, ~]=nmf_optimized(V2,K(2),[],[]); 
            else
                W2 = 1+rand(F, K(2));
            end
            
            if supervised(1) && supervised(2)
                fixedInds = 1:sum(K);
            elseif supervised(1)
                fixedInds = 1:K(1);
            elseif supervised(2)
                fixedInds = (K(1)+1):sum(K);
            else
                fixedInds = [];
            end
            [W, H]=nmf_optimized(V,K,[W1 W2],fixedInds); 
        else
            fixedInds = [];
            [W, H]=nmf_optimized(V,K,[],fixedInds); 
        end
    case 'Euclidean'
        opt=1;
        beta_loss=2;
        %[W, H]=nmf_general(V,K,beta_loss,MAXITER);
    case 'KL'
        opt=1;
        beta_loss=1;
        %[W, H]=nmf_general(V,K,beta_loss,MAXITER);
    case 'IS'
        opt=1;
        beta_loss=0;
        %[W, H]=nmf_general(V,K,beta_loss,MAXITER);
    case 'Euclidean_Sparse'
        opt=2;
        beta_loss=2;
        %[W,H]=sparse_nmf(V,K,sparsity,beta_loss,MAXITER);
    case 'KL_Sparse'
        opt=2;
        beta_loss=1;
        %[W,H]=sparse_nmf(V,K,sparsity,beta_loss,MAXITER);
    case 'IS_Sparse'
        opt=2;
        beta_loss=0;
        %[W,H]=sparse_nmf(V,K,sparsity,beta_loss,MAXITER);
    case 'KL_2'
        %% NMF
        if supervised(1)||supervised(2)
            if supervised(1)
                [W1, ~] = nmf(V1, K(1), [], MAXITER,[]);
            else
                W1 = 1+rand(F, K(1));
            end
            if supervised(2)
                [W2, ~] = nmf(V2, K(2), [], MAXITER,[]);
            else
                W2 = 1+rand(F, K(2));
            end
            
            if supervised(1) && supervised(2)
                fixedInds = 1:sum(K);
            elseif supervised(1)
                fixedInds = 1:K(1);
            elseif supervised(2)
                fixedInds = (K(1)+1):sum(K);
            else
                fixedInds = [];
            end
            
            [W, H]   = nmf(V, K, [W1 W2], MAXITER, fixedInds);
        else
            fixedInds = [];
            [W, H]   = nmf(V, K, [], MAXITER, fixedInds);
            %[W, H] = nmf(V, K, MAXITER);
        end
    case 'MU'
        [W, H, ~] = is_nmf2D_mu(V,MAXITER,K,Tau,Phi);
    case 'SAGE'
        [W, H, ~] = is_nmf2D_em(V,MAXITER,K,Tau,Phi);
    case 'als'
        [W,H] = nnmf(V,K,'algorithm','als');
    case 'mult'
        [W,H] = nnmf(V,K,'algorithm','mult');
end

if opt==1
    if supervised(1)||supervised(2)
        if supervised(1)
            [W1, ~]=nmf_general(V1,K(1),beta_loss,MAXITER,[],[]);
        else
            W1 = 1+rand(F, K(1));
        end
        if supervised(2)
            [W2, ~]=nmf_general(V2,K(2),beta_loss,MAXITER,[],[]);
        else
            W2 = 1+rand(F, K(2));
        end
        
        if supervised(1) && supervised(2)
            fixedInds = 1:sum(K);
        elseif supervised(1)
            fixedInds = 1:K(1);
        elseif supervised(2)
            fixedInds = (K(1)+1):sum(K);
        else
            fixedInds = [];
        end
        [W, H]=nmf_general(V,K,beta_loss,MAXITER,fixedInds,[W1 W2]);
    else
        fixedInds = [];
        [W, H]=nmf_general(V,K,beta_loss,MAXITER,fixedInds,[]);
    end
elseif opt==2
    if supervised(1)||supervised(2)
        if supervised(1)
            [W1,~]=sparse_nmf(V1,K(1),sparsity,beta_loss,MAXITER,[],[]);
        else
            W1 = 1+rand(F, K(1));
        end
        if supervised(2)
            [W2,~]=sparse_nmf(V2,K(2),sparsity,beta_loss,MAXITER,[],[]);
        else
            W2 = 1+rand(F, K(2));
        end
        
        if supervised(1) && supervised(2)
            fixedInds = 1:sum(K);
        elseif supervised(1)
            fixedInds = 1:K(1);
        elseif supervised(2)
            fixedInds = (K(1)+1):sum(K);
        else
            fixedInds = [];
        end
        [W,H]=sparse_nmf(V,K,sparsity,beta_loss,MAXITER,fixedInds,[W1 W2]);
    else
        fixedInds = [];
        [W,H]=sparse_nmf(V,K,sparsity,beta_loss,MAXITER,fixedInds,[]);
    end
end

%% Reconstruction
switch reconstruction
    case 'Synthesis'
        %% ISTFT / RECONSTRUCTION METHOD 1 (SYNTHESIS)
        % get the mixture phase
        phi = angle(X);
        
        for i=1:K
            
            XmagHat = W(:,i)*H(i,:);
            
            % create upper half of frequency before istft
            XmagHat = [XmagHat; conj( XmagHat(end-1:-1:2,:))];
            
            % Multiply with phase
            XHat = XmagHat.*exp(1i*phi);
            switch TF
                case 'STFT'
                    %xhat1(:,i)=real(istft(XHat,'Window',hann(WINDOWSIZE),'OverlapLength',HOPSIZE,'FFTLength',FFTSIZE));
                    xhat1(:,i) = real(invmyspectrogram(XHat,HOPSIZE))';
                case 'Q-transform'
                    xhat1(:,i) = real(icqt(XHat,g,fshifts))';
            end 
        end
        % sound(xhat1(:,1),fs)
        % sound(xhat1(:,2),fs)
        % sound(xhat1(:,3),fs)
    case 'Filtering'
        %% ISTFT / RECONSTRUCTION METHOD 2 (FILTERING)
        % get the mixture phase
        phi = angle(X);
        if 1==1%supervised(1)||supervised(2)
            % get the mixture phase
            c = [1 cumsum(K)];
            for i=1:length(K)
                % create masking filter
                Mask =  W(:,c(i):c(i+1))*H(c(i):c(i+1),:)./(W*H);
                
                % filter
                XmagHat = V.*Mask;
                
                % create upper half of frequency before istft
                XmagHat = [XmagHat; conj( XmagHat(end-1:-1:2,:))];
                
                % create upper half of frequency before istft
                switch TF
                    case 'STFT'
                         xhat1(:,i) = real(invmyspectrogram(XmagHat.*exp(1i*phi),HOPSIZE))';
                    case 'Q-transform'
                        xhat1(:,i) = real(icqt(XmagHat.*exp(1i*phi),g,fshifts))';
                end
                
                %     sound(xhat(:,i),fs)
            end
        else
            for i=1:K
                
                % create masking filter
                Mask =  (W(:,i)*H(i,:)./(W*H));
                
                % filter
                XmagHat = V.*Mask;
                
                % create upper half of frequency before istft
                XmagHat = [XmagHat; conj( XmagHat(end-1:-1:2,:))];
                
                switch TF
                    case 'STFT'
                        xhat1(:,i) = real(invmyspectrogram(XmagHat.*exp(1i*phi),HOPSIZE))';
                    case 'Q-transform'
                        xhat1(:,i) = real(icqt(XmagHat.*exp(1i*phi),g,fshifts))';
                end
                
                
            end
            % sound(xhat1(:,1),fs)
            % sound(xhat1(:,2),fs)
            % sound(xhat1(:,3),fs)
        end
    case 'Cochleagram'
        for idx = 1:K
            Rec(:,:,idx) = isp_nmf2d_rec(W, H, Tau,Phi, idx);
        end
        for k = 1:size(Rec,3)
            mask = logical(Rec(:,:,k)==max(Rec,[],3));
            r(k,:) = synthesisFast(x,mask,fRange);
            r(k,:) = r(k,:)./max(abs(r(k,:)));
            %     wavwrite(r,Fs,sprintf([wavfile, '%d.wav'], k));
        end
        %         figure
        %         subplot(311),plot(x);
        %         xlabel('Samples');
        %         ylabel('Amplitude');
        %         title('Mixture');
        %         subplot(312),plot(r(1,:));
        %         xlabel('Samples');
        %         ylabel('Amplitude');
        %         title('Estimated source one');
        %         subplot(313),plot(r(2,:));
        %         xlabel('Samples');
        %         ylabel('Amplitude');
        %         title('Estimated source two');
        %         sound(r(1,:),Fs);
        %         sound(r(2,:),Fs);
        xhat1=r';
end

end

%% NMF function
function [W, H] = nmf(V, K, W, MAXITER, fixedInds)

F = size(V,1);
T = size(V,2);

rand('seed',0)
if isempty(W)
    W = 1+rand(F, sum(K));
end
H = 1+rand(sum(K), T);

inds = setdiff(1:sum(K),fixedInds);
ONES = ones(F,T);

for i=1:MAXITER
    
    % update activations
    H = H .* (W'*( V./(W*H+eps))) ./ (W'*ONES);
    
    % update dictionaries
    W(:,inds) = W(:,inds) .* ((V./(W*H+eps))*H(inds,:)') ./(ONES*H(inds,:)');
end

% normalize W to sum to 1
sumW = sum(W);
W = W*diag(1./sumW);
H = diag(sumW)*H;
end

% function [W, H] = nmf(V, K, MAXITER)
%
% F = size(V,1);
% T = size(V,2);
%
% rand('seed',0)
% W = 1+rand(F, K);
% % W = W./repmat(sum(W),F,1);
% H = 1+rand(K, T);
%
% ONES = ones(F,T);
%
% for i=1:MAXITER
%     % update activations
%     H = H .* (W'*( V./(W*H+eps))) ./ (W'*ONES);
%
%     % update dictionaries
%     W = W .* ((V./(W*H+eps))*H') ./(ONES*H');
%
% end
%
%
% % normalize W to sum to 1
% sumW = sum(W);
% W = W*diag(1./sumW);
% H = diag(sumW)*H;
% end

function output=valid_r(V,r)
if r==0 %might need to be isempty instead
    output=false;
elseif r>=min(size(V))
    output=false;
else
    output=true;
end
end

function output=valid_beta(beta_loss)
if beta_loss==0
    output=true;
elseif beta_loss==1
    output=true;
elseif beta_loss==2
    output=true;
else
    output=false;
end
end


function [W, H]=nmf_general(V,r,beta_loss,max_iter,fixedInds,W)
if valid_r(V,r)==false
    disp('Invalid inner dimension R')
    return
elseif valid_beta(beta_loss)==false
    disp('Invalid beta loss' )
    return
end

% initialize base matrix W and activation matrix H
% W=rand(size(V,1),r);
% H=rand(r,size(V,2));
F = size(V,1);
T = size(V,2);
if isempty(W)
    W = 1+rand(F, sum(r));
end
H = 1+rand(sum(r), T);

inds = setdiff(1:sum(r),fixedInds);

WH=W*H;

for i=1:max_iter
    %H updates
    num=W'*(WH.^(beta_loss-2).*V);
    dom=W'*WH.^(beta_loss-1);
    
    H=H.*(num./dom);
    %WH=W*H;
    
    %num=(WH.^(beta_loss-2).*V)*H';
    %dom=WH.^(beta_loss-1)*H';
    
    %W=W.*(num./dom);
    
    num=(WH.^(beta_loss-2).*V)*H(inds,:)';
    dom=WH.^(beta_loss-1)*H(inds,:)';
    
    W(:,inds)=W(:,inds).*(num./dom);
    
    WH=W*H;
end
% normalize W to sum to 1
sumW = sum(W);
W = W*diag(1./sumW);
H = diag(sumW)*H;
end


function output=valid_sparsity(sparsity)
if sparsity<0
    output=false;
else
    output=true;
end
end



function [W,H]=matrix_normalization(W,H)
F=size(W,1);
T=size(H,2);
Wn=ones(F)*W.^2;
Hn=ones(T,F)*W.^2;

W=W./Wn.^0.5;
H=H.*(Hn.^0.5)';
end

function [W,H]=sparse_nmf(V,r,sparsity,beta_loss,max_iter,fixedInds,W)
if valid_r(V,r)==false
    disp('Invalid inner dimension R')
    return
elseif valid_beta(beta_loss)==false
    disp('Invalid beta loss' )
    return
elseif valid_sparsity(sparsity)==false
    disp('Invalid sparsity')
end


F = size(V,1);
T = size(V,2);

rand('seed',0)
if isempty(W)
    W = 1+rand(F, sum(r));
end
H = 1+rand(sum(r), T);

inds = setdiff(1:sum(r),fixedInds);

% initialize base matrix W and activation matrix H
%W=rand(size(V,1),r); %need to confirm this type of random
%H=rand(r,size(V,2));
[W,H]=matrix_normalization(W,H);
WH=W*H;
for i=1:max_iter
    %H updates
    num=W'*(WH.^(beta_loss-2).*V);
    dom=W'*WH.^(beta_loss-1)+sparsity;
    
    H=H.*(num./dom);
    %WH=W*H;
    
    % W update
    %     F=size(W,1);
    %     temp=WH.^(beta_loss-1)*H'.*W;
    %     temp=ones(F)*temp;
    %     temp=temp.*W;
    %     negdp=(WH.^(beta_loss-2).*V)*H'+temp;
    %
    %     temp=(WH.^(beta_loss-2).*V)*H'.*W;
    %     temp=ones(F)*temp;
    %     temp=temp.*W;
    %     posdp=WH.^(beta_loss-1)*H'+temp;
    %
    %     W=W.*(negdp./posdp);
    F=size(W,1);
    temp=WH.^(beta_loss-1)*H(inds,:)'.*W(:,inds);
    temp=ones(F)*temp;
    temp=temp.*W(:,inds);
    negdp=(WH.^(beta_loss-2).*V)*H(inds,:)'+temp;
    
    temp=(WH.^(beta_loss-2).*V)*H(inds,:)'.*W(:,inds);
    temp=ones(F)*temp;
    temp=temp.*W(:,inds);
    posdp=WH.^(beta_loss-1)*H(inds,:)'+temp;
    
    W(:,inds)=W(:,inds).*(negdp./posdp);
    
    [W,H]=matrix_normalization(W,H);
    WH=W*H;
end
end

%% NMF with different optimizer does not work
function [W, H]=nmf_optimized(V,K,W,fixedInds)
F = size(V,1);
T = size(V,2);

rand('seed',0)
if isempty(W)
    W = 1+rand(F, sum(K));
end

H = 1+rand(sum(K), T);
inds = setdiff(1:sum(K),fixedInds);

Wfixed=W(:,fixedInds); 

Wvary=W(:,inds);
sizeWvary=size(Wvary); 
totalWvary=numel(Wvary);
sizeH=size(H);
totalH=numel(H); 
x0=[reshape(Wvary,[1,totalWvary]) reshape(H,[1,totalH])]; 
lb=zeros(size(x0));
options = optimoptions(@fmincon,'MaxIterations',1); 
[x,fval,exitflag,output] = fmincon(@(x)(myfun(x,V,totalWvary,sizeWvary,sizeH,Wfixed)),x0,[],[],[],[],lb,[],[],options); 

Wvary=reshape(x(1:totalWvary),sizeWvary);
H=reshape(x(totalWvary+1:end),sizeH);
W=[Wfixed Wvary];

end

function f=myfun(x,V,totalWvary,sizeWvary,sizeH,Wfixed)

Wvary=reshape(x(1:totalWvary),sizeWvary);
H=reshape(x(totalWvary+1:end),sizeH);
W=[Wfixed Wvary]; 

f=0.5*sum(sum((V-W*H).^2)); 

end 


%% NMCF function
function [W, H] = nmcf(V,V1, K, MAXITER,fixedInds,lambda)
F = size(V,1);
T = size(V,2);

F2 = size(V1,1);
T2 = size(V1,2);

rand('seed',0)
W = 1+rand(F, sum(K));

H = 1+rand(sum(K), T);

inds = setdiff(1:sum(K),fixedInds);
ONES = ones(F,T);
ONES_fixed=ones(F2,T2);

% for i=1:MAXITER
%     % update activations
%     %H = H .* (W'*( V./(W*H+eps))) ./ (W'*ONES);
%     H = H .* (W'*( V./(W*H+eps))) ./ (W'*ONES); 
%     H(fixedInds,:) =lambda*(H(fixedInds,:)  .* (W(:,fixedInds)'*( V1./(W*H+eps))) ./ (W(:,fixedInds)'*ONES_fixed));
%     
%     % update dictionaries
%     %W(:,inds) = W(:,inds) .* ((V./(W*H+eps))*H(inds,:)') ./(ONES*H(inds,:)');
%     %W(:,inds) = W(:,inds) .* ((V./(W*H+eps))*H(inds,:)') ./(ONES*H(inds,:)');
%     W = W .* ((V./(W*H+eps))*H') ./(ONES*H');
%     W(:,fixedInds) = lambda*(W(:,fixedInds) .* ((V1./(W*H+eps))*H(fixedInds,:)') ./(ONES_fixed*H(fixedInds,:)'));  
% end

for i=1:MAXITER
    % update activations
    %H = H .* (W'*( V./(W*H+eps))) ./ (W'*ONES);
    H = H .* (W'*( V./(W*H+eps))) ./ (W'*ONES); 
    H(fixedInds,:) =H(fixedInds,:)+lambda*(H(fixedInds,:)  .* (W(:,fixedInds)'*( V1./(W*H+eps))) ./ (W(:,fixedInds)'*ONES_fixed)-H(fixedInds,:));
    
    % update dictionaries
    %W(:,inds) = W(:,inds) .* ((V./(W*H+eps))*H(inds,:)') ./(ONES*H(inds,:)');
    %W(:,inds) = W(:,inds) .* ((V./(W*H+eps))*H(inds,:)') ./(ONES*H(inds,:)');
    W = W .* ((V./(W*H+eps))*H') ./(ONES*H');
    W(:,fixedInds) = W(:,fixedInds)+lambda*(W(:,fixedInds) .* ((V1./(W*H+eps))*H(fixedInds,:)') ./(ONES_fixed*H(fixedInds,:)')-W(:,fixedInds));  
end

% normalize W to sum to 1
sumW = sum(W);
W = W*diag(1./sumW);
H = diag(sumW)*H;
end





