function predict_protein_p53

load p53_Xseqm_16fs_3wins % Xseqm Yseq seq

load predict_MPM_MoRFs_weight_16fs_3wins % a b prenormal threds maxs mins 
thred = 0  %threds(3)
Xseqm = getminmaxnor(Xseqm,maxs,mins);
PreRe = getPre(Xseqm,a,b,thred,prenormal);

% figure the prediction results
plot(1:length(seq),PreRe,'b-'); hold on
plot(1:length(seq),zeros(1,length(seq)),'k');

% figure the MoRF regions
yreg = [15 29; 33 60; 367 388]; % three regions -- 1
ymax = 0.35;
for i = 1:3
    ni = yreg(i,2)-yreg(i,1) +1;
    plot(yreg(i,1):yreg(i,2),ones(1,ni).*ymax,'r-'); hold on
end
yreg_x = yreg(:);
for i = 1:6
    ni = 100;
    step = ymax*2/ni;
plot(ones(1,ni+1).*yreg_x(i),-ymax:step:ymax,'r-'); hold on
end
xlabel('the residue index');
ylabel('prediction result');   

% fix TPR, calculate FPR and ACC
[AUC_train mapo]   = getAUC_MPM(Xseqm,Yseq,a,b);
TPR_taget = [0.222 0.254 0.389];
[TPR_get FPR_get ACC_get] = getTPRFPR(mapo,TPR_taget);
disp('test_[TPR_taget; TPR_get; FPR_get; ACC_get]: ')
disp([TPR_taget;TPR_get;FPR_get;ACC_get])

end


function PreRe = getPre(data,a,b,thred,prenormal)
n = size(data,1); % the number of samples 
pre_initial = data*a - b*ones(n,1);  
pre_scale = pre_initial ./ prenormal;  
PreRe = pre_scale + thred;

end

function Xm = getminmaxnor(X,maxs,mins)
% min-max
% X -- Ntr*Nfeas
[Ntr Nfeas] = size(X);
Xm = zeros(Ntr,Nfeas);
for i = 1:Nfeas
    Xi = X(:,i);
    mini = mins(i);
    maxi = maxs(i);
    Xmi = (Xi - mini)./(maxi-mini);
    Xm(:,i) = Xmi;
end
end

function [AUC maporder] = getAUC_MPM(xdata,ydata,a,b)
% xdata : Ntr*Nf;  w: Nf*1;
% ydata : Ntr*1  % disÎª1£»ordÎª-1
N = length(ydata);
mapdata = xdata*a  - b*ones(N,1);  
% min-max,[0,1]
epsio = 10^(-6); % add epsio, (0,1)
mapmin = min(mapdata) - epsio;
mapmax = max(mapdata) + epsio;
mapunit = (mapdata - mapmin)./(mapmax-mapmin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mapmatix = [mapunit ydata]; %Ntr*2
maporder = sortrows(mapmatix,1); % Sort in ascending order in the first column

loc_positive = find(maporder(:,2) == 1);
N_p = length(loc_positive); %the number of positive samples
N_n = length(ydata) - N_p;  %the number of negtive samples
AUC = (sum(loc_positive) - (N_p*(N_p+1)/2)) / (N_p*N_n);

end

function [TPR_get FPR_get ACC_get] = getTPRFPR(maporder,TPR_taget)
[TPR FPR ACC] = getdiffTPR(maporder);
N_taget = length(TPR_taget);
TPR_get = zeros(1,N_taget);
FPR_get = zeros(1,N_taget);
ACC_get = zeros(1,N_taget);
for i = 1:N_taget
    tageti = TPR_taget(i);
    [val loc] = min(abs(TPR-tageti));
    TPR_get(i) = TPR(loc);
    FPR_get(i) = FPR(loc);
    ACC_get(i) = ACC(loc);
end
end

function [TPR FPR ACC] = getdiffTPR(maporder)
% calculate all values of TPR,FPR,ACC
loc_p = find(maporder(:,2) == 1);
loc_n = find(maporder(:,2) == -1);
N_p = length(loc_p); 
N_n = length(loc_n); 

scores = maporder(:,1);
TPR = [];
FPR = [];
ACC = [];
for i = 1:N_p
    thredloci = loc_p(i);
    thredi = scores(thredloci);
    [TPR(i,1) FPR(i,1) ACC(i,1)] = calthred(scores,thredi,loc_p,loc_n,N_p,N_n);
end
end

function [TPR FPR ACC] = calthred(scores,thred,loc_p,loc_n,N_p,N_n)
ytest = scores > thred;  

TP = sum(ytest(loc_p)==1);
TN = sum(ytest(loc_n)==0);
%FN = N_p - TP;
FP = N_n - TN;

TPR = TP/N_p;
FPR = FP/N_n;
ACC = (TP+TN)/(N_p+N_n);
end

