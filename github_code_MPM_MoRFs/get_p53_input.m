function get_p53_input

%%%%%%%%%%%%%%%%%%%%%% p53 sequence %%%%%%%%%%%%%%%%%%%%%%%%%%%
seq1 = 'meepqsdpsvepplsqetfsdlwkllpennvlsplpsqamddlmlspddieqwftedpgp';
seq2 = 'deaprmpeaaprvapapaaptpaapapapswplsssvpsqktyqgsygfrlgflhsgtak';
seq3 = 'svtctyspalnkmfcqlaktcpvqlwvdstpppgtrvramaiykqsqhmtevvrrcphhe';
seq4 = 'rcsdsdglappqhlirvegnlrveylddrntfrhsvvvpyeppevgsdcttihynymcns';
seq5 = 'scmggmnrrpiltiitledssgnllgrnsfevrvcacpgrdrrtekenlrkkgephhelp';
seq6 = 'pgstkralpnntssspqpkkkpldgeyftlqirgrerfemfrelnealelkdaqagkepg';
seq7 = 'gsrahsshlkskkgqstsrhkklmfktegpdsd';
seq = [seq1 seq2 seq3 seq4 seq5 seq6 seq7];
seq = upper(seq);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%    p53_Y: N*1    %%%%%%%%%%%%%%%%%%%%%%
N = length(seq);
Yseq = zeros(N,1)-1;
yreg = [15 29; 33 60; 367 388]; % Èı¶ÎÇøÓòÎª1
for i = 1:size(yreg)
    Yseq(yreg(i,1):yreg(i,2)) = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%   p53_input  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load aaindexdata.mat  % head mapmatix %544*20
aalocs = [54    71   211   262   339   346   350   351   532   535 536   537   541];
mapmatix = mapmatix(aalocs,:);
wins = [10 45 90];
Xseqm = []; %
for i = 1:length(wins)
    win = wins(i);
    Xseqi = getdatamean(seq,win,mapmatix); %N * Nfea
    Xseqm = [Xseqm Xseqi];
end

save p53_Xseqm_16fs_3wins Xseqm Yseq seq

end




function Xseqi = getdatamean(seq,win,mapmatix)
Xseqi = [];  %N * Nfea
Ns = length(seq);
winadd = fix(win/2); 
Nindex = size(mapmatix,1);

seqadd = [zeros(1,winadd) seq zeros(1,winadd)];
seq_remark = [zeros(1,winadd) map_remark_seq(seq) zeros(1,winadd)];
seq_deleage = [zeros(1,winadd) map_deleage_seq(seq) zeros(1,winadd)];
aaout = aaindex_matix(seq,mapmatix);  %13*Ns 
seq_aaindex = [zeros(Nindex,winadd) aaout zeros(Nindex,winadd)] ; %13*(Ns+win-1)

N_features = 16; % feature number   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seq_result = zeros(N_features,Ns+win-1)   ;   % Store the cumulative sum of each residue when the window moves
seq_num = zeros(1,Ns+win-1);

for is =  1:Ns
    seq_aa = seqadd(is:is+win-1); %the region of the is-th window
    seq_num(is:is+win-1) = seq_num(is:is+win-1) + 1; %Add up the number of times each residue moves in the window
    seq380 = map_380(seq_aa);
    d=2;
    TK1 = mean(seq_remark(is:is+win-1));
    TK2 = mean(seq_deleage(is:is+win-1));
    TK3 = TopoEn(seq380,d);
    TKaaindex = mean(seq_aaindex(:,is:is+win-1)'); %1*13
    TK = [TK1 TK2 TK3 TKaaindex];
    nk = length(TK);
    for fi = 1:nk
        seq_result(fi,is:is+win-1) = seq_result(fi,is:is+win-1) + TK(fi);
        % Accumulate the sum of each feature when the window moves
    end
end
TK_matrix = seq_result;
for fi = 1:size(TK_matrix,1)
    TK_matrix(fi,:) = TK_matrix(fi,:)./seq_num;  
end
for is = 1:Ns
    Xseqi(is,:) = TK_matrix(:,is+winadd)';  
end

end

function aaout = aaindex_matix(data,mapmatix)
%mapmatix -- Overall AAindex mapping,544*20
%aaout -- The matrix of the mapped sequence, 544*Ndata
[Nindex Nele] = size(mapmatix);   
map = zeros(Nindex,Nele);
for aai = 1:Nindex
    map(aai,:) = mapmatix(aai,:)./max(abs(mapmatix(aai,:)));  %-1~1
end

Ndata = length(data);
aaout = zeros(Nindex,Ndata); 
for n = 1:Ndata
    switch data(n)
        case 'A'
            aaout(:,n) = map(:,1);
        case 'R'
            aaout(:,n) = map(:,2);
        case 'N'
            aaout(:,n) = map(:,3);
        case 'D'
            aaout(:,n) = map(:,4);
        case 'C'
            aaout(:,n) = map(:,5);
        case 'Q'
            aaout(:,n) = map(:,6);
        case 'E'
            aaout(:,n) = map(:,7);
        case 'G'
            aaout(:,n) = map(:,8);
        case 'H'
            aaout(:,n) = map(:,9);
        case 'I'
            aaout(:,n) = map(:,10);
        case 'L'
            aaout(:,n) = map(:,11);
        case 'K'
            aaout(:,n) = map(:,12);
        case 'M'
            aaout(:,n) = map(:,13);
        case 'F'
            aaout(:,n) = map(:,14);
        case 'P'
            aaout(:,n) = map(:,15);
        case 'S'
            aaout(:,n) = map(:,16);
        case 'T'
            aaout(:,n) = map(:,17);
        case 'W'
            aaout(:,n) = map(:,18);
        case 'Y'
            aaout(:,n) = map(:,19);
        case 'V'
            aaout(:,n) = map(:,20);
        otherwise
            aaout(:,n) = 0;
    end
    
end

end

function topo = TopoEn(seq,d)
%If the sequence length is larger than sup, take 1 as the step size, and take the window successively
Ns = length(seq); 
n0 = fix(log2(Ns)/log2(d))-1; %Calculates the maximum possible word length that the sequence length can count
%Using the bottom formula, calculate logd(Ns)
sup = d^n0 + n0 - 1; %Compute the minimum sequence length required to have a word length of n0
lin = d^(n0+1) + (n0+1) - 1; % Calculates the longest sequence length of the word length n0
while(Ns > lin )
    n0 = n0 + 1;
    sup = d^n0 + n0 - 1; 
    lin = d^(n0+1) + (n0+1) - 1; 
end
while(Ns < sup)
    n0 = n0 - 1;
    sup = d^n0 + n0 - 1; 
    lin = d^(n0+1) + (n0+1) - 1; 
end
wadd = Ns - sup + 1; % Total number of sliding windows
topo = 0;
for wi = 1:wadd
    mnum = seq(wi:n0+wi-1); %Stores the subsequence of the word length n0 that appears in the sequence
    for ns = 2+wi-1:sup+wi-1-(n0-1)
        nsub = seq(ns:ns+n0-1);   %A subsequence waiting for alignment
        Nmnum = size(mnum,1); %The number of columns in a subsequence matrix
        str_result = 0; %cumulative
        for nmn = 1:Nmnum
            nmnum = mnum(nmn,:);  %Subsequences that have already appeared
            midresult = diff(nsub,nmnum);
            str_result = str_result + midresult;
            %Compare the subsequence waiting for alignment with the labeled subsequence, equal to 1, different to 0, and accumulate
        end
        if str_result == 0
            mnum(Nmnum+1,:) = nsub;
        end
    end
    pwn = size(mnum,1); 
    topo = topo + (log2(pwn)/log2(d))/n0;
end
topo = topo/wadd;
end

function result = diff(a1,a2)
%same -- 1£¬different -- 0
N1 = length(a1);
aresult = 0;
for n = 1:N1      
    if a1(n) == a2(n)        
        aresult = aresult +1;        
    end    
end
if aresult == N1   
    result = 1;  % two vectors are exactly the same
else
    result = 0;
end
end


function   map_m = map_380(data_m)
hyd_unit = [0 0 0 0 0 0 0 0 0 1 1 0 0 1 0 0 0 1 1 1]; 
N = length(data_m);
map_m = [];
for i = 1:N
    switch data_m(i)
        case 'A'
            map_m(i) = hyd_unit(1);
        case 'R'
            map_m(i) = hyd_unit(2);
        case 'N'
            map_m(i) = hyd_unit(3);
        case 'D'
            map_m(i) = hyd_unit(4);
        case 'C'
            map_m(i) = hyd_unit(5);
        case 'Q'
            map_m(i) = hyd_unit(6);
        case 'E'
            map_m(i) = hyd_unit(7);
        case 'G'
            map_m(i) = hyd_unit(8);
        case 'H'
            map_m(i) = hyd_unit(9);
        case 'I'
            map_m(i) = hyd_unit(10);
        case 'L'
            map_m(i) = hyd_unit(11);
        case 'K'
            map_m(i) = hyd_unit(12);
        case 'M'
            map_m(i) = hyd_unit(13);
        case 'F'
            map_m(i) = hyd_unit(14);
        case 'P'
            map_m(i) = hyd_unit(15);
        case 'S'
            map_m(i) = hyd_unit(16);
        case 'T'
            map_m(i) = hyd_unit(17);
        case 'W'
            map_m(i) = hyd_unit(18);
        case 'Y'
            map_m(i) = hyd_unit(19);
        case 'V'
            map_m(i) = hyd_unit(20);
        otherwise
            map_m(i) = 0;
    end
end

end

function   map_m = map_remark_seq(data_m)
hyd_unit = [0.1739 -0.0537 -0.2141 0.2911 -0.5301 0.3088 0.5214 0.0149 0.1696 -0.2907 -0.3379 ...
            0.1984 -0.1113 -0.8434 -0.0558 0.2627 -0.1297 -1.3710 -0.8040 -0.2405]; 
N = length(data_m);
map_m = [];
for i = 1:N
    switch data_m(i)
        case 'A'
            map_m(i) = hyd_unit(1);
        case 'R'
            map_m(i) = hyd_unit(2);
        case 'N'
            map_m(i) = hyd_unit(3);
        case 'D'
            map_m(i) = hyd_unit(4);
        case 'C'
            map_m(i) = hyd_unit(5);
        case 'Q'
            map_m(i) = hyd_unit(6);
        case 'E'
            map_m(i) = hyd_unit(7);
        case 'G'
            map_m(i) = hyd_unit(8);
        case 'H'
            map_m(i) = hyd_unit(9);
        case 'I'
            map_m(i) = hyd_unit(10);
        case 'L'
            map_m(i) = hyd_unit(11);
        case 'K'
            map_m(i) = hyd_unit(12);
        case 'M'
            map_m(i) = hyd_unit(13);
        case 'F'
            map_m(i) = hyd_unit(14);
        case 'P'
            map_m(i) = hyd_unit(15);
        case 'S'
            map_m(i) = hyd_unit(16);
        case 'T'
            map_m(i) = hyd_unit(17);
        case 'W'
            map_m(i) = hyd_unit(18);
        case 'Y'
            map_m(i) = hyd_unit(19);
        case 'V'
            map_m(i) = hyd_unit(20);
        otherwise
            map_m(i) = 0;
    end
end
end

function   map_m = map_deleage_seq(data_m)
hyd_unit = [-0.275 -0.179 0.479 0.4645 -0.1255 -0.055 -0.2745 0.6675 0.135 -0.515 -0.4385 ...
            -0.0495 -0.4765 -0.497 1.117 0.2965 0.145 -0.257 0.0825 -0.7055]; 
N = length(data_m);
map_m = [];
for i = 1:N
    switch data_m(i)
        case 'A'
            map_m(i) = hyd_unit(1);
        case 'R'
            map_m(i) = hyd_unit(2);
        case 'N'
            map_m(i) = hyd_unit(3);
        case 'D'
            map_m(i) = hyd_unit(4);
        case 'C'
            map_m(i) = hyd_unit(5);
        case 'Q'
            map_m(i) = hyd_unit(6);
        case 'E'
            map_m(i) = hyd_unit(7);
        case 'G'
            map_m(i) = hyd_unit(8);
        case 'H'
            map_m(i) = hyd_unit(9);
        case 'I'
            map_m(i) = hyd_unit(10);
        case 'L'
            map_m(i) = hyd_unit(11);
        case 'K'
            map_m(i) = hyd_unit(12);
        case 'M'
            map_m(i) = hyd_unit(13);
        case 'F'
            map_m(i) = hyd_unit(14);
        case 'P'
            map_m(i) = hyd_unit(15);
        case 'S'
            map_m(i) = hyd_unit(16);
        case 'T'
            map_m(i) = hyd_unit(17);
        case 'W'
            map_m(i) = hyd_unit(18);
        case 'Y'
            map_m(i) = hyd_unit(19);
        case 'V'
            map_m(i) = hyd_unit(20);
        otherwise
            map_m(i) = 0;
    end
end
end



