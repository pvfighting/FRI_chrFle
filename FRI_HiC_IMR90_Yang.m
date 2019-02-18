clc;
clear;

%%  cHIP_seq data  for chr17 all seems to be good!
%dA=importdata('./ATAC-seq/GM12878_ATACseq_50k_AllReps_ZINBA_pp08.bed');
dA=importdata('./Dnase/GM12878_ENCSR000EJD_rep1_1_rep1_2_rep1_3_se_bwa_biorep_filtered_peaks.bed');

tic;
for chr_num = 20:20
if chr_num < 23    
 chr=strcat('chr',int2str(chr_num));
else
 chr='chrX';
end

%chr='chr19';

res_str='5kb';

nres=5000;
span = 2e5/nres;  

% It seems results for MAPQG0 is better than MAPQGE30!! Need to check out! 

% chr_hic=strcat('./IMR90/',res_str,'_resolution_intrachromosomal/',chr,'/MAPQGE30/',chr,'_',res_str,'.RAWobserved');
% chr_vcnorm=strcat('./IMR90/',res_str,'_resolution_intrachromosomal/',chr,'/MAPQGE30/',chr,'_',res_str,'.VCnorm');

chr_hic=strcat('./GM12878/',res_str,'_resolution_intrachromosomal/',chr,'/MAPQGE30/',chr,'_',res_str,'.RAWobserved');
chr_vcnorm=strcat('./GM12878/',res_str,'_resolution_intrachromosomal/',chr,'/MAPQGE30/',chr,'_',res_str,'.VCnorm');

%% Hi-C data
vdiag=load(chr_hic);
vnm=load(chr_vcnorm);

L=max(vdiag(:,2)/nres)+1;

uij_raw=zeros(L,L);
uij=zeros(L,L);
[lr, lc]=size(vdiag);
for it=1:lr
    i=vdiag(it,1)/nres+1;
    j=vdiag(it,2)/nres+1;
    uij_raw(i,j)=vdiag(it,3);
    uij_raw(j,i)=vdiag(it,3);
end

for i=1:L
    for j=1:L
        vcscale=vnm(i)*vnm(j);
        if(abs(vcscale)>1e-10)
            uij(i,j)=uij_raw(i,j)/vcscale;
        end
    end
end

%%   The process used in GNM code, seems not good as vnm_mat terms can be very close to zero. 

%     vnm = vnm(1:L);
%     vnm_mat = vnm * vnm';
%     mat = uij_raw./ vnm_mat;
%     mat(isnan(mat)) = 0;
%     mat(isinf(mat)) = 0;
%     uij=mat;

%% Same process as in GNM paper (rewrite to be more clear)
bed=str2double(dA.textdata(strcmp(dA.textdata,chr),2:3));
signal=dA.data(strcmp(dA.textdata,chr));  %%%%%%%%%%% not use count, instead use the signal information

startIDX = ceil((bed(:,1)+1)/nres);  
endIDX = ceil(bed(:,2)/nres);

series = zeros(1, max(endIDX));
%count = zeros(1, max(endIDX));
for i = 1:length(startIDX)
    series(startIDX(i)) = series(startIDX(i)) + signal(i);
   % count(startIDX(i)) = count(startIDX(i)) + 1;
end

series = series/norm(series);    %%%%%%%%%% signal data are first normalized  
dnasevec_org=smooth(series, span);        %%%%%%%%%% signal data are then averaged by smooth function, span=2e5/nres; 

dnasevec=zeros(L,1);
ndnas=min(L,size(dnasevec_org));
dnasevec(1:ndnas)=dnasevec_org(1:ndnas);

%% FRI model

idx_hic = find(diag(uij) ~= 0);  %%% select nonzero regions in Hi-C
Hg = uij(idx_hic,idx_hic);
[~, nL]=size(Hg);

 sigma=0.0001;
 np0=1;
%for sigma=0.00010:0.00010:0.0005
    %for np0=0.6:0.1:1.4
hfri=zeros(nL,1);
for i=1:nL
   temp=0; 
   for j=1:nL
       if(j~=i && Hg(i,j)>0)
         dis=(1/Hg(i,j));   
        % temp=temp+exp(-(dis/sigma)^np0);
        temp=temp+1/(1+(dis/sigma)^np0);
       end
   end
   hfri(i)=1/temp;
end

flex=zeros(L,1);    %% Put the zero terms back to the flexibility region (Obviously not very reasonable) 
                    %  But Spearman correlation value enlarges. However, PCC is not enlarged
                     
flex(idx_hic)=hfri;

mutlen = min([L, length(dnasevec)]);  %%%%%%  Hi-C, DNase and ATAC may varies about +/-1 in their length 

x=flex(1:mutlen);
y=dnasevec(1:mutlen);


tt=corr(x,y,'type', 'Spearman');
fprintf('%d chromosom %f \n',chr_num,tt);

p=polyfit(x,y,1);
plot(y);hold on;
ylim([0,0.02]);
xlabel('Locus index (5kb)')
ylabel('Normalized DNase')
xlim([0,13000]);
%plot(x,'g','linewidth',1.5);hold on;
plot(p(1)*x+p(2)); hold on;
set(legend('DNase','FRI'),'FontSize',12);
set(text(150,0.019,'Chromosome 20'),'FontSize',12);
%tt1=corr(hfri,dnasevec(idx_hic),'type', 'Spearman');  %% Spearman correlation value for nonzero flexibility regions. 
end
disp(['toc:',num2str(toc)]);

