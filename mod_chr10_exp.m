
clear;
clc;

%---------------------------------------------------------------------%
%--------Code for performing normal mode analysis using aFRI----------%
%---The conformations obtained from GEM code is loaded as atom data---%
%---------------------------------------------------------------------%

atom =load('conformation1.txt');
atom(:,1)=atom(:,1)/50; %Scale the data by 50.
atom(:,2)=atom(:,2)/50;
atom(:,3)=atom(:,3)/50;

loc=atom;
[num, ntt]=size(atom);
hess = zeros(3*num,3*num);
disData=zeros(num,num);
eta=26;   
eta2=eta*eta;
%---------------------------------------------------------------------%
%Construction of the Hessian Matrix for EVD
%---------------------------------------------------------------------%

for i=1:num
    ix1=loc(i,1);
    iy1=loc(i,2);
    iz1=loc(i,3);
    
    mu_xxj= 0.0; 
    mu_yyj= 0.0; 
    mu_zzj= 0.0; 

    mu_xyj= 0.0;
    mu_xzj= 0.0;
    mu_yzj= 0.0;
  
    for j=1:num
    hess3 = zeros(3,3);    
    ix2=loc(j,1);   
    iy2=loc(j,2);     
    iz2=loc(j,3);
    
    bx=ix1-ix2;
    by=iy1-iy2;
    bz=iz1-iz2;     

    dis2=bx^2+by^2+bz^2;
    dist=sqrt(dis2);
    disData(i,j)=dist;

        if(i~=j) 
        mu_xxj= (2*exp(-(dis2)/eta2))/eta2 - (exp(-(dis2)/eta2)*4*bx^2)/eta2^2;
        mu_yyj= (2*exp(-(dis2)/eta2))/eta2 - (exp(-(dis2)/eta2)*4*by^2)/eta2^2;
        mu_zzj= (2*exp(-(dis2)/eta2))/eta2 - (exp(-(dis2)/eta2)*4*bz^2)/eta2^2;
       
        mu_xyj= -(exp(-(dis2)/eta2)*4*bx*by)/eta2^2;
        mu_xzj= -(exp(-(dis2)/eta2)*4*bx*bz)/eta2^2;
        mu_yzj= -(exp(-(dis2)/eta2)*4*by*bz)/eta2^2;

        hess3(1,1)=(mu_yyj*mu_zzj-mu_yzj*mu_yzj);
        hess3(1,2)=(mu_xzj*mu_yzj-mu_xyj*mu_zzj);
        hess3(1,3)=(mu_xyj*mu_yzj-mu_xzj*mu_yyj);

        hess3(2,1)=(mu_yzj*mu_xzj-mu_xyj*mu_zzj);
        hess3(2,2)=(mu_xxj*mu_zzj-mu_xzj*mu_xzj);
        hess3(2,3)=(mu_xzj*mu_xyj-mu_xxj*mu_yzj);

        hess3(3,1)=(mu_xyj*mu_yzj-mu_yyj*mu_xzj);
        hess3(3,2)=(mu_xyj*mu_xzj-mu_xxj*mu_yzj);
        hess3(3,3)=(mu_xxj*mu_yyj-mu_xyj*mu_xyj); 


        hess(3*(i-1)+1,3*(j-1)+1)=hess3(1,1);
        hess(3*(i-1)+1,3*(j-1)+2)=hess3(1,2);
        hess(3*(i-1)+1,3*(j-1)+3)=hess3(1,3);

        hess(3*(i-1)+2,3*(j-1)+1)=hess3(2,1);
        hess(3*(i-1)+2,3*(j-1)+2)=hess3(2,2);
        hess(3*(i-1)+2,3*(j-1)+3)=hess3(2,3);

        hess(3*(i-1)+3,3*(j-1)+1)=hess3(3,1);
        hess(3*(i-1)+3,3*(j-1)+2)=hess3(3,2);
        hess(3*(i-1)+3,3*(j-1)+3)=hess3(3,3);

        hess(3*(i-1)+1,3*(i-1)+1) = hess(3*(i-1)+1,3*(i-1)+1)-hess3(1,1);
        hess(3*(i-1)+1,3*(i-1)+2) = hess(3*(i-1)+1,3*(i-1)+2)-hess3(1,2);
        hess(3*(i-1)+1,3*(i-1)+3) = hess(3*(i-1)+1,3*(i-1)+3)-hess3(1,3);

        hess(3*(i-1)+2,3*(i-1)+1) = hess(3*(i-1)+2,3*(i-1)+1)-hess3(2,1);
        hess(3*(i-1)+2,3*(i-1)+2) = hess(3*(i-1)+2,3*(i-1)+2)-hess3(2,2);
        hess(3*(i-1)+2,3*(i-1)+3) = hess(3*(i-1)+2,3*(i-1)+3)-hess3(2,3);

        hess(3*(i-1)+3,3*(i-1)+1) = hess(3*(i-1)+3,3*(i-1)+1)-hess3(3,1);
        hess(3*(i-1)+3,3*(i-1)+2) = hess(3*(i-1)+3,3*(i-1)+2)-hess3(3,2);
        hess(3*(i-1)+3,3*(i-1)+3) = hess(3*(i-1)+3,3*(i-1)+3)-hess3(3,3);
        end
    end
end

% ------------------Hessian Matrix Decompostion ----------------------- %
% The singular value decomposition is used for matrix decomposition                         
% --------------------------------------------------------------------- %

[U,S,V] = svd(hess);
ireg=15;
nv=ireg; 
nm3=3*num;
vect=zeros(num,3,nv);
s_diag=zeros(3*num,1);

for i=1:nm3
    s_diag(i)=S(nm3+1-i,nm3+1-i); 
end

for i=1:num
    for nc=1:nv
      vect(i,1,nc)=U(3*i-2,nm3+1-nc);
      vect(i,2,nc)=U(3*i-1,nm3+1-nc);
      vect(i,3,nc)=U(3*i  ,nm3+1-nc);
    end
end

% File writing for visualising the modes in VMD
Iresids=zeros(num);  
chainnam=blanks(num);
for i=1:num   
   Iresids(i)=i;
   chainnam(i)='A';
end

fileID = fopen('aFRI_modes.nmd','w');
fprintf(fileID,'isovalue of FRI surface\n');
fprintf(fileID,'name EMD\n');
fprintf(fileID,'resids');
fprintf(fileID,'%5d',Iresids(1:num));
fprintf(fileID,'\n');
fprintf(fileID,'chainids');
for i=1:num
fprintf(fileID,'%5s',chainnam(i));
end
fprintf(fileID,'\n');
fprintf(fileID,'coordinates');
for i=1:num
fprintf(fileID,'%8.3f %8.3f %8.3f',loc(i,1:3));
end

nmd=5;
pm=[1 2 3 11 15];
for it=1:nmd
        fprintf(fileID,'\n');
        fprintf(fileID,'mode %5d',it);
        for i=1:num
        fprintf(fileID,'%8.3f %8.3f %8.3f',vect(i,1:3,pm(it)));
        end
end
fclose(fileID);

