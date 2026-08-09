function [q, h]=threeD_datum()

% mosaic_mask='L';
Pline='1239';
Eline='0969';

dir = strcat('P:\s\japan\TIFF\ALOStiff_401\DA', Pline);
dir = strcat(dir,'.tif');
[PALSAR,PALSARMAP] = geotiffread(dir);
dir = strcat('P:\s\japan\TIFF\ALOStiff_401\M', Pline);
dir = strcat(dir,'.tif');
[MPALSAR] = geotiffread(dir);

dir = strcat('P:\s\japan\TIFF\ALOStiff_401\DASAR', Eline);
dir = strcat(dir,'.tif');
[ENVISAT,ENVISATMAP] = geotiffread(dir);
dir = strcat('P:\s\japan\TIFF\ALOStiff_401\M', Eline);
dir = strcat(dir,'.tif');
[MENVISAT] = geotiffread(dir);

dir = strcat('P:\s\japan\TIFF\ALOStiff_401\ILOS', Pline);
dir = strcat(dir,'.tif');
[PALSARILOS] = geotiffread(dir);

dir = strcat('P:\s\japan\TIFF\ALOStiff_401\ILOS', Eline);
dir = strcat(dir,'.tif');
[ENVISATILOS] = geotiffread(dir);

dir = strcat('P:\s\japan\TIFF\ALOStiff_401\precision', Pline);
dir = strcat(dir,'.tif');
[PALSARPRECISION] = geotiffread(dir);
dir = strcat('P:\s\japan\TIFF\ALOStiff_401\precisionmai', Pline);
dir = strcat(dir,'.tif');
[PALSARMAIPRECISION] = geotiffread(dir);

dir = strcat('P:\s\japan\TIFF\ALOStiff_401\precision', Eline);
dir = strcat(dir,'.tif');
[ENVISATPRECISION] = geotiffread(dir);
dir = strcat('P:\s\japan\TIFF\ALOStiff_401\precisionmai', Eline);
dir = strcat(dir,'.tif');
[ENVISATMAIPRECISION] = geotiffread(dir);

k=1;
for i=1:20:PALSARMAP.RasterSize(1,1)
    for j=1:20:PALSARMAP.RasterSize(1,2)
        if ~isnan(PALSAR(i,j)) & ~isnan(MPALSAR(i,j))
            SNP(k,1)=((i-1)*PALSARMAP.RasterSize(1,2)+j);
            k=k+1;
        end   
    end
end

pos=zeros(0,0);
prc=zeros(0,0);
arc=zeros(0,0);
k=1;
for i=1:size(SNP,1)
       lat=0;
       lon=0;
       i_r=floor(SNP(i)/PALSARMAP.RasterSize(1,2));
       res=SNP(i)-i_r*PALSARMAP.RasterSize(1,2);
       lat=PALSARMAP.Latlim(1,2)+PALSARMAP.DeltaLat*i_r;
       if 0~=res
            lat=lat+PALSARMAP.DeltaLat/2;
            lon=res*PALSARMAP.DeltaLon+PALSARMAP.Lonlim(1,1)-PALSARMAP.DeltaLon/2;
       else
            lat=lat-PALSARMAP.DeltaLat/2;
            lon=PALSARMAP.Lonlim(1,2)-PALSARMAP.DeltaLon/2;
       end
       
       if lat<=ENVISATMAP.Latlim(1,2) & lat>=ENVISATMAP.Latlim(1,1) & lon<=ENVISATMAP.Lonlim(1,2) & lon>=ENVISATMAP.Lonlim(1,1)  
            Z=floor((lat-ENVISATMAP.Latlim(1,2))/ENVISATMAP.DeltaLat)+1;
            U=floor((lon-ENVISATMAP.Lonlim(1,1))/ENVISATMAP.DeltaLon)+1;
            
            if ~isnan(ENVISAT(Z,U)) & ~isnan(MENVISAT(Z,U))
                            
                            i_r=floor(SNP(i)/PALSARMAP.RasterSize(1,2));
                            res=SNP(i)-i_r*PALSARMAP.RasterSize(1,2);
                            if 0~=res
                                row=i_r+1;
                                col=res;
                            elseif res==0
                                row=i_r;
                                col=PALSARMAP.RasterSize(1,2);
                            end
                            
                            prc(k,1)=row;prc(k,2)=col;
                            arc(k,1)=Z;arc(k,2)=U;
                            pos(k,1)=lat;pos(k,2)=lon;
                            k=k+1;                           
            end        
       end      
end
clear SNP PALSARMAP ENVISATMAP;
Ya=zeros(size(prc,1),1);
Aa=zeros(size(prc,1),3*size(prc,1));
Pa=zeros(size(prc,1),1);
for k=1:size(prc,1)
   Aa(k,3*(k-1)+1)=-sin(PALSARILOS(prc(k,1),prc(k,2))/180*pi)*cos(35/18*pi);
   Aa(k,3*(k-1)+2)=sin(PALSARILOS(prc(k,1),prc(k,2))/180*pi)*sin(35/18*pi);
   Aa(k,3*(k-1)+3)=cos(PALSARILOS(prc(k,1),prc(k,2))/180*pi);
   Ya(k,1)=PALSAR(prc(k,1),prc(k,2));
   % Pa(k)=1/(PALSARPRECISION(prc(k,1),prc(k,2))*PALSARPRECISION(prc(k,1),prc(k,2)));
   Pa(k)=1/0.00001139188;
end
clear PALSAR PALSARILOS PALSARPRECISION;
Yd=zeros(size(prc,1),1);
Ad=zeros(size(prc,1),3*size(prc,1));
Pd=zeros(size(prc,1),1);
for k=1:size(prc,1)
   Ad(k,3*(k-1)+1)=-sin(ENVISATILOS(arc(k,1),arc(k,2))/180*pi)*cos(19/18*pi);
   Ad(k,3*(k-1)+2)=sin(ENVISATILOS(arc(k,1),arc(k,2))/180*pi)*sin(19/18*pi);
   Ad(k,3*(k-1)+3)=cos(ENVISATILOS(arc(k,1),arc(k,2))/180*pi);                                              
   Yd(k,1)=ENVISAT(arc(k,1),arc(k,2));
   % Pd(k)=1/(ENVISATPRECISION(arc(k,1),arc(k,2))*ENVISATPRECISION(arc(k,1),arc(k,2)));
   Pd(k)=1/0.00000465918;
end
clear ENVISAT ENVISATILOS ENVISATPRECISION;
Yma=zeros(size(prc,1),1);
Ama=zeros(size(prc,1),3*size(prc,1));
Pma=zeros(size(prc,1),1);
for k=1:size(prc,1)
   Ama(k,3*(k-1)+1)=sin(35/18*pi);
   Ama(k,3*(k-1)+2)=cos(35/18*pi);
   Ama(k,3*(k-1)+3)=0;
   Yma(k,1)=MPALSAR(prc(k,1),prc(k,2));
   % Pma(k)=1/(PALSARMAIPRECISION(prc(k,1),prc(k,2))*PALSARMAIPRECISION(prc(k,1),prc(k,2)));
   Pma(k)=1/0.01168540938;
end
clear MPALSAR PALSARMAIPRECISION;
Ymd=zeros(size(prc,1),1);
Amd=zeros(size(prc,1),3*size(prc,1));
Pmd=zeros(size(prc,1),1);
for k=1:size(prc,1)                       
   Amd(k,3*(k-1)+1)=sin(19/18*pi);
   Amd(k,3*(k-1)+2)=cos(19/18*pi);
   Amd(k,3*(k-1)+3)=0;                         
   Ymd(k,1)=MENVISAT(arc(k,1),arc(k,2));
   % Pmd(k)=1/(ENVISATMAIPRECISION(arc(k,1),arc(k,2))*ENVISATMAIPRECISION(arc(k,1),arc(k,2)));
   Pmd(k)=1/0.71935237628;
end
clear MENVISAT ENVISATMAIPRECISION arc;

s=0;
T=0;
for k=1:size(prc,1)
   if ~isnan(Pa(k))
    s=s+Pa(k);
    T=T+1;
   end
end
for k=1:size(prc,1)
   if isnan(Pa(k))
    Pa(k)=s/T;
   end
end

s=0;
T=0;
for k=1:size(prc,1)
   if ~isnan(Pd(k))
    s=s+Pd(k);
    T=T+1;
   end
end
for k=1:size(prc,1)
   if isnan(Pd(k))
    Pd(k)=s/T;
   end
end

s=0;
T=0;
for k=1:size(prc,1)
   if ~isnan(Pma(k))
    s=s+Pma(k);
    T=T+1;
   end
end
for k=1:size(prc,1)
   if isnan(Pma(k))
    Pma(k)=s/T;
   end
end

s=0;
T=0;
for k=1:size(prc,1)
   if ~isnan(Pmd(k))
    s=s+Pmd(k);
    T=T+1;
   end
end
for k=1:size(prc,1)
   if isnan(Pmd(k))
    Pmd(k)=s/T;
   end
end

dir = strcat('P:\s\japan\TIFF\ALOStiff_401\GPS', Pline);
dir = strcat(dir,'e.tif');
[GPSE] = geotiffread(dir);
dir = strcat('P:\s\japan\TIFF\ALOStiff_401\pre', Pline);
dir = strcat(dir,'e.tif');
[PREGPSE] = geotiffread(dir);
dir = strcat('P:\s\japan\TIFF\ALOStiff_401\GPS', Pline);
dir = strcat(dir,'n.tif');
[GPSN] = geotiffread(dir);
dir = strcat('P:\s\japan\TIFF\ALOStiff_401\pre', Pline);
dir = strcat(dir,'n.tif');
[PREGPSN] = geotiffread(dir);
dir = strcat('P:\s\japan\TIFF\ALOStiff_401\GPS', Pline);
dir = strcat(dir,'u.tif');
[GPSU] = geotiffread(dir);
dir = strcat('P:\s\japan\TIFF\ALOStiff_401\pre', Pline);
dir = strcat(dir,'u.tif');
[PREGPSU] = geotiffread(dir);
Yg=zeros(3*size(prc,1),1);
Pg=zeros(3*size(prc,1),1);
for k=1:size(prc,1)               
   % Pg(3*(k-1)+1)=1/PREGPSE(prc(k,1),prc(k,2));
   % Pg(3*(k-1)+2)=1/PREGPSN(prc(k,1),prc(k,2));
   % Pg(3*(k-1)+3)=1/PREGPSU(prc(k,1),prc(k,2));
   Pg(3*(k-1)+1)=1/0.00210794697;
   Pg(3*(k-1)+2)=1/0.00210794697;
   Pg(3*(k-1)+3)=1/(3*0.00210794697);
   Yg(3*(k-1)+1)=GPSE(prc(k,1),prc(k,2));
   Yg(3*(k-1)+2)=GPSN(prc(k,1),prc(k,2));
   Yg(3*(k-1)+3)=GPSU(prc(k,1),prc(k,2));
end

clear GPSE PREGPSE GPSN PREGPSN GPSU PREGPSU;

AaTBa=zeros(size(prc,1)*3,6);
Ba=zeros(size(prc,1),6);
for  k=1:size(prc,1)
    Ba(k,1:6)=[1 pos(k,1) pos(k,2) pos(k,1)*pos(k,2) pos(k,1)*pos(k,1) pos(k,2)*pos(k,2)];
end
AaTBa=Aa'*diag(Pa)*Ba;
BaTBa=Ba'*diag(Pa)*Ba;

AdTBd=zeros(size(prc,1)*3,6);
Bd=zeros(size(prc,1),6);
for  k=1:size(prc,1)
    Bd(k,1:6)=[1 pos(k,1) pos(k,2) pos(k,1)*pos(k,2) pos(k,1)*pos(k,1) pos(k,2)*pos(k,2)];
end
AdTBd=Ad'*diag(Pd)*Bd;
BdTBd=Bd'*diag(Pd)*Bd;

AmaTBma=zeros(size(prc,1)*3,3);
Bma=zeros(size(prc,1),3);
for  k=1:size(prc,1)
    Bma(k,1:3)=[1 pos(k,1) pos(k,2)];
end
AmaTBma=Ama'*diag(Pma)*Bma;
BmaTBma=Bma'*diag(Pma)*Bma;

AmdTBmd=zeros(size(prc,1)*3,3);
Bmd=zeros(size(prc,1),3);
for  k=1:size(prc,1)
    Bmd(k,1:3)=[1 pos(k,1) pos(k,2)];
end
AmdTBmd=Amd'*diag(Pmd)*Bmd;
BmdTBmd=Bmd'*diag(Pmd)*Bmd;

B=[AaTBa AdTBd AmaTBma AmdTBmd];
clear AaTBa AdTBd AmaTBma AmdTBmd;

Q=diag(1./Pa)+Aa*diag(1./Pg)*Aa';
P=zeros(size(Pa,1),1);
for i=1:size(Pa,1)
    P(i)=1/Q(i,i);
end
Pthetaa=Ba'*diag(P)*Ba;
thetaa=pinv(Pthetaa)*Ba'*diag(P)*(Ya-Aa*Yg);

clear Q P;

Q=diag(1./Pd)+Ad*diag(1./Pg)*Ad';
P=zeros(size(Pd,1),1);
for i=1:size(Pd,1)
    P(i)=1/Q(i,i);
end
Pthetad=Bd'*diag(P)*Bd;
thetad=pinv(Pthetad)*Bd'*diag(P)*(Yd-Ad*Yg);
clear Q P;

Q=diag(1./Pma)+Ama*diag(1./Pg)*Ama';
P=zeros(size(Pma,1),1);
for i=1:size(Pma,1)
    P(i)=1/Q(i,i);
end
Pthetama=Bma'*diag(P)*Bma;
thetama=pinv(Pthetama)*Bma'*diag(P)*(Yma-Ama*Yg);
clear Q P;

Q=diag(1./Pmd)+Amd*diag(1./Pg)*Amd';
P=zeros(size(Pmd,1),1);
for i=1:size(Pmd,1)
    P(i)=1/Q(i,i);
end
Pthetamd=Bmd'*diag(P)*Bmd;
thetamd=pinv(Pthetamd)*Bmd'*diag(P)*(Ymd-Amd*Yg);
clear Q P;

BaTBa=BaTBa+Pthetaa;
BdTBd=BdTBd+Pthetad;
BmaTBma=BmaTBma+Pthetama;
BmdTBmd=BmdTBmd+Pthetamd;

datum=zeros(18,18);
datum(1:6,1:6)=BaTBa;
datum(7:12,7:12)=BdTBd;
datum(13:15,13:15)=BmaTBma;
datum(16:18,16:18)=BmdTBmd;

clear BaTBa BdTBd BmaTBma BmdTBmd;

N=zeros(size(prc,1)*3,size(prc,1)*3);
N=Aa'*diag(Pa)*Aa;
N=N+Ad'*diag(Pd)*Ad;
N=N+Ama'*diag(Pma)*Ama;
N=N+Amd'*diag(Pmd)*Amd;
N=N+diag(Pg); 

for k=1:size(prc,1)
    N(3*(k-1)+1:3*k,3*(k-1)+1:3*k)=inv(N(3*(k-1)+1:3*k,3*(k-1)+1:3*k));
end

Q22=pinv(datum-B'*N*B);

Q12=-N*B*Q22;

Q11=N-Q12*B'*N;

clear N;

for i=1:size(Q11,1)/3
std(i,1)=sqrt(abs(Q11(3*(i-1)+1,3*(i-1)+1)));
end
for i=1:size(Q11,1)/3
std(i,2)=sqrt(abs(Q11(3*(i-1)+2,3*(i-1)+2)));
end
for i=1:size(Q11,1)/3
std(i,3)=sqrt(abs(Q11(3*(i-1)+3,3*(i-1)+3)));
end

L=zeros(size(prc,1)*3+18,1);
L(1:size(prc,1)*3,1)=Aa'*diag(Pa)*Ya;
L(1:size(prc,1)*3,1)=L(1:size(prc,1)*3,1)+Ad'*diag(Pd)*Yd;
L(1:size(prc,1)*3,1)=L(1:size(prc,1)*3,1)+Ama'*diag(Pma)*Yma;
L(1:size(prc,1)*3,1)=L(1:size(prc,1)*3,1)+Amd'*diag(Pmd)*Ymd;
L(1:size(prc,1)*3,1)=L(1:size(prc,1)*3,1)+diag(Pg)*Yg;

L(size(prc,1)*3+1:size(prc,1)*3+6,1)=Ba'*diag(Pa)*Ya+Pthetaa*thetaa;
L(size(prc,1)*3+7:size(prc,1)*3+12,1)=Bd'*diag(Pd)*Yd+Pthetad*thetad;

L(size(prc,1)*3+13:size(prc,1)*3+15,1)=Bma'*diag(Pma)*Yma+Pthetama*thetama;
L(size(prc,1)*3+16:size(prc,1)*3+18,1)=Bmd'*diag(Pmd)*Ymd+Pthetamd*thetamd;

% clear Aa Pa Ya Ad Pd Yd Ama Pma Yma Amd Pmd Ymd Pg Yg;

X1=Q11*L(1:size(prc,1)*3,1)+Q12*L(1+size(prc,1)*3:18+size(prc,1)*3,1);

X2=Q12'*L(1:size(prc,1)*3,1)+Q22*L(1+size(prc,1)*3:18+size(prc,1)*3,1);

 X=[X1;X2];
 Xd=[X(1:size(prc,1)*3,1);X(size(prc,1)*3+7:size(prc,1)*3+12,1)];
 Xma=[X(1:size(prc,1)*3,1);X(size(prc,1)*3+13:size(prc,1)*3+15,1)];
 Xmd=[X(1:size(prc,1)*3,1);X(size(prc,1)*3+16:size(prc,1)*3+18,1)];
 e=[Ya-[Aa Ba]*X(1:size(prc,1)*3+6,1);Yd-[Ad Bd]*Xd;Yma-[Ama Bma]*Xma;Ymd-[Amd Bmd]*Xmd;Yg-X(1:size(prc,1)*3, 1)];
% 
% h=zeros(5,1);
% q=zeros(5,1);
% N0=zeros(18, 18);
% N0(1:6, 1:6)=Pthetaa;
% N0(7:12, 7:12)=Pthetad;
% N0(13:15, 13:15)=Pthetama;
% N0(16:18, 16:18)=Pthetamd;
% 
% for i=1:5
% N_i=zeros(size(prc,1)*3+18, size(prc,1)*3+18);
% if i==1
% N_i(1:size(prc,1)*3, 1:size(prc,1)*3)=Aa'*diag(Pa)*Aa;
% N_i(1:size(prc,1)*3, size(prc,1)*3+1:size(prc,1)*3+6)=Aa'*diag(Pa)*Ba;
% N_i(size(prc,1)*3+1:size(prc,1)*3+6, 1:size(prc,1)*3)=N_i(1:size(prc,1)*3, size(prc,1)*3+1:size(prc,1)*3+6)';
% N_i(size(prc,1)*3+1:size(prc,1)*3+6, size(prc,1)*3+1:size(prc,1)*3+6)=Ba'*diag(Pa)*Ba;
% M=diag(Pa);
% elseif i==2
% N_i(1:size(prc,1)*3, 1:size(prc,1)*3)=Ad'*diag(Pd)*Ad;
% N_i(1:size(prc,1)*3, size(prc,1)*3+7:size(prc,1)*3+12)=Ad'*diag(Pd)*Bd;
% N_i(size(prc,1)*3+7:size(prc,1)*3+12, 1:size(prc,1)*3)=N_i(1:size(prc,1)*3, size(prc,1)*3+7:size(prc,1)*3+12)';
% N_i(size(prc,1)*3+7:size(prc,1)*3+12, size(prc,1)*3+7:size(prc,1)*3+12)=Bd'*diag(Pd)*Bd;
% M=diag(Pd);
% elseif i==3
% N_i(1:size(prc,1)*3, 1:size(prc,1)*3)=Ama'*diag(Pma)*Ama;
% N_i(1:size(prc,1)*3, size(prc,1)*3+13:size(prc,1)*3+15)=Ama'*diag(Pma)*Bma;
% N_i(size(prc,1)*3+13:size(prc,1)*3+15, 1:size(prc,1)*3)=N_i(1:size(prc,1)*3, size(prc,1)*3+13:size(prc,1)*3+15)';
% N_i(size(prc,1)*3+13:size(prc,1)*3+15, size(prc,1)*3+13:size(prc,1)*3+15)=Bma'*diag(Pma)*Bma;
% M=diag(Pma);
% elseif i==4
% N_i(1:size(prc,1)*3, 1:size(prc,1)*3)=Amd'*diag(Pmd)*Amd;
% N_i(1:size(prc,1)*3, size(prc,1)*3+16:size(prc,1)*3+18)=Amd'*diag(Pmd)*Bmd;
% N_i(size(prc,1)*3+16:size(prc,1)*3+18, 1:size(prc,1)*3)=N_i(1:size(prc,1)*3, size(prc,1)*3+16:size(prc,1)*3+18)';
% N_i(size(prc,1)*3+16:size(prc,1)*3+18, size(prc,1)*3+16:size(prc,1)*3+18)=Bmd'*diag(Pmd)*Bmd;
% M=diag(Pmd);
% elseif i==5
% N_i(1:size(prc,1)*3, 1:size(prc,1)*3)=diag(Pg);
% M=diag(Pg);
% end    
% 
% h(i,1) = size(prc,1) - trace(Q11*N_i(1:size(prc,1)*3, 1:size(prc,1)*3));
% h(i,1) = h(i,1) - trace(Q12*N_i(size(prc,1)*3+1:size(prc,1)*3+18, 1:size(prc,1)*3)); 
% h(i,1) = h(i,1) - trace(Q12'*N_i(1:size(prc,1)*3, size(prc,1)*3+1:size(prc,1)*3+18) + Q22*N_i(size(prc,1)*3+1:size(prc,1)*3+18, size(prc,1)*3+1:size(prc,1)*3+18));

%h0 = trace((Q12'*N_i(1:size(prc,1)*3, 1:size(prc,1)*3)+Q22*N_i(size(prc,1)*3+1:size(prc,1)*3+18, 1:size(prc,1)*3))*Q12*N0);
%h0 = h0 + trace((Q12'*N_i(1:size(prc,1)*3, size(prc,1)*3+1:size(prc,1)*3+18) + Q22*N_i(size(prc,1)*3+1:size(prc,1)*3+18, size(prc,1)*3+1:size(prc,1)*3+18))*Q22*N0);
% h0 = 0;
% if i~=5 
%     q(i,1) = e((i-1)*size(prc,1)+1:i*size(prc,1), 1)'*M*e((i-1)*size(prc,1)+1:i*size(prc,1), 1)-h0;
% elseif i==5
%     q(i,1) = e(4*size(prc,1)+1:7*size(prc,1), 1)'*M*e(4*size(prc,1)+1:7*size(prc,1), 1)-h0;
% end
% 
% if i~=5 
%     h(i,1) = 4*size(prc,1)/7;
% elseif i==5
%     h(i,1) = 12*size(prc,1)/7;
% end
% 
% disp(q(i,1)/h(i,1));
% 
% end



% clear Q11 X1 X2 L ;

% for i=1:size(pos,1)
%     E(i,1)=X(3*(i-1)+1);N(i,1)=X(3*(i-1)+2);U(i,1)=X(3*(i-1)+3);
% end
% lon=pos(:,2); pos(:,2)=pos(:,1); pos(:,1)=lon; 
% ENU=[E N U];
% clear lon Aa Ba Ama Bma Ya Yd Ymd Yma Pma  Amd Pmd Pg Yg E N U;
% 
% 
% dir = strcat('P:\s\japan\TIFF\ALOStiff_401\DA', Pline);
% dir = strcat(dir,'.tif');
% [PALSAR,PALSARMAP] = geotiffread(dir);
% A = NaN(size(PALSAR,1), size(PALSAR,2));
% clear PALSAR;
% 
% for i=1:size(prc,1)
%     A(prc(i,1),prc(i,2))=X1(3*(i-1)+1,1);
% end
% dir = strcat('P:\s\japan\TIFF\ALOStiff_401\', mosaic_mask);
% dir = strcat(dir, Pline);
% dir = strcat(dir,'e.tif');
% geotiffwrite(dir,A, PALSARMAP);
% 
% for i=1:size(prc,1)
%     A(prc(i,1),prc(i,2))=X1(3*(i-1)+2,1);
% end
% dir = strcat('P:\s\japan\TIFF\ALOStiff_401\', mosaic_mask);
% dir = strcat(dir, Pline);
% dir = strcat(dir,'n.tif');
% geotiffwrite(dir,A, PALSARMAP);
% 
% for i=1:size(prc,1)
%     A(prc(i,1),prc(i,2))=X1(3*(i-1)+3,1);
% end
% dir = strcat('P:\s\japan\TIFF\ALOStiff_401\', mosaic_mask);
% dir = strcat(dir, Pline);
% dir = strcat(dir,'u.tif');
% geotiffwrite(dir,A, PALSARMAP);
% 
% dir = strcat('P:\s\japan\TIFF\ALOStiff_401\', mosaic_mask);
% dir = strcat(dir, Pline);
% dir = strcat(dir, 'prc.txt');
% fileID = fopen(dir, 'w');
% for i=1:size(prc,1)
%     fprintf(fileID, '%d %d\n', prc(i,1), prc(i,2));
% end
% fclose(fileID);
% 
% dir = strcat('P:\s\japan\TIFF\ALOStiff_401\', mosaic_mask);
% dir = strcat(dir, Pline);
% dir = strcat(dir, 'datum.txt');
% fileID = fopen(dir, 'w');
% for i=1:size(18,1)
%     fprintf(fileID, '%f\n', X2(i,1));
% end
% fclose(fileID);
% 
% dir = strcat('P:\s\japan\TIFF\ALOStiff_401\', mosaic_mask);
% dir = strcat(dir, Pline);
% dir = strcat(dir, 'Q.txt');
% fileID = fopen(dir, 'w');
% for k=1:size(prc,1)
%     fprintf(fileID, '%f %f %f\n', Q11(3*(k-1)+1,3*(k-1)+1), Q11(3*(k-1)+1,3*(k-1)+2), Q11(3*(k-1)+1,3*(k-1)+3));
%     fprintf(fileID, '%f %f %f\n', Q11(3*(k-1)+2,3*(k-1)+1), Q11(3*(k-1)+2,3*(k-1)+2), Q11(3*(k-1)+2,3*(k-1)+3));
%     fprintf(fileID, '%f %f %f\n', Q11(3*(k-1)+3,3*(k-1)+1), Q11(3*(k-1)+3,3*(k-1)+2), Q11(3*(k-1)+3,3*(k-1)+3));
% end
% fclose(fileID);
% 
% dir = strcat('P:\s\japan\TIFF\ALOStiff_401\', mosaic_mask);
% dir = strcat(dir, Pline);
% dir = strcat(dir, 'Qdatum.txt');
% fileID = fopen(dir, 'w');
% for k=1:18
%     fprintf(fileID, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n', Q22(k,1), Q22(k,2), Q22(k,3), Q22(k,4), Q22(k,5), Q22(k,6), ...
%         Q22(k,7), Q22(k,8), Q22(k,9), Q22(k,10), Q22(k,11), Q22(k,12), Q22(k,13), Q22(k,14), Q22(k,15), Q22(k,16), Q22(k,17), Q22(k,18));
% end
% fclose(fileID);

end