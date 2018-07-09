close all
clear all

R0=20;       %mean radius. unit: Angstroms
PDI=0.1;       %unit: relative polydispersity
volF=0.1;
contrast=6e-6;    %contrast,angstrom^-1, scattering length density difference

Filename='testPolydisperseGaussianSphere.dat';

Q=[0.001:0.001:0.8];

N=10;

deltaR=R0*PDI;
SplitN=100;
R=[R0/N:deltaR/SplitN:R0+5*deltaR];

FGaussian=1/sqrt(2*pi*deltaR^2)*exp(-(R-R0).^2/2/deltaR^2);

Ni=FGaussian*deltaR/SplitN;
%trapz(R,FGaussian)
%sum(FGaussian*deltaR/SplitN)
sum(Ni)


Vi=4*pi*R.^3/3;
AverageVolume=sum(Ni.*Vi);


IQN=zeros(size(Q));
%IQ=sphereFormFactor(Q,R0);

for i=1:length(R)
    IQN=IQN+Ni(i)*Vi(i)^2*contrast^2*sphereFormFactor(Q,R(i));
    %IQ=1*AverageVolume^2*contrast^2*sphereFormFactor(Q,R0);
    %IQ=sphereFormFactor(Q,R(500));
end

IQ=volF/AverageVolume*IQN;

IQ=IQ*1e8;  %convert from A^-1 to cm^-1
loglog(Q,IQ,'b.-')

FQ=zeros(size(Q));

for i=1:length(R)
    FQ=FQ+Ni(i)*Vi(i)*contrast*3*(sin(Q*R(i)) - Q*R(i).*cos(Q*R(i))) ./ (Q*R(i)).^3;
end

beta=FQ.^2./IQN;
figure
semilogx(Q,beta,'g*-');

%stop

NormIQN=IQN/IQN(1);

FQ2=beta.*IQ;

normQ=Q*R0*2;
Sq=HSS_SQ(volF,normQ);

figure
loglog(Q,IQ,'bo-');
figure
plot(Q,beta,'g*-');
figure
plot(Q,Sq,'b.-');
figure
Sq_eff=1+beta.*(Sq-1);
plot(Q,Sq_eff,'gd-');

fileID=fopen(Filename,'w');

FileString=strcat('Filename:',Filename, '**R0=', num2str(R0),'**PDI=', num2str(PDI),...
    '**contrast',num2str(contrast), '**volF=',num2str(volF));
Outputformat='File format : Q, <FQ>, <FQ^2>, PQ, betaQ, SQ, SQ_Eff';
fprintf(fileID,'%s\n', FileString);
fprintf(fileID,'%s\n', Outputformat);
for i=1:length(Q)
    %file format : Q, <FQ>, <FQ^2>, PQ, betaQ, SQ, SQ_Eff
    fprintf(fileID,'%e\t%e\t%e\t%e\t%e\t%e\t%e\n',...
        Q(i),FQ(i),IQN(i),IQ(i),beta(i),Sq(i),Sq_eff(i));
end

fclose(fileID)

