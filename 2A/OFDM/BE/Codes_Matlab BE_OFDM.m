clear all
close all 
N=16;
Fe=N;
bits=randi([0,1],1,5000*N);
sym=2*bits-1;
M=reshape(sym,16,5000);




%%%% Chaine de transmission OFDM sans canal
% emission 1
M1=reshape(sym,16,5000);
L1=M1(1,:);
M1=[L1;zeros(N-1,5000)];
mat_sortie1=ifft(M1);
mat_emis1=reshape(mat_sortie1,1,N*5000);
[DSP1 f1]=pwelch(mat_emis1,[],[],[],Fe,'twosided');
Freq=linspace(0,Fe,length(DSP1));
figure
plot(Freq,10*log(DSP1));
title("Densité spectrale du signal emis: une seule porteuse activé")
xlabel("Frequences")
ylabel("DSP")


% emission 2
L2=M(1,:);
L3=M(2,:);
M2=[L2;L3;zeros(N-2,5000)];

mat_sortie2=ifft(M2);
mat_emis2=reshape(mat_sortie2,1,N*5000);
[DSP2 f2]=pwelch(mat_emis2,[],[],[],Fe,'twosided');
Freq2=linspace(0,Fe,length(DSP2));
figure
plot(Freq2,10*log(DSP2)')
title("Densité spectrale du signal emis: deux porteuses activé")
xlabel("Frequences")
ylabel("DSP")

% emission 3
 L=zeros(N,5000);
for i=1:8
    L(4+i,:)=M(4+i,:);
end
LL=[];
for i=1:8
    LL=[LL; M(4+i,:)];
end

M3=L;
mat_sortie3=ifft(M3);
mat3=reshape(mat_sortie3,1,N*5000);
[DSP3 f]=pwelch(mat3,[],[],[],Fe,'twosided');
Freq3=linspace(0,Fe,length(DSP3));
figure
plot(Freq3,10*log(DSP3)')

title("Densité spectrale du signal emis: 8 porteuses centrales activé")
xlabel("Frequences")
ylabel("DSP")

% reception 1
mat_recu1=mat_emis1;
mat_entre1=reshape(mat_recu1,N,5000);
M_recu1=fft(mat_entre1);
a1=real(M_recu1);
a1=sign(a1);
a1=a1(1,:)==M1(1,:);
a1=sum(sum(a1));
TEB1=1-a1/(5000*1)


% reception 2
mat_recu2=fft(mat_sortie2);
a2=real(mat_recu2);
a2=[a2(1,:);a2(2,:)];
a2=sign(a2)==[M2(1,:);M2(2,:)];
a2=sum(sum(a2));
TEB2=1-a2/(5000*2)


% reception 3
mat_recu3=fft(mat_sortie3);
a3=real(mat_recu3);
L=[];
for i=1:8
    L=[L ;a3(4+i,:)];
end

a3=sign(L)==LL;
a3=sum(sum(a3));

TEB3=1-a3/(5000*8)





%%%% Chaine OFDM avec canal Muti-trajets
Ts=1;
Te=Ts;
Fe=1/Te;
N=16;
%%% Implémentation sans intervalle de garde 

%le filtre du canal h
a0=0.227;
a1=0.46;
a2=0.688;
a3=0.46;
a4=0.227;

figure
subplot(2,1,1)
h=[a0 a1 a2 a3 a4];
h=[h zeros(1,11)];
H=fft(h);
plot(linspace(0,Fe,length(H)),abs(H))
title("Module de la reponse frequence du canal")
subplot(2,1,2)
plot(linspace(0,Fe,length(H)),angle(H))
title("Phase de la reponse frequenciel du canal ")

%%%emission avec 16 porteur sans intervalle de garde
%emission 
M=reshape(sym,16,5000);
mat_sortie=ifft(M);
mat_emis=reshape(mat_sortie,1,N*5000);

[DSP f]=pwelch(mat_emis,[],[],[],Fe,'twosided');
Freq=linspace(0,Fe,length(DSP));
figure
plot(Freq,10*log(DSP)');
title("Densité spectrale du signal emis avant le canal")
xlabel("Frequences")
ylabel("DSP")

%canal
signal_recu=filter(h,1,mat_emis);
[DSP f]=pwelch(signal_recu,[],[],[],Fe,'twosided');
Freq=linspace(0,Fe,length(DSP));
figure
plot(Freq,10*log(DSP)');
title("Densité spectrale du signal recu après le canal")
xlabel("Frequences")
ylabel("DSP")

%reception 
mat_recu=signal_recu;
mat_entre=reshape(mat_recu,N,5000);
M_recu=fft(mat_entre);
R=real(M_recu);
sym_decide=sign(R);
symboles1=M_recu(1,:);
symboles16=M_recu(16,:);
figure

plot(real(symboles16),imag(symboles16),'*')
hold on 
plot(real(symboles1),imag(symboles1),'*')
title("Les constellations sur la porteuse 1 et 16 sans intervalle de garde")
A=sym_decide==M;
s=sum(sum(A));
teb=1-s/(16*5000)


%%%emission avec intervalle de garde composé des zeros
%emission
mat_zeros=zeros(4,5000);
M=reshape(sym,16,5000);
mat_sortie=ifft(M);
mat_IG=[mat_zeros ; mat_sortie];
mat_emis=reshape(mat_IG,1,(4+N)*5000);

%canal
signal_recu=filter(h,1,mat_emis);


%reception 
mat_recu=signal_recu;
mat_entre=reshape(mat_recu,N+4,5000);
mat_entre=mat_entre(5:20,:);
M_recu=fft(mat_entre);
a=real(M_recu);
aa=sign(a);

symboles1=M_recu(1,:);
symboles16=M_recu(16,:);
figure
plot(real(symboles16),imag(symboles16),'*')
hold on 
plot(real(symboles1),imag(symboles1),'*')

title("Constellations sur les porteuses 1 et 16 avec intervalle de garde")
a=aa==M;
a=sum(sum(a));
teb=1-a/(16*5000)




%%%emission avec prefix cyclique sans égalisation
%emission
M=reshape(sym,16,5000);
mat_sortie=ifft(M);
mat_pre_cyc=[mat_sortie(13:16,:); mat_sortie];
mat_emis=reshape(mat_pre_cyc,1,(4+N)*5000);

%canal
signal_recu=filter(h,1,mat_emis);

%reception 
mat_recu=signal_recu;
mat_entre=reshape(mat_recu,N+4,5000);
mat_entre=mat_entre(5:20,:);
M_recu=fft(mat_entre);
a=real(M_recu);
aa=sign(a);
symboles1=M_recu(1,:);
symboles16=M_recu(16,:);

figure
plot(real(symboles16),imag(symboles16),'*')
hold on 
plot(real(symboles1),imag(symboles1),'*')

title("constellation sur les porteuses 1 et 16 avec prefixe cyclique sans égalisation")
a=aa==M;
a=sum(sum(a));
teb=1-a/(16*5000)

 
%%% préfixe cyclique avec égalisation  

%emission
M=reshape(sym,16,5000);
mat_sortie=ifft(M);
mat_pre_cyc=[mat_sortie(13:16,:); mat_sortie];
mat_emis=reshape(mat_pre_cyc,1,(4+N)*5000);

%canal
signal_recu=filter(h,1,mat_emis);

%reception 
mat_recu=signal_recu;
mat_entre=reshape(mat_recu,N+4,5000);
mat_entre=mat_entre(5:20,:);
M_recu=fft(mat_entre);

b_zf=M_recu./H.'; % (methode  ZF: zero forcing )
b_ml=M_recu.*H'; %methode ML: maximum de vraisemblance

a=real(b_zf);
aa=sign(a);

symboles1_zf=b_zf(1,:);
symboles2_zf=b_zf(16,:);

symboles1_ml=b_ml(1,:);
symboles2_ml=b_ml(16,:);





a=aa==M;
a=sum(sum(a));
teb=1-a/(16*5000)

figure
subplot (2,1,1)
plot(real(symboles1_zf),imag(symboles1_zf),'*')
hold on 
plot(real(symboles2_zf),imag(symboles2_zf),'*')
axis([-5  5 -5 5])
title("Constellation avec prefixe cyclique et égalisation ZF")

subplot (2,1,2)
plot(real(symboles1_ml),imag(symboles1_ml),'*')
hold on 
plot(real(symboles2_ml),imag(symboles2_ml),'*')
axis([-5  5 -5 5])
title("Constellation avec prefixe cyclique et égalisation ML")



%%%%  Erreur d'horloge
h=[a0 a1 a2 a3 a4];
h=[h zeros(1,11)];
H=fft(h);
N=16;


%%% CAS 2: 

%emission
M=reshape(sym,16,5000);
mat_sortie1=ifft(M);
mat_pre_cyc=[mat_sortie1(9:16,:); mat_sortie1];
mat_emis=reshape(mat_pre_cyc,1,(8+N)*5000);

% canal
signal_recu1=filter(h,1,mat_emis);

% reception 
mat_recu1=signal_recu1;
mat_entre1=reshape(mat_recu1,N+8,5000);
mat_entre1=mat_entre1(6:21,:);
M_recu1=fft(mat_entre1);

b=M_recu./H.';  %(methode  ZF: zero forcing )
%b=M_recu1.*H'; %methode ML: maximum de vraisemblance
a=real(b);
aa=sign(a);
symboles=[M_recu1(10,:)];


a=aa==M;
a=sum(sum(a));
teb=1-a/(16*5000)
%scatterplot(symboles)
figure
symboles1=M_recu(1,:);
symboles16=M_recu(16,:);


plot(real(symboles16),imag(symboles16),'*')
hold on 
plot(real(symboles1),imag(symboles1),'*')
title("Constellations pour le cas2")


%%% CAS 1

% emission
M=reshape(sym,16,5000);
mat_sortie=ifft(M);
mat_pre_cyc=[mat_sortie(9:16,:); mat_sortie];
mat_emis=reshape(mat_pre_cyc,1,(8+N)*5000);

% canal
signal_recu=filter(h,1,mat_emis);

%reception 
mat_recu=signal_recu;
mat_entre=reshape(mat_recu,N+8,5000);
mat_entre=mat_entre(2:17,:);
M_recu=fft(mat_entre);

b=M_recu./H.'; % (methode  ZF: zero forcing )
%b=M_recu.*H'; % (methode ML: maximum de vraisemblance)
a=real(b);
aa=sign(a);

a=aa==M;
a=sum(sum(a));
teb=1-a/(16*5000)



symboles1=M_recu(1,:);
symboles16=M_recu(16,:);

figure
plot(real(symboles16),imag(symboles16),'*')
hold on 
plot(real(symboles1),imag(symboles1),'*')
title("Constellation pour le cas 1")


%%% CAS 3

% emission
M=reshape(sym,16,5000);
mat_sortie=ifft(M);
mat_pre_cyc=[mat_sortie(9:16,:); mat_sortie];
mat_emis=reshape(mat_pre_cyc,1,(8+N)*5000);

% canal
signal_recu=filter(h,1,mat_emis);

% reception 
mat_recu=signal_recu;
signal_recu=[signal_recu(3:end) signal_recu(1) signal_recu(2)];
mat_entre=reshape(signal_recu,N+8,5000);
mat_entre=mat_entre(9:24,:);
M_recu=fft(mat_entre);

b=M_recu./H.'; % (methode  ZF: zero forcing )
%b=M_recu.*H'; % (methode ML: maximum de vraisemblance)



a=real(b);
aa=sign(a);
symboles=b(10,:);



a=aa==M;
a=sum(sum(a));
teb=1-a/(16*5000)


symboles1=M_recu(1,:);
symboles16=M_recu(16,:);

figure
plot(real(symboles16),imag(symboles16),'*')
hold on 
plot(real(symboles1),imag(symboles1),'*')

title("Constellation pour le cas 3")





%%%estimation du canal dans le CAS 2


Sym_pilote=M(:,1);
Sym_recu=M_recu1(:,1);
H_estime=Sym_recu./Sym_pilote;

%b=M_recu1.*conj(H_estime); %methode ML: maximum de vraisemblance
b=M_recu1./H_estime; %methode ZF: Zero Forcing


a=real(b);
aa=sign(a);
symboles=[b(1,:)];


figure

plot(real(symboles),imag(symboles),'*')
axis([-1 1 -5 5])

%scatterplot(symboles)
title("Consellation pour le cas 2 après l'estimation du canal")
