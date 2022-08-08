function Noise=fGn_Davies(H,sigma,L)
%�ú��������� CEM ����ģ���������Ϊ N �ķ�����˹����
%H----- ����Ƶ��ָ��
%sigma -----�����ľ�����
%L-------���ݳ���

N=L-1;
k=0:N;
% K=sqrt(2*N-2);%K Ϊ����Ҷ�任�ĳ�������

rGH=sigma^2*(abs(k+1).^(2*H)-2*abs(k).^(2*H)+abs(k-1).^(2*H))/2;
s=[rGH rGH(N:-1:2)];

% k=0:2*N-1;
% n=k;
% WN=exp(-j*2*pi/(2*N));
% nk=n'*k;
% WNnk=WN.^nk;
% S1=s*WNnk;
S=fft(s);
% S=real(S);
% S=(2*N)*ifft(s);


tempX=[normrnd(0,sqrt(2),1,1) normrnd(0,1,1,N-1) normrnd(0,sqrt(2),1,1)];
tempY=[0 normrnd(0,1,1,N-1) 0];
tempZ=tempX+j*tempY;
Z=[tempZ conj(tempZ(N:-1:2))];

X=N^(1/2)*ifft(Z.*(S.^(1/2)));
Noise=X(1:L);