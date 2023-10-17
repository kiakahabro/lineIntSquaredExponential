clc
clear

addpath build/

n = 6;
ntrials = 10000;
T1 = NaN(ntrials,1);
T2 = T1;
T3 = T1;
T4 = T1;
T5 = T1;
error2 = T1;
error3 = T1;
error4 = T1;
error5 = T1;


% rng(1)

for i = 1:ntrials
u = rand(n,1);
w = rand(n,1);
x = w+1e-8*randn(n,1);
V = eye(n);
L1 = norm(w);
L2 = norm(x);


%Using Matlab's integral2
z3   = u;
d1   = z3'*V*z3;
d2   = z3'*V*w;
d3   = z3'*V*x;
d4   = w'*V*x;
d5   = w'*V*w;
d6   = x'*V*x;

tic
I1 = L1*L2*integral2(@(t,v) exp(-0.5*(d1+2*t*d2-2*v*d3-2*t.*v.*d4+t.*t*d5+v.*v*d6)),...
    0,1,0,1,'AbsTol',eps,'RelTol',eps);
T1(i) = toc;


% using 2D simpsosn rules

order = 10;  %Higher order => higher accuracy
tic
I2 = intTwoSimps(u,w,x,V,order,order);
T2(i) = toc;
error2(i) = abs(I1 - I2);

order = 200;  %Higher order => higher accuracy
tic
I2 = intTwoSimps(u,w,x,V,order,order);
T3(i) = toc;
error3(i) = abs(I1 - I2);


% bivariate normal method

tic
i4 = L1*L2*intTwo_BN(u,w,x,V);
T4(i) = toc;
error4(i) = abs(I1-i4);

% proposed method
tic
OUT2 = intTwoK(u,w,x,V); % OUT presently contains [result,abserr,nevals,info]

T5(i) = toc;
error5(i) = abs(I1-OUT2(1));

i
end


%%
% mean(I1)
disp(['Integral2: t = ' num2str(mean(T1)) 's'])
disp(['2D simpsons order 10: t = ' num2str(mean(T2)) 's, error mean = ' num2str(mean(error2))])
disp(['2D simpsons order 200: t = ' num2str(mean(T3)) 's, error mean = ' num2str(mean(error3))])
disp(['bivariate normal: t = ' num2str(mean(T4)) 's, error mean = ' num2str(mean(error4))])
disp(['alternate: t = ' num2str(mean(T5)) 's, error mean = ' num2str(mean(error5))])

%%
fig = figure(1);
clf
subplot 221
hist(error2,100)
title('2D Simpsons order 10')
axis square
set(gca,'fontsize',14)

subplot 222
hist(error3,100)
title('2D Simpsons rule order 200')
axis square
set(gca,'fontsize',14)

subplot 223
bins = linspace(0,10,100);
hist(error4,bins)
title('Bivariate Normal Method')
axis square
set(gca,'fontsize',14)
xlim([0,10])

subplot 224
bins = linspace(0,5e-14,100);
hist(error5,bins)
title('Proposed method')
set(gca,'fontsize',14)
axis square
xlim([0,max(bins)])

