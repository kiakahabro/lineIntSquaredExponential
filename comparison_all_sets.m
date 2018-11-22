
n = 6;
ntrials = 10000;




uppertri = logical(triu(ones(n,n)));

nparamsets = 8;

matI2time = NaN(nparamsets,ntrials);
simps10time = matI2time;
simps200time = matI2time;
newMethtime = matI2time;
simps10err = matI2time;
simps200err = matI2time;
newMetherr = matI2time;
BNtime = matI2time;
BNerr = matI2time;


for p = 1:nparamsets
    
for i = 1:ntrials
if p == 1           % set 1             % standard set
    u = rand(n,1);
    w = rand(n,1);
    x = rand(n,1);
    V = eye(n);
    sV = chol(V);
elseif p==2         % set 2             % almsot colinear
    u = rand(n,1);
    w = rand(n,1);
    x = w+1e-8*rand(n,1);
    V = eye(n);
    sV = chol(V);
elseif p ==3        % set 3         % diagonal randon V 
    u = rand(n,1);
    w = rand(n,1);
    x = rand(n,1);
    V = diag(rand(n,1));
    sV = chol(V);
elseif p==4             % set 4        % w = zeros()
    u = rand(n,1);
    w = zeros(n,1);
    x = rand(n,1);
    V = eye(n);
    sV = chol(V);
elseif p ==5        % set 5         % scaled nicely
    u = 10*rand(n,1);
    w = 10*rand(n,1);
    x = 10*rand(n,1);
    V = diag(0.01*rand(n,1));
    sV = chol(V);
    
elseif p==6         % set 6     % scaled poorly
    u = 10*rand(n,1);
    w = 10*rand(n,1);
    x = 10*rand(n,1);
    V = diag(10*rand(n,1));
    sV = chol(V);
    
elseif p==7         % set 7     % first line close to zero
    u = rand(n,1);
    w = 1e-8*rand(n,1);
    x = rand(n,1);
    V = eye(n);
    sV = chol(V);
elseif p==8    % both lines close to zero
    u = rand(n,1);
    w = 1e-8*rand(n,1);
    x = 1e-8*rand(n,1);
    V = eye(n);
    sV = chol(V);
else
   error('invalid param settings') 
    
end



%Using Matlab's integral2
z3   = u;
d1   = z3'*V*z3;
d2   = z3'*V*w;
d3   = z3'*V*x;
d4   = w'*V*x;
d5   = w'*V*w;
d6   = x'*V*x;
L1 = norm(w);
L2 = norm(x);

tic
I1 = L1*L2*integral2(@(t,v) exp(-0.5*(d1+2*t*d2-2*v*d3-2*t.*v.*d4+t.*t*d5+v.*v*d6)),...
    0,1,0,1,'AbsTol',eps,'RelTol',eps);
matI2time(p,i) = toc;

if isnan(I1)
    keyboard
end

% using 2D simpsosn rule

order = 10;  %Higher order => higher accuracy
tic
I2 = intTwoSimps(u,w,x,V,order,order);
simps10time(p,i) = toc;
simps10err(p,i) = abs(I1 - I2);

order = 200;  
tic
I2 = intTwoSimps(u,w,x,V,order,order);
simps200time(p,i) = toc;
simps200err(p,i) = abs(I1 - I2);

% Bivariate normal method
tic
i4 = L1*L2*intTwo_BN(u,w,x,V);
BNtime(p,i) = toc;
BNerr(p,i) = abs(I1-i4);

if isnan(i4)
    keyboard
end


% using new method
tic
OUT = intTwoK(u,w,x,V);            %% OUT presently contains [result,abserr,nevals,info]
newMethtime(p,i) = toc;
newMetherr(p,i) = abs(I1 - OUT(1));

end

p
end


%% plot resutls
figure(1)
clf
subplot 211
plot(log10(mean(matI2time,2)))
hold on
plot(log10(mean(simps10time,2)),'-x')
plot(log10(mean(simps200time,2)),'-o')
plot(log10(mean(BNtime,2)),'--')
plot(log10(mean(newMethtime,2)),'-s')
hold off
title('1-norm run times')
ylabel('log10 of time')
xlabel('param set')
legend('matlab integral2','simpsons order 10','simpsons order 200','BN','new method')
set(gca,'XTick',[1:nparamsets])

subplot 212
plot(log10(mean(simps10err,2)),'-x')
hold on
plot(log10(mean(simps200err,2)),'-o')
plot(log10(mean(BNerr,2)),'-x')
plot(log10(mean(newMetherr,2)),'-s')
hold off
title('1-norm error')
ylabel('log10 of error')
xlabel('param set')
legend('simpsons order 10','simpsons order 200','BN','new method')
set(gca,'XTick',[1:nparamsets])


