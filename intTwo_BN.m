function i3 = intTwo_BN(u,w,x,V)

% Double Integral of exp(-0.5*(u+t*w-s*x).?*V*(u+t*w-s*x))
% between 0 and 1 on both integrals

f   = u.'*V*u;
d   = -u.'*V*w;
e   = u.'*V*x;
b   = -w.'*V*x;
a   = w.'*V*w;
c   = x.'*V*x;

% fcn = @(t,v) exp(-0.5*(a*t.*t -b*t.*v-b*v.*t+c*v.*v + 2*d*t - 2*e*v + f));


% %Backup code
% z3   = u;
% d1   = z3'*V*z3;
% d2   = z3'*V*w;
% d3   = z3'*V*x;
% d4   = w'*V*x;
% d5   = w'*V*w;
% d6   = x'*V*x;

% i2 = integral2(@(t,v) exp(-0.5*(d1+2*t*d2-2*v*d3-2*t.*v.*d4+t.*t*d5+v.*v*d6)),0,1,0,1,'AbsTol',eps,'RelTol',eps);
% %
%i4 = integral2(@(t,v) exp(-0.5*(a*t.*t + b*t.*v + b*v.*t + c*v.*v - 2*d*t - 2*e*v + f)),0,1,0,1,'AbsTol',eps,'RelTol',eps);
% %
% % i3 = simp2D(fcn,0,1,0,1,1,1);
% %
% % if abs(i2-i3) > 1e-4
% %     keyboard
% % end

% r=triu(qr(sqrtm(V)*[w -x]));
% r=r(1:2,1:2);


% integral2(@(t,v) exp(-0.5*((r(1,1)*t+r(1,2)*v - eta(1)).*(r(1,1)*t+r(1,2)*v - eta(1)) + (r(2,2)*v-eta(2)).*(r(2,2)*v-eta(2)) + f-eta.?*eta)),0,1,0,1,?AbsTol?,eps,?RelTol?,eps) - i2
%
% %(1/(r(1,1)*r(2,2)))*integral2(@(t,v) exp(-0.5*((r(1,1)*t+r(1,2)*v - eta(1)).*(r(1,1)*t+r(1,2)*v - eta(1)) + (r(2,2)*v-eta(2)).*(r(2,2)*v-eta(2)) + f-eta.?*eta)),0,r(1,1)+r(1,2),0,r(2,2),?AbsTol?,eps,?RelTol?,eps) - i2
%
% (exp(-(f+eta.?*eta)/2)/(r(1,1)*r(2,2)))*integral2(@(t,v) exp(-0.5*(t.*t + v.*v)),-eta(1),r(1,1)+r(1,2)-eta(1),-eta(2),r(2,2)-eta(2),?AbsTol?,eps,?RelTol?,eps) - i2
%
%
% integral2(@(t,v) exp(-0.5*(a*t.*t -b*t.*v-b*v.*t+c*v.*v + 2*d*t - 2*e*v + f)),0,1,0,1,?AbsTol?,eps,?RelTol?,eps) - i2

%First we determine if a*c - b*b = 0



% abs(r(2,2)) %% cant just check R22, need to check R22*R11?
R = triu(qr([a, b;b c]));
abc = abs(det(R));
if abs(abc) > eps*10


    %     abc =  a*c-b*b;
    c11 =  c/abc;
    c12 =  -b/abc;
    c22 =  a/abc;

    mu1 = c11*d + c12*e;
    mu2 = c12*d + c22*e;

    %i3 = 2*pi*sqrt(c11*c22-c12*c12)*exp(-0.5*f+0.5*(mu1*a*mu1+2*mu1*mu2*b+c*mu2*mu2));
    i30 = 2*pi*sqrt(c11*c22-c12*c12)*exp(-0.5*f+0.5*(d*d*c11 + 2*d*e*c12 + e*e*c22));
    ibvn = bvn(-mu1/sqrt(c11)  ,  (1-mu1)/sqrt(c11)  ,  -mu2/sqrt(c22)  ,  (1-mu2)/sqrt(c22)  ,  c12/sqrt(c11*c22));
    i3 = i30*ibvn;
    %
    %     xl = -mu1/sqrt(c11);
    %     xu = (1-mu1)/sqrt(c11);
    %     yl = -mu2/sqrt(c22);
    %     yu = (1-mu2)/sqrt(c22);
    %
    %     rho = c12/sqrt(c11*c22);
    %
    %     ibvn2 = (1/(2*pi*sqrt(1-rho^2)))*integral2(@(t,v) exp(-(t.^2-2*rho*t.*v+v.^2)/(2*(1-rho^2))),xl,xu,yl,yu,'AbsTol',eps,'RelTol',eps)
    %
    %

    %     if abs(i4-i3) > 1e-5
    %         keyboard
    %     end

    %    eta = r\(r.'\[-d;e]);
    %    sc    = sqrt(1+(r(1,2)/r(2,2))^2);
    %    isc11 = abs(r(1,1))/sc;
    %    isc22 = abs(r(2,2));
    %
    %    mu1 = eta(1);
    %    mu2 = eta(2);
    %
    %    i3 = (2*pi/(abs(r(1,1)*r(2,2))))*exp(-0.5*f+0.5*(mu1*a*mu1+2*mu1*mu2*b+c*mu2*mu2));
    %    i3 = i3*bvn(-mu1*isc11,(1-mu1)*isc11,-mu2*isc22,(1-mu2)*isc22,sign(r(1,1))*r(1,2)/sqrt(r(2,2)^2+r(1,2)^2));


    %    if abs(i2-i3)>1e-2
    %        keyboard
    %    end

else
    d = -d;
    b = -b;
    if abs(a) < eps*10
        %keyboard
        if abs(c) > eps*10
            i3 = sqrt(pi/2/c)*exp((e^2/c-f)/2)*(erf(sqrt(c/2)-e/sqrt(2*c))-erf(-e/sqrt(2*c)));
        else %a = c = 0
            i3 = exp(-f/2);
        end
    elseif abs(c) < eps*10
        i3 = sqrt(pi/2/a)*exp((d^2/a-f)/2)*(erf(sqrt(a/2)+d/sqrt(2*a))-erf(+d/sqrt(2*a)));
    else %a=b=c ~= 0
        %determine g s.t. a = g*c
        g   = sqrt(a/c);

        sa  = sqrt(a/2);
        spi = sqrt(pi);
        if abs(g-1) < 100*eps % need to check this case
            i3  = (2*spi/a)*exp(-0.5*f)*(sa*erf(sa)+exp(-a/2)/spi-1/spi);

        else
            cn = (exp(-0.5*f)*spi/2/g/sa);
            i3 = (1/sa)*(sa*erf(sa)+exp(-a/2)/spi - ((1-g)*sa*erf((1-g)*sa)+exp(-((sa-g*sa))^2)/spi));
            i4 = (1/sa)*(1/spi-(-g*sa*erf(-g*sa) + exp(-(g*sa)^2)/spi));
            i3 = cn*(i3-i4);

            %            beta = b/a;
            %            i5 = (exp(-f/2)/beta) * integral2(@(t,v) exp(-(1/2/a)*(t-v).^2),0,beta,0,1,?AbsTol?,eps,?RelTol?,eps);
            %
            %            if abs(i2-i3)>1e-2
            %                keyboard
            %            end
        end
    end
end