function [V,pol,money,payment,size,grid_m,grid_s]=...
    SolveModel(beta,eta_h,eta_l,R,sigma,kappa,phi,mu_unc,sigma_unc,N,T)

rng(20);
% Grids
nmm     = 300;
nss     = 300;
grid_s  = linspace(0,250,nss);
grid_m  = linspace(0,250,nmm);

% Unconditional and conditional distributions for s and s'

ds                 = (grid_s(end)-grid_s(1))/(nss-1);

ucd_distr        = zeros(nss,1);
ucd_distr(1)     = logncdf(grid_s(1)+ds/2,mu_unc,sigma_unc);
for i=2:nss-1
    ucd_distr(i) = logncdf(grid_s(i)+ds/2,mu_unc,sigma_unc)-logncdf(grid_s(i)-ds/2,mu_unc,sigma_unc);
end
ucd_distr(nss)   = 1-logncdf(grid_s(nss)-ds/2,mu_unc,sigma_unc);


cd_distr            = zeros(nss,nss);
cd_distr(:,1)       = normcdf((grid_s(1)+ds/2-grid_s)/sigma);
cd_distr(:,2:end-1) = normcdf((grid_s(2:end-1)+ds/2-grid_s')/sigma)-...
    normcdf((grid_s(2:end-1)-ds/2-grid_s')/sigma);
cd_distr(:,end)     = 1-normcdf((grid_s(end)-ds/2-grid_s)/sigma);

% Value functions
Vtilde  = zeros(nmm,nss);
Vadj    = zeros(nmm,nss);
Vnadj   = zeros(nmm,nss);
V       = zeros(nmm,nss);
TV      = zeros(nmm,nss);
TVA     = zeros(nmm,nss);
TVN     = zeros(nmm,nss);
Vcard   = zeros(nmm,nss);
Vcash   = zeros(nmm,nss);

% Setup VFI
dist=10;
tol=10^(-4);
iter=0;
maxiter=1000;

% VFI
while dist>tol && iter<maxiter
    iter=iter+1;
    V=TV;
    Vtilde=R.*grid_m'+min(V*cd_distr, min(V*cd_distr)+eta_l);
    for m=1:nmm
        for s=1:nss
            if s<=m
            TVA(m,s)=min(kappa+R*grid_m(m)+beta*(Vtilde(m,:)*ucd_distr)...
                    ,R*(grid_m(m))+beta*(Vtilde(m-s+1,:)*ucd_distr));
            TVN(m,s)= R*(grid_m(m))+beta*(Vtilde(m-s+1,:)*ucd_distr);
            else
            TVA(m,s)=kappa+R*grid_m(m)+beta*(Vtilde(m,:)*ucd_distr);         
            TVN(m,s)=eta_h+R*grid_m(m)+min(beta*(Vtilde*ucd_distr));
            end
        end
    end
    TV=(phi.*TVA+(1-phi).*TVN);
    dist=norm(TV-V)/norm(V);
end
post_adj_optpos_raw=zeros(nmm,nss);
for m=1:nmm
    for s=1:nss
        Vnadj(m,s) = V(m,:)*cd_distr(s,:)';
        [Vadj(m,s),post_adj_optpos_raw(m,s)]  = min(V*cd_distr(s,:)'+eta_l);
        if s>m
            Vcard(m,s)=kappa+R*grid_m(m)+beta*(Vtilde(m,:)*ucd_distr);
            Vcash(m,s)=eta_h+R*grid_m(m)+min(beta*(Vtilde*ucd_distr));
        else
            Vcard(m,s)=kappa+R*grid_m(m)+beta*(Vtilde(m,:)*ucd_distr);
            Vcash(m,s)=R*(grid_m(m))+beta*(Vtilde(m-s+1,:)*ucd_distr);
        end
    end
end
[~,forced_adj]      = min(Vtilde*ucd_distr);
pol.post_forced_adj = grid_m(forced_adj);
pol.adj_yn          = zeros(nmm,nss);
pol.adj_size        = zeros(nmm,nss);
post_adj_optpos     = zeros(nmm,nss);
for m=1:nmm
    for s=1:nss
        if Vadj(m,s)<Vnadj(m,s)
            pol.adj_yn(m,s)      = 1;
            post_adj_optpos(m,s) = post_adj_optpos_raw(m,s);
            pol.adj_size(m,s)    = grid_m(post_adj_optpos(m,s))-grid_m(m);
        end
    end
end

pol.cash=zeros(nmm,nss);
for m=1:nmm
    for s=1:m
        if Vcash(m,s)<Vcard(m,s)
            pol.cash(m,s)=1;
        end
    end
end

money.morning    = zeros(N,T);
money.afternoon  = zeros(N,T);
money.aft_adj    = zeros(N,T);
money.path       = zeros(N,2*T);
money.adjust     = zeros(N,T);
size.morning     = zeros(N,T);
size.afternoon   = zeros(N,T);
payment.method   = zeros(N,T);
payment.cash     = zeros(N,T);
payment.cashless = zeros(N,T);
payment.all      = zeros(N,T);
payment.type     = zeros(N,T);    %1 voluntary choice, 0 forced cash/card

for n=1:N
    % Initialize day one: random initial money and signal
    money.morning(n,1)  = randsample(grid_m, 1, true);
    size.morning(n,1)   = randsample(grid_s, 1, true, ucd_distr);
    
    % Afternoon day one: adjusted money and true size
    
    % Adjusted wallet
    money.afternoon(n,1) = money.morning(n,1)+...
        pol.adj_size(find(grid_m==money.morning(n,1)),...
        find(grid_s==size.morning(n,1)));
    money.adjust(n,1)=money.afternoon(n,1)-money.morning(n,1);
    % Realized expenditure size
    size.afternoon(n,1) = randsample(grid_s,1,true,...
        cd_distr(find(grid_s==size.morning(n,1)),:));
    
    % Other days : same structure
    for t=2:T
        
        % Signal and true size of day t>2
        size.morning(n,t)  = randsample(grid_s, 1, true, ucd_distr);
        size.afternoon(n,t)= randsample(grid_s,1,true,cd_distr(find(grid_s==size.morning(n,t)),:));
        
        % Which kind of shop visited previous afternoon?
        
        
        if rand<=phi                             % Accepted cashless
            
            % If accepted cashless, cash usage policy function needed.
            money.morning(n,t)   = money.afternoon(n,t-1)-...
                size.afternoon(n,t-1)*pol.cash(find(abs(grid_m-money.afternoon(n,t-1))<tol),...
                find(abs(grid_s-size.afternoon(n,t-1))<tol));
            
            % Get payment choice of previous day
            payment.method(n,t-1)=pol.cash(find(abs(grid_m-money.afternoon(n,t-1))<tol),...
                find(abs(grid_s-size.afternoon(n,t-1))<tol));
            payment.type(n,t-1)=(money.afternoon(n,t-1)>size.afternoon(n,t-1));
        else                                     % Did not accept cashless
            payment.type(n,t-1)=0;
            % If had enough cash, simply used that.
            if money.afternoon(n,t-1)>=size.afternoon(n,t-1)
                money.morning(n,t)= money.afternoon(n,t-1)-size.afternoon(n,t-1);
            % Otherwise, forced withdrawal.           
            else
                money.morning(n,t)= pol.post_forced_adj;
            end
            % Of course paid with cash.
            payment.method(n,t-1)=1;
        end 
        
        % Afternoon adjustment in day t
        money.afternoon(n,t) = money.morning(n,t)+...
            pol.adj_size(find(abs(grid_m-money.morning(n,t))<tol),find(abs(grid_s-size.morning(n,t))<tol));
        money.adjust(n,t)=money.afternoon(n,t)-money.morning(n,t);
    end
    % Payment choice of last day 
    if rand<phi
    payment.method(n,T)=pol.cash(find(abs(grid_m-money.afternoon(n,T))<tol),...
            find(abs(grid_s-size.afternoon(n,T))<tol)); 
    else
        payment.method(n,T)=1;
    end
    for t=1:T
        payment.all(n,t)          = size.afternoon(n,t);
        if payment.method(n,t)==1
            payment.cash(n,t)     = size.afternoon(n,t);
        else
            payment.cashless(n,t) = size.afternoon(n,t);
        end
        if money.afternoon(n,t)~=money.morning(n,t)
            money.aft_adj(n,t)=money.afternoon(n,t);
        end
    end
    for t=1:2:2*T-1
        money.path(n,t)   = money.morning(n,(t+1)/2);
        money.path(n,t+1) = money.afternoon(n,(t+1)/2);
    end    
end
