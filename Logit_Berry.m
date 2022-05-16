clear; clear global;

%Read data

filename = 'Dataset.xlsx';
sheet = 1;
xlRange = 'A2:G351';
[num,text,raw] = xlsread(filename, sheet, xlRange);
%columns:
% 1: Market ID
% 2: Product ID
% 3: Characteristics 1: dummy variable
% 4: Characteristics 2: cts variable 
% 5: Characteristics 3: cts variable
% 6: Price
% 7: Market share

global alpha
global data
global residuals
global marginalcost
global shares
global markets
global products
global obs


data=num;
markets=70;
products=max(data(:,2));
obs=size(data,1);
X=data(:,3:6);

%%%%%%%%%%%%%%%%%PARAMETERS


%First, I calculate the outside good market share

outside_share=zeros(70,1);
for m=1:markets
    
    s0=0;
    
    for i=1:obs
    if data(i,1)==m
    s0=s0+data(i,7);
    end
    end
    
    outside_share(m,1)=1-s0;
end

%Second, using the logit formula, I calculate the mean utility of each product in each market

delta=zeros(obs,3); %Market, product, mean utility

for m=1:markets
        
    for i=1:obs
    if data(i,1)==m
        
    delta(i,1)=m;
    delta(i,2)=data(i,2);
    delta(i,3)=log(data(i,7))-log(outside_share(m,1));
    
    end
    end
    
end

%Third, I calculate the mean utility for each product across markets, and substract that average from the
%mean utility; like a fixed-effect estimator when using panel data

delta_avg=zeros(products,1);
X_avg=zeros(products,4);

for p=1:products
    
    d=0;
    x=zeros(1,4);
    n=0;
    for i=1:obs
    if delta(i,2)==p
        
    d=d+delta(i,3);
    x=x+data(i,3:6);
    n=n+1;
    end
    end
    
    delta_avg(p,1)=d/n;
    X_avg(p,:)=x/n;
    
end

%Fourth, I calculate the demeaned deltas and covariates for each
%observation

delta_bar=delta; %Market, product, mean utility

for p=1:products
    
    for i=1:obs
    if delta(i,2)==p
        
    delta_bar(i,3)=delta_bar(i,3)-delta_avg(p,1);
    X(i,:)=X(i,:)-X_avg(p,:);
    end
    end
    
end

%Lastly, we note that the first three covariates (PPO, network, and satisfaction)
%don't change accross observartions, so they are unidentifiable with this
%method. Thus, I regress delta_bar on the demeaned price x

y=delta_bar(:,3);
x=X(:,4);

alpha=dot(y,x)/dot(x,x);
e=y-x*alpha;
error_sq=dot(e,e);
var=error_sq/(obs-2);
V=var/dot(x,x);
SE_alpha=sqrt(V);


%%%%%%%%%%%%%%%%%ELASTICITIES


%Now, having the estimate for the price, we can calculate the share
%elasticities with respect to prices

prices=zeros(markets,products);
shares=zeros(markets,products);

for m=1:markets
    for p1=1:products
    for i=1:obs
        if (data(i,1)==m)&(data(i,2)==p1)
            prices(m,p1)=data(i,6);
            shares(m,p1)=data(i,7);
        end    
    end
    end
end


elasticities=zeros(markets,products,products);

for m=1:markets
    for p1=1:products
    for p2=1:products

        if p1==p2       
   elasticities(m,p1,p2)=alpha*prices(m,p1)*(1-shares(m,p1));
        else if prices(m,p1)>0
   elasticities(m,p1,p2)=-alpha*prices(m,p2)*shares(m,p2);        
            end
        end
        
    end
    end 
end

av_elasticities=zeros(products,products);


for p1=1:products
for p2=1:products
    
count=0;

for m=1:markets
    
if elasticities(m,p1,p2)==0
else
count=count+1;
av_elasticities(p1,p2)=av_elasticities(p1,p2)+elasticities(m,p1,p2);
end
end

av_elasticities(p1,p2)=av_elasticities(p1,p2)/count;

end


end

%%%%%%%%%%%%%%%%% I now recover firms' mark-ups for each market-product using the demand elasticty formula and lerner index

d_shares=zeros(markets,products);
markups=zeros(markets,products);

for m=1:markets
    for p=1:products
    for i=1:obs
        if (data(i,1)==m)&(data(i,2)==p)
            d_shares(m,p)=alpha*(1-data(i,7))*data(i,7);
            markups(m,p)=-data(i,7)/d_shares(m,p);    
        end    
    end
    end
end


%%%%%%%%%%%%%%%%% Now estimate the marginal cost of each product

mc=prices-markups;

marginalcost=zeros(obs,1);

for m=1:markets
    for p=1:products
    for i=1:obs
        
        if (data(i,1)==m)&(data(i,2)==p)
            marginalcost(i)=mc(m,p);    %marginal cost for market-product
        end    
    end
    end
end

y=marginalcost;
x=ones(obs,4);
x(:,2:4)=data(:,3:5);

sum1=zeros(4,4);
sum2=zeros(4,1);

for i=1:obs
    
sum1=sum1+transpose(x(i,:))*x(i,:);
sum2=sum2+transpose(x(i,:))*y(i); 

end

beta=inv(sum1)*sum2;    %linear regression of mc on product variates and an intercept
e=y-x*beta;
error_sq=dot(e,e);
var=error_sq/(obs-2);
V=var*inv(sum1);
SE_beta1=sqrt(V(1,1));
SE_beta2=sqrt(V(2,2));
SE_beta3=sqrt(V(3,3));
SE_beta4=sqrt(V(4,4));

%%%%%%%%%%%%%%%%% Simulate the effect of subsidizing the price and compute new Bertrand-NE

%First, I calculate the mean utility without the price: tsi plus x*b


residuals=zeros(markets,products);
for m=1:markets
    for p=1:products
    for i=1:obs
        if (data(i,1)==m)&(data(i,2)==p)
            residuals(m,p)=delta(i,3)-data(i,6)*alpha;
        end    
    end
    end
end

global subsidy
global m1;
global m2;
solution=zeros(obs,1);

for m=1:70
% Compute NE for each market. 
    %This can can be optimized to the price market by market rather than the full vector  
m1=m;    
m2=m;
subsidy=250;

options = optimset('Display','final','TolFun',1e-8,'TolX',1e-4,'MaxIter',100000);
p=data(:,6); % Starting Values 
test = dist(p);
p_star = fminsearch('dist',p,options);

for i=1:obs   
    if data(i,1)==m 
     
    solution(i)=p_star(i);
    
    end  
end

end

m1=1;    
m2=70;
[dif,newshares,newmarkups]=dist2(solution);

new_s=zeros(obs,1);
for m=1:markets
    for p=1:products
    for i=1:obs
        
        if (data(i,1)==m)&(data(i,2)==p)
            new_s(i)=newshares(m,p);    
        end    
    end
    end
end


new_s0=zeros(70,1);
for m=1:markets
    
    s0=0;
    
    for i=1:obs
    if data(i,1)==m
    s0=s0+new_s(i);
    end
    end
    
    new_s0(m)=1-s0;
end

mean(new_s0)


