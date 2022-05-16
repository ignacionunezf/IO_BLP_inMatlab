clear;  clear global;

%Read data

filename = 'Dataset.xlsx';
sheet = 1;
xlRange = 'A2:G351';
[num,text,raw] = xlsread(filename, sheet, xlRange);

data=num;
markets=70;
products=max(data(:,2));
obs=size(data,1);
X=zeros(obs,5);
X(:,1:4)=data(:,3:6);
%columns:
% 1: Market ID
% 2: Product ID
% 3: Characteristics 1: dummy variable
% 4: Characteristics 2: cts variable 
% 5: Characteristics 3: cts variable
% 6: Price
% 7: Market share

%We proceed just like in the logits case. First, we calculate the outside good market share

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

%Second, we calculate the within market share
within_share=zeros(obs,1);

for i=1:obs
    
   within_share(i)=data(i,7)/(1-outside_share(data(i,1),1));
   X(i,5)= log(within_share(i));
    
end

%Like before, I define the mean utility for each product

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

%Third, I get rid of the unobserved tsi_j by calculating the average delta
%for each product across markets, and substrating that average from the
%previous deltas; like a fixed-effect estimator when using panel data

delta_avg=zeros(products,1);
X_avg=zeros(products,5);

for p=1:products
    
    d=0;
    x=zeros(1,5);
    n=0;
    for i=1:obs
    if delta(i,2)==p
        
    d=d+delta(i,3);
    x(1,1:4)=x(1,1:4)+data(i,3:6);
    x(1,5)=x(1,5)+log(within_share(i));
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
%don't change accross observartionsm, so they are unidentifiable with this
%method. Thus, I regress delta_bar on the demeaned price x

y=delta_bar(:,3);
x=X(:,4:5);

sum1=zeros(2,2);
sum2=zeros(2,1);

for i=1:obs
    
sum1=sum1+transpose(X(i,4:5))*X(i,4:5);
sum2=sum2+transpose(X(i,4:5))*y(i); 

end

alpha=inv(sum1)*sum2;
e=y-x*alpha;
error_sq=dot(e,e);
var=error_sq/(obs-2);
V=var*inv(sum1);
SE_alpha1=sqrt(V(1,1));
SE_alpha2=sqrt(V(2,2));

%Now, having the estimate for the price, we can calculate the share
%elasticities with respect to prices. I'll do this numerically


for m=1:markets
        
    for i=1:obs
    if data(i,1)==m
        
    delta(i,3)=log(data(i,7))-log(outside_share(m,1))-alpha(2)*log(within_share(i));
    
    end
    end
    
end


%First, I calculate the tsi plus x*b

prices=zeros(markets,products);
shares=zeros(markets,products);
new_shares=zeros(markets,products);
deltas=zeros(markets,products);
residuals=zeros(markets,products);

for m=1:markets
    for p=1:products
    for i=1:obs
        if (data(i,1)==m)&(data(i,2)==p)
            prices(m,p)=data(i,6);
            shares(m,p)=data(i,7);
            deltas(m,p)=delta(i,3);
            residuals(m,p)=delta(i,3)-data(i,6)*alpha(1);
        end    
    end
    end
end

elasticities=zeros(markets,products,products);
d_shares=zeros(markets,products,products);
deltas_plus=deltas;
step=0.00001;

for m=1:markets
   for p1=1:products
       
   deltas_plus=deltas;
   deltas_plus(m,p1)=residuals(m,p1)+prices(m,p1)*alpha(1)*(1+step);
   Dg=0;
   
   for p2=1:products
       if shares(m,p2)>0&shares(m,p1)>0
   Dg=Dg+exp(deltas_plus(m,p2)/(1-alpha(2)));
       end
   end
   Dg_sigma=Dg^alpha(2);
   Dg_sum=Dg^(1-alpha(2))+1;
   
   for p2=1:products
       if shares(m,p2)>0&shares(m,p1)>0
   new_shares(m,p2)=exp(deltas_plus(m,p2)/(1-alpha(2)))/(Dg_sigma*Dg_sum); 
       end
   end
   
   for p2=1:products
       if shares(m,p2)>0&shares(m,p1)>0
   elasticities(m,p1,p2)=(new_shares(m,p2)-shares(m,p2))/(step*shares(m,p2));
   d_shares(m,p1,p2)=(new_shares(m,p2)-shares(m,p2))/(step*prices(m,p1));
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

%%%%%%%%%%%%%%%%%MARK-UPS

markups=zeros(markets,products);

for m=1:markets
    for p=1:products
    for i=1:obs
        if (data(i,1)==m)&(data(i,2)==p)
            markups(m,p)=-data(i,7)/d_shares(m,p,p);    
        end    
    end
    end
end

%%%%%%%%%%%%%%%%% Marginal Costs

mc=prices-markups;
marginalcost=zeros(obs,1);


for m=1:markets
    for p=1:products
    for i=1:obs
        
        if (data(i,1)==m)&(data(i,2)==p)
            marginalcost(i)=mc(m,p);    
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

beta=inv(sum1)*sum2;
e=y-x*beta;
error_sq=dot(e,e);
var=error_sq/(obs-2);
V=var*inv(sum1);
SE_beta1=sqrt(V(1,1));
SE_beta2=sqrt(V(2,2));
SE_beta3=sqrt(V(3,3));
SE_beta4=sqrt(V(4,4));


