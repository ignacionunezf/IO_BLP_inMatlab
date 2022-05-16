function [ret,newshares,markups]=dist2(input);

global alpha
global data
global residuals
global marginalcost
global shares
global markets
global products
global obs
global subsidy
global m1
global m2


deltas=zeros(markets,products);
newshares=zeros(markets,products);
prices=zeros(markets,products);


for m=m1:m2
    for p1=1:products
    for i=1:obs
        if (data(i,1)==m)&(data(i,2)==p1)
            prices(m,p1)=input(i)-subsidy/1000;
        end    
    end
    end
end


for m=m1:m2
   for p1=1:products
 deltas(m,p1)=residuals(m,p1)+prices(m,p1)*alpha;
   end
end

for m=m1:m2
    
   x=0;
   for p1=1:products
       if shares(m,p1)>0
   x=x+exp(deltas(m,p1));
       end
   end 
     
   for p1=1:products   
       if shares(m,p1)>0
   newshares(m,p1)=exp(deltas(m,p1))/(1+x);
       end
   end
    
end

s=zeros(obs,1);

for m=m1:m2
    for p1=1:products
    for i=1:obs        
        if (data(i,1)==m)&(data(i,2)==p1)
            s(i)=newshares(m,p1);    
        end    
    end
    end
end

d_shares=zeros(obs,1);
markups=zeros(obs,1);

    for i=1:obs
            d_shares(i)=alpha*(1-s(i))*s(i);
            markups(i)=-s(i)/d_shares(i);    
    end


prices2=marginalcost+markups;

ret=0;

for i=1:obs
    if (data(i,1)>=m1)&(data(i,1)<=m2)
    ret=ret+abs(prices2(i)-input(i));
    end
end


end

