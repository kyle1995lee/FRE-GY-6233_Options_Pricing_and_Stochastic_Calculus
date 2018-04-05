%Calculation of Asian options and using the variance reduction technique
%2018/04/04 Kyle

%Set Variables
times=1000000;m=12;K=100;T=1;r=0.01;sigma=0.2;
%times:simulation times, m:number of prices in a year, T:year, 
%r:interest rates, sigma:volatility, assume known
S = NaN(12,1);S(1)=110;Sarith=S(1);Sgeo=S(1);Smean=0;
C1=NaN(times,1);C2=NaN(times,1);
C1total=0;C2total=0;
temp=0;temp2=0;temp_sc1=0;temp_sc2=0;
z=NaN(12,1);temp_exp=NaN(times,1);

%Calculating Arithmetic mean and geometric mean
%Solution of Black Scholes using Monte Carlo simulations
for i=1:times
    %Setting starting point
    Sarith=S(1);
    Sgeo=S(1);
    %create an array of standard normal numbers to simulate and plugging it
    %into the black scholes formula, save all in S() till maturity
    for j=1:m
        z(j)=normrnd(0,1);
        S(j+1)=S(j)*exp((r-0.5*sigma^2)*(1/m)+sigma*sqrt(1/m)*z(j));
    end
    Sarithmean=mean(S);
    Sgeomean=geomean(S);
    %C1:Asian Call option price using arithmetic mean
    %C2:Asian Call option price using geometric mean
    C1(i)=exp(-r*T)*max(Sarithmean-K,0);
    C2(i)=exp(-r*T)*max(Sgeomean-K,0);
    C1total=C1total+C1(i);
    C2total=C2total+C2(i);
end
%Average call price
C1avg=C1total/times;
C2avg=C2total/times;

%Calculating b, a number to apply the control variates method
for i=1:times
    temp=(C1(i)-C1avg)*(C2(i)-C2avg)+temp;
    temp2=(C2(i)-C2avg)*(C2(i)-C2avg)+temp2;
end
b=temp/temp2;

%Plug everything into Black Scholes to get Tbar, sigma squared, and delta
Tbar=(1/m)*sum(1/12:1/12:1);
temp_var_norm=0;
for i=1:m
    temp_var_norm=temp_var_norm+(2*i-1)*((m+1-i)/m);
end
sigma_sq=(sigma^2/(m^2*T))*temp_var_norm;
delta=0.5*sigma^2-0.5*sigma_sq;
%getting the d+ in N(d+) in the Black Scholes formula
d=(log(S(1)/K)+(r-delta+0.5*sigma_sq)*Tbar)/(sqrt(sigma_sq*Tbar));
%Achieve the expected value to implement the method
expected=exp(-delta*Tbar-r*(T-Tbar))*S(1)*normcdf(d)-exp(-r*T)*K*normcdf(d-sqrt(sigma_sq*Tbar));

%getting the Ybar(b)'s average
for i=1:times
    temp_exp(i)=C1(i)-b*(C2(i)-expected);
end
y_avg=sum(temp_exp)/times;

%calculating standard errors for Arithmetic and Geometric means
for i=1:times
    temp_sc1=temp_sc1+(temp_exp(i)-y_avg)^2;
    temp_sc2=temp_sc2+(C2(i)-C2avg)^2;
end
samp_std_dev1=sqrt(temp_sc1/(times-1));
samp_std_dev2=sqrt(temp_sc2/(times-1));
std_err_1=samp_std_dev1/sqrt(times);    
std_err_2=samp_std_dev2/sqrt(times);
%calculating correlation between Arithmatic and Geometric
temp_cov=cov(C1,C2);
corr=temp_cov(2,1)/sqrt(var(C1)*var(C2));
