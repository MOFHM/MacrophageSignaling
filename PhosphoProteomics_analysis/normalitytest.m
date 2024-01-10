function Results = normalitytest1(x)


%Kolmogorov-Smirnov
n=length(x);
i=1:n;
y=sort(x);
fx=normcdf(zscore(y));
dplus=max(abs(fx-i/n));
dminus=max(abs(fx-(i-1)/n));
Dn=max(dplus,dminus);
KSz=sqrt(n)*Dn;
s=-20:1:20;
a=(-1).^s.*exp(-2*(s.*KSz).^2); 
pvalue=1-sum(a);
Results(1,1)=KSz;
Results(1,2)=pvalue;


%Anderson-Darling
adj=1+0.75/n+2.25/(n^2);
i=1:n;
ui=normcdf(zscore(y),0,1);
oneminusui=sort(1-ui);
lastt=(2*i-1).*(log(ui)+log(oneminusui));
asquare=-n-(1/n)*sum(lastt);
AD=asquare*adj;
 
if AD<=0.2
    pvalue=1-exp(-13.436+101.14*AD-223.73*AD^2);
elseif AD<=0.34
    pvalue=1-exp(-8.318+42.796*AD-59.938*AD^2);
elseif AD<=0.6
    pvalue=exp(0.9177-4.279*AD-1.38*AD^2);
elseif AD<=153.467
    pvalue=exp(1.2937*AD-5.709*AD+0.0186*AD^2);
else
    pvalue=0;
end
Results(2,1)=AD;
Results(2,2)=pvalue;

%Cramer-Von Mises
adj=1+0.5/n;
i=1:n;
fx=normcdf(zscore(y),0,1);
gx=(fx-((2*i-1)/(2*n))).^2;
CvMteststat=(1/(12*n))+sum(gx);
AdjCvM=CvMteststat*adj;
 
if AdjCvM<0.0275
    pvalue=1-exp(-13.953+775.5*AdjCvM-12542.61*(AdjCvM^2));
elseif AdjCvM<0.051
    pvalue=1-exp(-5.903+179.546*AdjCvM-1515.29*(AdjCvM^2));
elseif AdjCvM<0.092
    pvalue=exp(0.886-31.62*AdjCvM+10.897*(AdjCvM^2));
elseif AdjCvM>=0.093
    pvalue=exp(1.111-34.242*AdjCvM+12.832*(AdjCvM^2));
end
Results(3,1)=AdjCvM;
Results(3,2)=pvalue;

%Shapiro-Wilk
a=[];
i=1:n;
mi=norminv((i-0.375)/(n+0.25));
u=1/sqrt(n);
m=mi.^2;
 
a(n)=-2.706056*(u^5)+4.434685*(u^4)-2.07119*(u^3)-0.147981*(u^2)+0.221157*u+mi(n)/sqrt(sum(m));
a(n-1)=-3.58263*(u^5)+5.682633*(u^4)-1.752461*(u^3)-0.293762*(u^2)+0.042981*u+mi(n-1)/sqrt(sum(m));
a(1)=-a(n);
a(2)=-a(n-1);
eps=(sum(m)-2*(mi(n)^2)-2*(mi(n-1)^2))/(1-2*(a(n)^2)-2*(a(n-1)^2));
a(3:n-2)=mi(3:n-2)./sqrt(eps);
    ax=a.*y;
    KT=sum((x-mean(x)).^2);
    b=sum(ax)^2;
    SWtest=b/KT;
mu=0.0038915*(log(n)^3)-0.083751*(log(n)^2)-0.31082*log(n)-1.5861;
sigma=exp(0.0030302*(log(n)^2)-0.082676*log(n)-0.4803);
z=(log(1-SWtest)-mu)/sigma;
pvalue=1-normcdf(z,0,1);
Results(4,1)=SWtest;
Results(4,2)=pvalue;

%Jarque-Bera
E=skewness(y);
B=kurtosis(y);
JBtest=n*((E^2)/6+((B-3)^2)/24);
pvalue=1-chi2cdf(JBtest,2);
Results(5,1)=JBtest;
Results(5,2)=pvalue;

%D'Agostino
beta2=(3*(n^2+27*n-70)*(n+1)*(n+3))/((n-2)*(n+5)*(n+7)*(n+9));
wsquare=-1+sqrt(2*(beta2-1));
delta=1/sqrt(log(sqrt(wsquare)));
alfa=sqrt(2/(wsquare-1));
 
expectedb2=(3*(n-1))/(n+1);
varb2=(24*n*(n-2)*(n-3))/(((n+1)^2)*(n+3)*(n+5));
sqrtbeta=((6*(n^2-5*n+2))/((n+7)*(n+9)))*sqrt((6*(n+3)*(n+5))/(n*(n-2)*(n-3)));
A=6+(8/sqrtbeta)*(2/sqrtbeta+sqrt(1+4/(sqrtbeta^2)));
 
squarerootb=skewness(y);
Y=squarerootb*sqrt(((n+1)*(n+3))/(6*(n-2)));
zsqrtbtest=delta*log(Y/alfa+sqrt((Y/alfa)^2+1));
 
b2=kurtosis(y);
zet=(b2-expectedb2)/sqrt(varb2);
ztestb2=((1-2/(9*A))-((1-2/A)/(1+zet*sqrt(2/(A-4))))^(1/3))/sqrt(2/(9*A));
 
DAPtest=zsqrtbtest^2+ztestb2^2;
 
pvalue=1-chi2cdf(DAPtest,2);
Results(6,1)=DAPtest;
Results(6,2)=pvalue;

Results(7,1) = skewness(x);
Results(7,2) = skewness(x);
Results(8,1) = kurtosis(x);
Results(8,2) = kurtosis(x);

%D'Angostino 2
alpha = 0.05;
x = sort(x);
x = x(:);
n = length(x);
[c,v]=hist(x,x);  %record of data in a frequency table form
nc=find(c~=0);
c=[v(nc) c(nc)'];
x = c(:,1);
f = c(:,2);
s1 = f'*x;
s2 = f'*x.^2;
s3 = f'*x.^3;
s4 = f'*x.^4;
SS = s2-(s1^2/n);
v = SS/(n-1);
k3 = ((n*s3)-(3*s1*s2)+((2*(s1^3))/n))/((n-1)*(n-2));
g1 = k3/sqrt(v^3);
k4 = ((n+1)*((n*s4)-(4*s1*s3)+(6*(s1^2)*(s2/n))-((3*(s1^4))/(n^2)))/((n-1)*(n-2)*(n-3)))-((3*(SS^2))/((n-2)*(n-3)));
g2 = k4/v^2;
eg1 = ((n-2)*g1)/sqrt(n*(n-1));  %measure of skewness
eg2 = ((n-2)*(n-3)*g2)/((n+1)*(n-1))+((3*(n-1))/(n+1));  %measure of kurtosis
A = eg1*sqrt(((n+1)*(n+3))/(6*(n-2)));
B = (3*((n^2)+(27*n)-70)*((n+1)*(n+3)))/((n-2)*(n+5)*(n+7)*(n+9));
C = sqrt(2*(B-1))-1;
D = sqrt(C);
E = 1/sqrt(log(D));
F = A/sqrt(2/(C-1));
Zg1 = E*log(F+sqrt(F^2+1));
G = (24*n*(n-2)*(n-3))/((n+1)^2*(n+3)*(n+5));
H = ((n-2)*(n-3)*g2)/((n+1)*(n-1)*sqrt(G)); %We thank A.M. Winkler for 
%having indicated this bug
%H = ((n-2)*(n-3)*abs(g2))/((n+1)*(n-1)*sqrt(G));
J = ((6*(n^2-(5*n)+2))/((n+7)*(n+9)))*sqrt((6*(n+3)*(n+5))/((n*(n-2)*(n-3))));
K = 6+((8/J)*((2/J)+sqrt(1+(4/J^2))));
L = (1-(2/K))/(1+H*sqrt(2/(K-4)));
Zg2 = (1-(2/(9*K))-L^(1/3))/sqrt(2/(9*K));
K2 = Zg1^2 + Zg2^2;  %D'Agostino-Pearson statistic
X2 = K2;  %approximation to chi-distribution
df = 2;  %degrees of freedom
P = 1-chi2cdf(X2,df);
Results(9,1) = X2;
Results(9,2) = P;


end
