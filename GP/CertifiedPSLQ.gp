\\ Implementation of the short PSLQ algorithm over Q
\\ Based on arXiv:1809.06063v1 [math.NT]
\\ Extending the PSLQ algorithm to algebraic integer relations by Matthew Skerritt and Paul Vrbik
\\ Code written by: Abhiram Kidambi
\\ Updated 13 May 2026

issqmat(M)={
return(type(M)=="t_MAT"&&matsize(M)[1]==matsize(M)[2]);
};

istrapezoidal(M)={
my(sz=matsize(M));
return(sz[1]>sz[2]);
};

diagonalLTM(M)={
my(sz=matsize(M),r);
r=min(sz[1],sz[2]);
return(vector(r,i,M[i,i]));
};

rowswap(M,m)={
my(tmp);
tmp=M[m,];
M[m,]=M[m+1,];
M[m+1,]=tmp;
return(M);
};

colswap(M,m)={
my(tmp);
tmp=M[,m];
M[,m]=M[,m+1];
M[,m+1]=tmp;
return(M);
};

nint(x)={
return(round(x));
};

selm(gamma,H)={
my(sz,ncols,vals,m);
sz=matsize(H);
ncols=sz[2];
vals=vector(ncols,i,gamma^i*abs(H[i,i]));
vecmax(vals,&m);
return(m);
};

buildH(n,s,y)={
my(H);
H=matrix(n,n-1);
for(i=1,n,
for(j=1,n-1,
if(i<j,
H[i,j]=0,
if(i==j,
H[i,j]=s[i+1]/s[i],
H[i,j]=-y[i]*y[j]/(s[j]*s[j+1])
)
);
);
);
return(H);
};

hermitered(H,A,B,{eps=0})={
my(n,q);
n=matsize(A)[1];
if(eps==0,
eps=10^(-default(realprecision)+20);
);
for(i=2,n,
forstep(j=i-1,1,-1,
if(abs(H[j,j])<eps,
next;
);
q=round(H[i,j]/H[j,j]);
if(q==0,
next;
);
for(k=1,j,
H[i,k]-=q*H[j,k];
);
A[i,]-=q*A[j,];
B[,j]+=q*B[,i];
);
);
return([H,A,B]);
};

fixH(H,m)={
my(a,b,r,c,s,h0,h1);
a=H[m,m];
b=H[m,m+1];
r=sqrt(a^2+b^2);
if(r==0,
return(H);
);
c=a/r;
s=b/r;
for(i=m,matsize(H)[1],
h0=H[i,m];
h1=H[i,m+1];
H[i,m]=c*h0+s*h1;
H[i,m+1]=-s*h0+c*h1;
);
return(H);
};

pslqtau(gamma)={
my(gcrit,tau);
gcrit=sqrt(4/3);
if(gamma<=gcrit,
error("pslqtau: gamma must be strictly larger than sqrt(4/3).");
);
tau=1/sqrt(1/gamma^2+1/4);
return(tau);
};

pslqiterbound(vec,Mrel,{gamma=1.2})={
my(n,tau,N);
n=#vec;
if(n<2,
error("pslqiterbound: dimension of vec must be at least 2.");
);
if(Mrel<=0,
error("pslqiterbound: relation norm bound Mrel must be positive.");
);
tau=pslqtau(gamma);
N=binomial(n,2)*log(gamma^(n-1)*Mrel)/log(tau);
return(max(1,ceil(N)));
};

pslq(vec,{itermax=1000},{flag=0},{gamma=sqrt(4/3)+0.01},{tol=0},{verbose=1})={
my(n,svec,yvec,H,A,B,m,iter);
my(res,normv,eps,residual,relation);
n=#vec;
iter=0;
if(n<2,
error("pslq: vector must have length at least 2.");
);
if(gamma<=sqrt(4/3),
error("pslq: gamma must be strictly larger than sqrt(4/3).");
);
if(tol==0,
tol=10^(-floor(0.75*default(realprecision)));
);
eps=10^(-default(realprecision)+20);
normv=sqrt(vec*vec~);
if(normv==0,
error("pslq: zero input vector.");
);
vec=vec/normv;
svec=vector(n,k,sqrt(sum(j=k,n,vec[j]^2)));
yvec=vec/svec[1];
svec=svec/svec[1];
H=buildH(n,svec,yvec);
A=matid(n);
B=matid(n);
while(iter<itermax,
if(verbose,
printf("%sPSLQ iteration %d/%d",Strchr(13),iter+1,itermax);
);
res=hermitered(H,A,B,eps);
H=res[1];
A=res[2];
B=res[3];
m=selm(gamma,H);
H=rowswap(H,m);
A=rowswap(A,m);
B=colswap(B,m);
if(m<n-1,
H=fixH(H,m);
);
for(j=1,n,
relation=B[,j];
residual=abs(vec*relation);
if(residual<tol,
if(verbose,
printf("%sRelation found at PSLQ iteration %d/%d, column %d\n",Strchr(13),iter+1,itermax,j);
print("Evaluation of relation = ",residual);
print("Relation = ");
print(relation);
);
if(flag==1,
return(1),
return(relation)
);
);
);
iter++;
);
if(verbose,
print("");
);
if(flag==1,
return(0),
error("pslq: no relation found within itermax.");
);
};

pslqgamma(vec,Mrel,{flag=0},{gamma=1.2},{tol=0},{verbose=1})={
my(itermax);
if(gamma<=sqrt(4/3),
error("pslqgamma: gamma must be strictly larger than sqrt(4/3).");
);
itermax=pslqiterbound(vec,Mrel,gamma);
return(pslq(vec,itermax,flag,gamma,tol,verbose));
};

pslqbounded(vec,Mrel,{flag=0},{gamma=1.2},{tol=0},{verbose=1})={
my(N);
N=pslqiterbound(vec,Mrel,gamma);
print("Dimension n       = ",#vec);
print("Relation norm M   = ",Mrel);
print("gamma             = ",gamma);
print("tau               = ",pslqtau(gamma));
print("Iteration bound N = ",N);
return(pslq(vec,N,flag,gamma,tol,verbose));
};

