#include <oxstd.h>
#include <oxprob.h>
#include <oxfloat.h>
#include <oxdraw.oxh>

dgp1(N, theta, k, c)
{

	/*	This first part creates a correlated Uniform Distribution
	 	Technique taken from
		https://wernerantweiler.ca/blog.php?item=2020-07-05
	*/
	decl C0, C1, E0, E1;
	decl rho = 0.5;
	
	decl A = (sqrt((49+rho)/(1+rho))-5)/2 ;
	decl w1= ranbeta(N,1,A,1) ;
	decl w2= ranbeta(N,1,A,1);  
	decl u1 = ranu(N,1);
	decl u2 = ranu(N,1);
	
	decl v1 = ranu(N,1) .< 0.5 .? fabs(w1 - u1) .: 1- fabs(1-w1-u1);
	decl v2 = ranu(N,1) .< 0.5 .? fabs(w2 - u2) .: 1- fabs(1-w2-u2);
	
	/* 	Finally we create the half-normal draws.
	*/
	
	C0 = fabs(quann(u1));
	C1 = fabs(quann(u2));
	E0 = fabs(quann(v1));
	E1 = fabs(quann(v2))+theta;

	return C0~C1~E0~E1;

}

dgp2(N, theta, k, c)
{

	/* Singh Maddala Distribution*/
	decl C0, C1, E0, E1;
	decl rho = 0.5;
	
	decl A = (sqrt((49+rho)/(1+rho))-5)/2 ;
	decl w1= ranbeta(N,1,A,1) ;
	decl w2= ranbeta(N,1,A,1);  
	decl u1 = ranu(N,1);
	decl u2 = ranu(N,1);
	
	decl v1 = ranu(N,1) .< 0.5 .? fabs(w1 - u1) .: 1- fabs(1-w1-u1);
	decl v2 = ranu(N,1) .< 0.5 .? fabs(w2 - u2) .: 1- fabs(1-w2-u2);
	
	 
	//decl k=1; // If k=1 then log-logistic distribution.
	//decl c=1;	// If c=1, then Pareto Type II.
	
	C0 = (((1-u1).^(-1/k) -1)).^(1/c);
	C1 = (((1-u2).^(-1/k) -1)).^(1/c);
	E0 = (((1-v1).^(-1/k) -1)).^(1/c);
	E1 = (((1-v2).^(-1/k) -1)).^(1/c)+theta;
	return C0~C1~E0~E1;

}

dgp3(N,theta,K,alpha)
{

	decl C0, C1, E0, E1;
	decl rho = 0.5;

	C0 = ranpareto(N,1,K,alpha);
	C1 = ranpareto(N,1,K,alpha);
	E0 = rho*rann(N,1) + sqrt(1-rho^2)* C0;
	E1 = rho*rann(N,1) + sqrt(1-rho^2)* C1+theta;
	return C0~C1~E0~E1;
}

iidbootstrap(C0,C1,E0,E1)
{

	decl C0star, C1star, E0star,E1star;
	decl N=rows(C0);
	decl rows0 = ranu(1,N);
	decl rows1 = ranu(1,N);
    C0star = C0[rows0*N][];
	E0star = E0[rows0*N][];
	C1star = C1[rows1*N][];
	E1star = E1[rows1*N][];
	return  C0star~C1star~ E0star~E1star;

}


wildbootstrap(C0,C1,E0,E1)
{

	decl C0star, C1star, E0star,E1star;
	decl N=rows(C0);
	decl wc0	= ranu(N,1) .> 0.5 .? 1 .: -1;
	decl we0	= ranu(N,1) .> 0.5 .? 1 .: -1;
	decl wc1	= ranu(N,1) .> 0.5 .? 1 .: -1;
	decl we1	= ranu(N,1) .> 0.5 .? 1 .: -1;	
	C0star = meanc(C0) + (C0 -quantilec(C0, 0.5) ).*wc0 ;
	E0star = meanc(E0) + (E0 -quantilec(E0, 0.5) ).*wc0 ;	// Note, you use wc0 in C_0 and E0 in the hope															// that will preserve correlation...
	C1star = meanc(C1) + (C1 -quantilec(C1, 0.5) ).*wc1 ;
	E1star = meanc(E1) + (E1 -quantilec(E1, 0.5) ).*wc1 ;
 	return  C0star~C1star~ E0star~E1star;
}

randinf(C0,C1,E0,E1)
{

	decl C0star, C1star, E0star,E1star;
	decl N=rows(C0);
	decl randomizer = rann(N*2,1) ;
	decl allC = (C0|C1)~randomizer;
	decl allE = (E0|E1)~randomizer;


	allC=sortbyc(allC, 1);
	allE=sortbyc(allE, 1);
	
	C0star =allC[0:N-1][0];
	E0star =allE[0:N-1][0];
	C1star =allC[N:][0];
	E1star =allE[N:][0];
	return  C0star~C1star~ E0star~E1star;
}

estimate(K, C_0, C_1, E_0, E_1)
{
	decl INB, vINB, CI_L, CI_U, stud;
	decl sqrtN= sqrt(rows(C_0));
	decl Cov_0 = variance(C_0~E_0);
	decl Cov_1 = variance(C_1~E_1);
	decl mC_0 = meanc(C_0);     decl vC_0 = Cov_0[0][0];
	decl mE_0 = meanc(E_0);     decl vE_0 = Cov_0[1][1];
	decl cov_0 = Cov_0[0][1];
	
	decl mC_1 = meanc(C_1);     decl vC_1 = Cov_1[0][0];
	decl mE_1 = meanc(E_1);     decl vE_1 = Cov_1[1][1];
	decl cov_1 = Cov_1[0][1];
	
	decl Delta_E = mE_1 - mE_0;
	decl Delta_C = mC_1 - mC_0;
	
	INB =  K*Delta_E - Delta_C;
	vINB = (K^2)*(vE_0 + vE_1) + (vC_0 + vC_1) - 2*K*(cov_0 + cov_1);
	
	stud = sqrtN*INB/sqrt(vINB);
	CI_L = INB - 1.96*sqrt(vINB)/sqrtN;
	CI_U = INB + 1.96*sqrt(vINB)/sqrtN;	

	return INB~vINB~ CI_L~ CI_U~ stud;
}


asyInf(theta, CI_L, CI_U, stud)
{

	decl ci_rej, p_rej;
	
	if( theta> CI_L && theta< CI_U){
		ci_rej= 1;
	}
	else if( theta< CI_L || theta > CI_U){
		ci_rej= 0;
	}
	// Two sided test. Arguably, we should be doing 1 sided tests (INB > 0?)
	if(fabs(stud)>1.96)
	{
		p_rej=1;
	}
	else{
		p_rej=0;
	}
	return ci_rej~p_rej;
}


bootFisher(theta, INB,vINB, B, K, type, C0, C1, E0, E1)
{
	/* theta = INB in bootstrap and the value under H0 in Fisher
	*/
	decl b, bootdata, C0s, C1s, E0s, E1s;
	decl estimates, INBb, vINBb, z, mINB, mSeINB, mZ;
	decl ci_l, ci_u, ci_rej, mRej,p_rej;

	mINB=<>;
	mSeINB=<>;
	mZ=<>;
	for(b = 0; b < B; ++b)
	{

	
		bootdata = randinf(C0,C1,E0,E1);

		C0s=bootdata[][0];
		C1s=bootdata[][1];
		E0s=bootdata[][2];
		E1s=bootdata[][3];
		
		estimates = estimate(K, C0s, C1s, E0s, E1s);	 //INB~vINB~ CI_L~ CI_U~ stud;

		INBb = estimates[0][0];
		vINBb = estimates[0][1];
		z = (INBb)/sqrt(vINBb);	  
		
		mINB=mINB|INBb;
		mSeINB=mSeINB|sqrt(vINBb);
		mZ = mZ| z ;		  
	}
	
	/* Get percentiles of the boostrap distribution for the studentized INB */
	decl c_alpha = (B+1)*0.05/2;
	decl c_beta = (B+1)*(1-(0.05/2)); 
	mZ = sortc(mZ);
	decl z_u = mZ[c_alpha -1]; // 1-a percentile goes to lower limit
	decl z_l = mZ[c_beta -1];  // a percentile goes to upper limit

	ci_l = INB - z_l*sqrt(vINB);
	ci_u = INB - z_u*sqrt(vINB);
	if( theta> ci_l && theta< ci_u){
		ci_rej= 1;
	}
	else if( theta< ci_l || theta > ci_u){
		ci_rej =  0;
	}
	
	mRej = meanc( fabs(mZ).>= fabs( INB/sqrt(vINB)) .? 1 .: 0 );
	p_rej= mRej <= 0.05.? 1 .: 0;

	return ci_rej~p_rej;
 
}


bootInf(theta, INB,vINB, B, K, type, C0, C1, E0, E1)
{
	/* theta = INB in bootstrap and the value under H0 in Fisher
	*/
	decl b, bootdata, C0s, C1s, E0s, E1s;
	decl estimates, INBb, vINBb, z, mINB, mSeINB, mZ;
	decl ci_l, ci_u, ci_rej, mRej,p_rej;

	mINB=<>;
	mSeINB=<>;
	mZ=<>;
	for(b = 0; b < B; ++b)
	{

		if(type ==0){
			bootdata = iidbootstrap(C0,C1,E0,E1);
		}
		
		else {
			bootdata = wildbootstrap(C0,C1,E0,E1);
		}

		C0s=bootdata[][0];
		C1s=bootdata[][1];
		E0s=bootdata[][2];
		E1s=bootdata[][3];
		
		estimates = estimate(K, C0s, C1s, E0s, E1s);	 //INB~vINB~ CI_L~ CI_U~ stud;

		INBb = estimates[0][0];
		vINBb = estimates[0][1];
		z = (INBb-INB)/sqrt(vINBb);	  
		
		mINB=mINB|INBb;
		mSeINB=mSeINB|sqrt(vINBb);
		mZ = mZ| z ;		  
	}
	
	/* Get percentiles of the boostrap distribution for the studentized INB */
	decl c_alpha = (B+1)*0.05/2;
	decl c_beta = (B+1)*(1-(0.05/2)); 
	mZ = sortc(mZ);
	decl z_u = mZ[c_alpha -1]; // 1-a percentile goes to lower limit
	decl z_l = mZ[c_beta -1];  // a percentile goes to upper limit

	ci_l = INB - z_l*sqrt(vINB);
	ci_u = INB - z_u*sqrt(vINB);
	if( theta> ci_l && theta< ci_u){
		ci_rej= 1;
	}
	else if( theta< ci_l || theta > ci_u){
		ci_rej =  0;
	}
	
	mRej = meanc( fabs(mZ).>= fabs( INB/sqrt(vINB)) .? 1 .: 0 );
	p_rej= mRej <= 0.05.? 1 .: 0;

	return ci_rej~p_rej;
 
}

iteration(K,N,theta,k,c,B, R, dgp)
{
	decl data;
	if(dgp==1){
		data = dgp1(N, theta, k, c);		//return C0~C1~E0~E1;
	}
	else if(dgp==2) {
		data = dgp2(N, theta,k,c);
	}
	else if(dgp==3) {
		data = dgp3(N, theta,k,c);
	}
	
	//

	decl C0 = data[][0];
	decl C1 = data[][1];
	decl E0 = data[][2];
	decl E1 = data[][3];
	
	decl estimates = estimate(K, C0, C1, E0, E1);
	decl INB = estimates[0][0];
	decl vINB = estimates[0][1];
	decl CI_L = estimates[0][2];
	decl CI_U = estimates[0][3];
	decl stud = estimates[0][4];;
	
	decl asyresults = asyInf(theta, CI_L, CI_U, stud);
	
	decl bootresults1 = bootInf(theta, INB,vINB, B, K, 0, C0, C1, E0, E1);

	decl bootresults2 = bootInf(theta, INB,vINB, B, K, 1, C0, C1, E0, E1);

	decl junk = bootFisher(theta, INB,vINB, B, K, 1, C0, C1, E0, E1);
	
	return(asyresults~bootresults1~bootresults2~junk);
}
main()
{

decl K = 1;
decl N=100;
decl theta = 0;
decl k = 2;
decl c= 0.5;
decl B= 199;
decl R= 10000;
decl dgp = 1; // 1= half normal; 2 = Singh Maddala; 3 = Pareto.
decl mRej=<>;


decl i, iter;
for(i=0; i<R; ++i){

	iter = iteration(K,N,theta,k,c,B, R, dgp) ;
	
	mRej=mRej|iter;
}

print(meanc(mRej));

}


