#include <oxstd.h>
#include <oxprob.h>
#include <oxfloat.h>
#include <oxdraw.oxh>

dgp1(N, theta, C_0, C_1, E_0, E_1)
{

/* Singh Maddala Distribution*/
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
	
	C_0[0] = fabs(quann(u1));
	C_1[0] = fabs(quann(u2));
	E_0[0] = fabs(quann(v1));
	E_1[0] = fabs(quann(v2))+theta;

}

dgp2(N, theta, k, c, C_0, C_1, E_0, E_1)
{

/* Singh Maddala Distribution*/
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
	
	C_0[0] = (((1-u1).^(-1/k) -1)).^(1/c);
	C_1[0] = (((1-u2).^(-1/k) -1)).^(1/c);
	E_0[0] = (((1-v1).^(-1/k) -1)).^(1/c);
	E_1[0] = (((1-v2).^(-1/k) -1)).^(1/c)+theta;

}

dgp4(N,theta,K,alpha, C_0, C_1, E_0, E_1)
{
decl rho = 0.5;

	C_0[0] = ranpareto(N,1,K,alpha);
	C_1[0] = ranpareto(N,1,K,alpha);
	E_0[0] = rho*rann(N,1) + sqrt(1-rho^2)* C_0[0];
	E_1[0] = rho*rann(N,1) + sqrt(1-rho^2)* C_1[0];
}

estimate(K, C_0, C_1, E_0, E_1, INB, vINB, CI_L, CI_U, stud)
{

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
	
	INB[0] = K*Delta_E - Delta_C;
	vINB[0] = (K^2)*(vE_0 + vE_1) + (vC_0 + vC_1) - 2*K*(cov_0 + cov_1);
	
	stud[0] = sqrt(rows(C_0))*INB[0]/sqrt(vINB[0]);
	CI_L[0] = INB[0] - 1.96*sqrt(vINB[0])/sqrt(rows(C_0));
	CI_U[0] = INB[0] + 1.96*sqrt(vINB[0])/sqrt(rows(C_0));	
	
}

iidbootstrap(N,B0, B1, C_0, E_0, C_1, E_1)
{
    
    decl B0star = B0[ranu(1,N) * N][];
    decl B1star = B1[ranu(1,N) * N][];
    C_0[0] = B0star[][0];
    E_0[0] = B0star[][1];

    C_1[0] = B1star[][0];
    E_1[0] = B1star[][1];
    
}



moonbootstrap(N,m,B0, B1, C_0, E_0, C_1, E_1)
{

	/* If N = m, we have the standard i.i.d. bootstrap */    

    decl B0star = B0[ranu(1,m) * N][];
    decl B1star = B1[ranu(1,m) * N][];
										 
    C_0[0] = B0star[][0];
    E_0[0] = B0star[][1];

    C_1[0] = B1star[][0];
    E_1[0] = B1star[][1];

}

wildbootstrap(N,B0, B1, C_0, E_0, C_1, E_1){

	decl c0	= B0[][0];
	decl e0	= B0[][1];
	decl c1	= B1[][0];
	decl e1	= B1[][1];

	decl wc0	= ranu(N,1) .> 0.5 .? 1 .: -1;
	decl we0	= ranu(N,1) .> 0.5 .? 1 .: -1;
	decl wc1	= ranu(N,1) .> 0.5 .? 1 .: -1;
	decl we1	= ranu(N,1) .> 0.5 .? 1 .: -1;
	
	C_0[0] = meanc(c0) + (c0 -quantilec(c0, 0.5) ).*wc0 ;
	E_0[0] = meanc(e0) + (e0 -quantilec(e0, 0.5) ).*wc0 ;	// Note, you use wc0 in C_0 and E0 in the hope															// that will preserve correlation...
	C_1[0] = meanc(c1) + (c1 -quantilec(c1, 0.5) ).*wc1 ;
	E_1[0] = meanc(e1) + (e1 -quantilec(e1, 0.5) ).*wc1 ;

}

randinf(N,B0, B1, C_0, E_0, C_1, E_1)
{

	decl randomizer = rann(N*2,1) ;
	
	decl B = (B0|B1)~randomizer;

	B=sortbyc(B, 2);
	
	C_0[0] =B[0:N-1][0];
	E_0[0] =B[0:N-1][1];
	C_1[0] =B[N:][0];
	E_1[0] =B[N:][1];
	
}

asymptotic(K, theta, C_0, C_1, E_0, E_1, inb, vinb, ci_l, ci_u, tstat, ci_rej, p_rej)
{

	decl  INB, vINB, CI_L, CI_U, stud;
	
	estimate(K, C_0, C_1, E_0, E_1, &INB, &vINB, &CI_L, &CI_U, &stud);
	inb[0]=INB;
	vinb[0]=vINB;
	ci_l[0]=CI_L;
	ci_u[0]=CI_U;
	tstat[0]=stud;
	
	if( theta> CI_L && theta< CI_U){
		ci_rej[0]= 1;
	}
	else if( theta< CI_L || theta > CI_U){
		ci_rej[0]= 0;
	}
	// Two sided test. Arguably, we should be doing 1 sided tests (INB > 0?)
	if(fabs(stud)>1.96)
	{
		p_rej[0]=1;
	}
	else{
		p_rej[0]=0;
	}
}

bs(K, B, type, mOn, N,theta, INB, vINB, C_0, C_1, E_0, E_1, ci_l, ci_u, ci_rej, p_rej)
{
	decl C_0_star, E_0_star, C_1_star, E_1_star,b;
	decl INB_b, vINB_b, CI_L_b, CI_U_b, stud_b, mRej;
	decl c_alpha = (B+1)*0.05/2;
	decl c_beta = (B+1)*(1-(0.05/2)); 
	decl mINB = <>;
	decl mSeINB = <>;
	decl mZ = <>;
	for(b = 0; b < B; ++b)
	{

		if(type ==0){
			iidbootstrap(N, C_0~E_0, C_1~E_1, &C_0_star, &E_0_star, &C_1_star, &E_1_star);
		}
		else if(type ==1){
			moonbootstrap(N, mOn, C_0~E_0, C_1~E_1, &C_0_star, &E_0_star, &C_1_star, &E_1_star);
		}
		else if (type==2){
			wildbootstrap(N,  C_0~E_0, C_1~E_1, &C_0_star, &E_0_star, &C_1_star, &E_1_star);
		}
		estimate(K, C_0_star, C_1_star, E_0_star, E_1_star, &INB_b, &vINB_b, &CI_L_b, &CI_U_b,&stud_b);
		mINB = mINB | INB_b;
		mSeINB = mSeINB | sqrt(vINB_b);
		mZ = mZ | ((INB_b - INB) /  sqrt(vINB_b));	  
	}

	/* Get percentiles of the boostrap distribution for the studentized INB */
	mZ = sortc(mZ);
	decl z_u = mZ[c_alpha -1]; // 1-a percentile goes to lower limit
	decl z_l = mZ[c_beta -1];  // a percentile goes to upper limit

	// Rejections null hypothesis.
	mRej = meanc( fabs(mZ).>= fabs( INB/sqrt(vINB)) .? 1 .: 0 );
	p_rej[0]= mRej <= 0.05.? 1 .: 0;
	// coverage
	 
	ci_l[0]	= INB - z_l*sqrt(vINB);
	ci_u[0] = INB - z_u*sqrt(vINB);
	if( theta> ci_l[0] && theta< ci_u[0]){
		ci_rej[0]= 1;
	}
	else if( theta< ci_l[0] || theta > ci_u[0]){
		ci_rej[0] =  0;
	}

}

fisher(K, B, N,theta, INB, vINB, C_0, C_1, E_0, E_1, ci_l, ci_u, ci_rej, p_rej)
{
  	decl C_0_star, E_0_star, C_1_star, E_1_star,b;
	decl INB_b, vINB_b, CI_L_b, CI_U_b, stud_b, mRej;
	decl c_alpha = (B+1)*0.05/2;
	decl c_beta = (B+1)*(1-(0.05/2)); 
	decl mINB = <>;
	decl mSeINB = <>;
	decl mZ = <>;
		for(b = 0; b < B; ++b)
		{
			randinf(N,  C_0~E_0, C_1~E_1, &C_0_star, &E_0_star, &C_1_star, &E_1_star);
			estimate(K, C_0_star, C_1_star, E_0_star, E_1_star, &INB_b, &vINB_b, &CI_L_b, &CI_U_b,&stud_b);
			mINB = mINB | INB_b;
		 
			mZ = mZ | (INB_b / sqrt(vINB_b));	  
		}

		/* Get percentiles of the boostrap distribution for the studentized INB */
	mZ = sortc(mZ);
	decl z_u = mZ[c_alpha -1]; // 1-a percentile goes to lower limit
	decl z_l = mZ[c_beta -1];  // a percentile goes to upper limit

	// Rejections null hypothesis.
	mRej = meanc( fabs(mZ).>= fabs( INB/sqrt(vINB)) .? 1 .: 0 );
	p_rej[0]= mRej <= 0.05.? 1 .: 0;
	// coverage
	 
	ci_l[0]	= INB - z_l*sqrt(vINB);
	ci_u[0] = INB - z_u*sqrt(vINB);
	if( theta> ci_l[0] && theta< ci_u[0]){
		ci_rej[0]= 1;
	}
	else if( theta< ci_l[0] || theta > ci_u[0]){
		ci_rej[0] =  0;
	}


}

main()
{

decl N;
decl b,r, C_0, C_1, E_0, E_1, INB, vINB, CI_L, CI_U, stud;
decl C_0_star, E_0_star, C_1_star, E_1_star;
decl INB_b, vINB_b, CI_L_b, CI_U_b, stud_b, mRej;

/* Parameters of the simulation */

decl K = 1 ; //20000;
decl B = 399;
decl R = 50000 ;
decl theta = 0;
decl k=1;  // If k=1 then log-logistic distribution.
decl c=2;  // If c=1, then Pareto Type II.
decl typeBoot, mOn;

decl mRES =<>;
for(N = 50; N<=200; N*=2){

	decl c_alpha = (B+1)*0.05/2;
	decl c_beta = (B+1)*(1-(0.05/2)); 

	decl mB_esti=<>;
	decl mR_asy = <>;
	decl mR_bot = <>;
	decl mR_fisher = <>;
	decl mP_asy=<>;
	decl mP_bot=<>;
	decl mP_fisher=<>;
	decl mR_iidbot =<>;
	decl mP_iidbot=<>;
	decl mR_moon70bot =<>;
	decl mP_moon70bot =<>;
	decl mR_moon90bot =<>;
	decl mP_moon90bot =<>;
	decl mR_moon80bot =<>;
	decl mP_moon80bot =<>;
	
	for(r=0; r<R; ++r)
	{

	//	dgp1(N, theta,&C_0, &C_1, &E_0, &E_1);

		dgp2(N, theta,k,c, &C_0, &C_1, &E_0, &E_1);
		//dgp4(N,theta,1,1 , &C_0, &C_1, &E_0, &E_1);	  //First param is minimum and second is alpha (more tail)

		
		decl ci_rej, p_rej;
		asymptotic(K, theta, C_0, C_1, E_0, E_1, &INB, &vINB, &CI_L, &CI_U, &stud, &ci_rej, &p_rej);
		mR_asy =mR_asy| ci_rej;
		mP_asy= mP_asy| p_rej;
		

		decl ci_rej_b, p_rej_b, ci_l_b, ci_u_b;
		mOn=0;
		typeBoot =0;
		bs(K, B, typeBoot, mOn, N,theta, INB, vINB, C_0, C_1, E_0, E_1, &ci_l_b, &ci_u_b, &ci_rej_b, &p_rej_b);
		mR_iidbot= mR_iidbot| ci_rej_b;
		mP_iidbot= mP_iidbot | p_rej_b;

		
		mOn=floor(0.9*N);
		typeBoot =1;
		bs(K, B, typeBoot, mOn, N,theta, INB, vINB, C_0, C_1, E_0, E_1, &ci_l_b, &ci_u_b, &ci_rej_b, &p_rej_b);
		mR_moon90bot= mR_moon90bot| ci_rej_b;
		mP_moon90bot= mP_moon90bot | p_rej_b;

		mOn=floor(0.8*N);
		typeBoot =1;
		bs(K, B, typeBoot, mOn, N,theta, INB, vINB, C_0, C_1, E_0, E_1, &ci_l_b, &ci_u_b, &ci_rej_b, &p_rej_b);
		mR_moon80bot= mR_moon80bot| ci_rej_b;
		mP_moon80bot= mP_moon80bot | p_rej_b;

		mOn=floor(0.7*N);
		typeBoot =1;
		bs(K, B, typeBoot, mOn, N,theta, INB, vINB, C_0, C_1, E_0, E_1, &ci_l_b, &ci_u_b, &ci_rej_b, &p_rej_b);
		mR_moon70bot= mR_moon70bot| ci_rej_b;
		mP_moon70bot= mP_moon70bot | p_rej_b;
		
		 
		typeBoot =2;
		bs(K, B, typeBoot, mOn, N,theta, INB, vINB, C_0, C_1, E_0, E_1, &ci_l_b, &ci_u_b, &ci_rej_b, &p_rej_b);
		mR_bot= mR_bot| ci_rej_b;
		mP_bot= mP_bot | p_rej_b;




		decl ci_rej_f, p_rej_f, ci_l_f, ci_u_f;
		fisher(K, B, N,theta, INB, vINB, C_0, C_1, E_0, E_1, &ci_l_f, &ci_u_f, &ci_rej_f, &p_rej_f);
		mR_fisher= mR_fisher| ci_rej_f;
		mP_fisher= mP_fisher | p_rej_f;

	}		

	print("%c",{"N", "R", "B", "Theta","k","c"},  N~R~B~theta~k~c);
	print("%c",{"Coverage", "Rejections"}, "%r",{"Asy", "iid", "moon 90", "moon 80", "moon 70", "Wild","Fisher"}, (meanc(mR_asy)~meanc(mP_asy))|
	(meanc(mR_iidbot)~meanc(mP_iidbot))|(meanc(mR_moon90bot)~meanc(mP_moon90bot))|
	(meanc(mR_moon80bot)~meanc(mP_moon80bot))|(meanc(mR_moon80bot)~meanc(mP_moon70bot)) |
	(meanc(mR_bot)~meanc(mP_bot))|(meanc(mR_fisher)~meanc(mP_fisher)));
  
	mRES= mRES|	N~R~B~ (
		(meanc(mR_asy)~meanc(mP_asy))~
	(meanc(mR_iidbot)~meanc(mP_iidbot))~(meanc(mR_moon90bot)~meanc(mP_moon90bot))~
	(meanc(mR_moon80bot)~meanc(mP_moon80bot))~(meanc(mR_moon80bot)~meanc(mP_moon70bot)) ~
	(meanc(mR_bot)~meanc(mP_bot))~(meanc(mR_fisher)~meanc(mP_fisher))
	);
}

print(mRES);

decl s=sprint("mRES","_k=",k,"_c=",c,"_theta=",theta,"_dgp=2");
print(s);
savemat(sprint(s,".csv"), mRES, {"N", "R", "B","CAsyy","Rasy", "Ciid","Riid", "Cmoon90", "Rmoon90",
"Cmoon80","Rmoon80", "Cmoon70","Rmoon70", "CWild","RWild","CFisher","RFisher"});
}

				 