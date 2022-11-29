#include <oxstd.h>
#include <oxprob.h>
#include <oxfloat.h>
#include <oxdraw.oxh>


dgp2(N, C_0, C_1, E_0, E_1)
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
	
	 
	decl k=1; // If k=1 then log-logistic distribution.
	decl c=1;	// If c=1, then Pareto Type II.
	
	C_0[0] = (((1-u1).^(-1/k) -1)).^(1/c);
	C_1[0] = (((1-u2).^(-1/k) -1)).^(1/c);
	E_0[0] = (((1-v1).^(-1/k) -1)).^(1/c);
	E_1[0] = (((1-v2).^(-1/k) -1)).^(1/c);

}

dgp4(N,K,alpha, C_0, C_1, E_0, E_1)
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

main()
{

decl N;
decl b,r, C_0, C_1, E_0, E_1, INB, vINB, CI_L, CI_U, stud;
decl C_0_star, E_0_star, C_1_star, E_1_star;
decl INB_b, vINB_b, CI_L_b, CI_U_b, stud_b, mRej;

/* Parameters of the simulation */

decl K = 1 ; //20000;
decl B = 199;
decl R = 10000 ;
decl theta = 0;

for(N = 20; N<1000; N*=2){

	decl c_alpha = (B+1)*0.05/2;
	decl c_beta = (B+1)*(1-(0.05/2)); 

	decl mB_esti=<>;
	decl mR_asy = <>;
	decl mR_bot = <>;
	decl mP_asy=<>;
	decl mP_bot=<>;

	for(r=0; r<R; ++r)
	{
	   	
		//dgp2(N, &C_0, &C_1, &E_0, &E_1);
		dgp4(N,1,1.5, &C_0, &C_1, &E_0, &E_1);	  //First param is minimum and second is alpha (more tail)

		if(r==0){
		println(minc(C_0)~maxc(C_0)~meanc(C_0)~minc(C_1)~maxc(C_1)~meanc(C_1));
		}
		// Asymptotic Method.
		
		estimate(K, C_0, C_1, E_0, E_1, &INB, &vINB, &CI_L, &CI_U, &stud);
		
		if( theta> CI_L && theta< CI_U){
		  	mR_asy = mR_asy| 1;
		}
		else if( theta< CI_L || theta > CI_U){
			mR_asy = mR_asy| 0;
		}
		// Two sided test. Arguably, we should be doing 1 sided tests (INB > 0?)
		if(fabs(stud)>1.96)
		{
			mP_asy=mP_asy|1;
		}
		else{
			mP_asy=mP_asy|0;
		}
			
		/* Part 3: IID Bootstrap */

		decl mINB = <>;
		decl mSeINB = <>;
		decl mZ = <>;
		for(b = 0; b < B; ++b)
		{
			randinf(N,  C_0~E_0, C_1~E_1, &C_0_star, &E_0_star, &C_1_star, &E_1_star);
			estimate(K, C_0_star, C_1_star, E_0_star, E_1_star, &INB_b, &vINB_b, &CI_L_b, &CI_U_b,&stud_b);
			mINB = mINB | INB_b;
			//mSeINB = mSeINB | sqrt(vINB_b);
			mZ = mZ | (INB_b / sqrt(vINB_b));	  
		}

		/* Get percentiles of the boostrap distribution for the studentized INB */
		mZ = sortc(mZ);
	   	decl z_u = mZ[c_alpha -1]; // 1-a percentile goes to lower limit
		decl z_l = mZ[c_beta -1];  // a percentile goes to upper limit

		// Rejections null hypothesis.
		mRej = meanc( fabs(mZ).>= fabs( INB/sqrt(vINB)) .? 1 .: 0 );
		mP_bot=mP_bot | (mRej <= 0.05.? 1 .: 0);

		// Studentized.	  Comment what follows, to line 281 to run the non-studentized version. 

		mB_esti = mB_esti | meanc(mINB); 
		decl ci_l_b	= INB - z_l*sqrt(vINB);
		decl ci_u_b = INB - z_u*sqrt(vINB);

		// Coverage of the boostrap CI
		
		if( theta> ci_l_b && theta< ci_u_b){
			mR_bot= mR_bot| 1;
		}
		else if( theta< ci_l_b || theta > ci_u_b){
			mR_bot = mR_bot| 0;
		}

	    

			

	}		

	

	print("%c",{"N", "R", "B", "Asymptotic","Bootstrap","R ASY", "R b"}, "%r",{"Coverage"}, N~R~B~meanc(mR_asy)~meanc(mR_bot)~meanc(mP_asy)~meanc(mP_bot)) ;
  



}



}





/* Semiparametric bootstrap

alphahat(N,Y, alphaHat,  y0, nDraws)
{
	decl k=sqrt(N);
	decl n=N-1;	
	decl logYnk1 = Y[n-k+1][];
	decl H = meanc(log(Y[n-k+2:][])-log(logYnk1));
	alphaHat[0]=1/H;
	//H = 0.95;// THis is either H or a value in 0-1
	decl ptail = H*k/N;								  
	y0[0] = Y[floor(N*(1-ptail[0]))][];
 	nDraws[0]= floor((1-ptail[0])*N);
}

spbootstrap(N,	C_0, C_0_aph, C_0_y0, C_0_nDraws,
				C_1, C_1_aph, C_1_y0, C_1_nDraws,
				E_0, E_0_aph, E_0_y0, E_0_nDraws,
				E_1, E_1_aph, E_1_y0, E_1_nDraws,
				C_0_star, E_0_star, C_1_star, E_1_star)
{
	
	decl ytrim = C_0[0:C_0_nDraws-1][];
	
	decl ystar = ytrim[C_0_nDraws*ranu(1,C_0_nDraws)][];
	
	decl tailstar= C_0_y0*( (1- ranu(N-C_0_nDraws,1)).^(-1/C_0_aph) );

	C_0_star[0]=	(ystar|tailstar);

	ytrim = C_1[0:C_1_nDraws-1][];
	ystar = ytrim[C_1_nDraws*ranu(1,C_1_nDraws)][];
	tailstar= C_1_y0*( (1- ranu(N-C_1_nDraws,1)).^(-1/C_1_aph) );
	C_1_star[0]=	(ystar|tailstar);

	ytrim = E_0[0:E_0_nDraws-1][];
	ystar = ytrim[E_0_nDraws*ranu(1,E_0_nDraws)][];
	tailstar= E_0_y0*( (1- ranu(N-E_0_nDraws,1)).^(-1/E_0_aph) );
	E_0_star[0]=	(ystar|tailstar);
	
	ytrim = E_1[0:E_1_nDraws-1][];
	ystar = ytrim[E_1_nDraws*ranu(1,E_1_nDraws)][];
	tailstar= E_1_y0*( (1- ranu(N-E_1_nDraws,1)).^(-1/E_1_aph) );
	E_1_star[0]=	(ystar|tailstar);
    
}

/*alphahat(N, Y)
{
	
	decl k=sqrt(N);
	decl n=N-1;
	
	decl sortY = sortc(Y);

	
	decl logYnk1 = sortY[n-k+1][];
	decl H = meanc(log(sortY[n-k+2:][])-log(logYnk1));
	decl alphaHat=1/H;
	//H = 0.5;// THis is either H or a value in 0-1
	decl ptail = H*k/N;
									  
	decl y0 = sortY[floor(N*(1-ptail))][];
	
 	decl noDraws= floor((1-ptail)*N);
	decl ytrim = sortY[0:noDraws-1][];
	decl ystar = ytrim[noDraws*ranu(1,noDraws)][];
	decl tailstar= y0*( (1- ranu(N-noDraws,1)).^(-1/alphaHat) );

	return(ystar|tailstar);

	
}

	 */

//

*/
