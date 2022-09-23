#include <oxstd.h>

#include <oxprob.h>

#include <oxfloat.h>



// THis is the code used in the paper for Health Economics.





dgp(N,K, C_0, C_1, E_0, E_1)

{


decl sd_e = 2;
decl sd_c =0.05;   //0.05 => skew 0.69; 0.1 => skew 1.007; 0.5 => skew 2.93

decl mu_c0 =log(50000);
decl mu_c1 =log(50000 + 30000);

decl mu_e0 = 5;
decl mu_e1 = mu_e0 +1;



decl logC_0 = sd_c*rann(N,1)+mu_c0;
C_0[0] = exp(logC_0);
decl logC_1 = sd_c*rann(N,1)+mu_c1;
C_1[0] = exp(logC_1);

E_0[0] = sd_e*rann(N,1) + mu_e0 + 0.3*(logC_0 - mu_c0);
E_1[0] = sd_e*rann(N,1) + mu_e1 + 0.3*(logC_1 - mu_c1);



}







estimate(K, C_0, C_1, E_0, E_1, INB, vINB, CI_L, CI_U)

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




main()

{



/* Part 1: DGP. Taken from Grieve's paper. We don't simulate exactly the same data as we

    produce relatively different patterns of skewness*/

decl N;

for(N = 10; N<100; N*=2){

	   

	decl K = 20000;

	decl B = 199;

	decl R = 100000 ;

	decl theta = -10000;

	decl c_alpha = (B+1)*0.05/2;

	decl c_beta = (B+1)*(1-(0.05/2)); print(c_alpha~c_beta);

	

	

	decl b,r;

	decl C_0, C_1, E_0, E_1, INB, vINB, CI_L, CI_U;

	decl C_0_star, E_0_star, C_1_star, E_1_star;

	decl INB_b, vINB_b, CI_L_b, CI_U_b;

	

	

	decl mB_esti=<>;

	decl mR_asy = <>;

	decl mR_bot = <>;

	for(r=0; r<R; ++r)

	{

			

			/* Part 1: Call DGP*/

			

			dgp(N,K, &C_0, &C_1, &E_0, &E_1);

			

			

			/* Part 2: Estimates*/

			

			estimate(K, C_0, C_1, E_0, E_1, &INB, &vINB, &CI_L, &CI_U);

			

			if( theta> CI_L && theta< CI_U){

	

				mR_asy = mR_asy| 1;

			}

			else if( theta< CI_L || theta > CI_U){

	

				mR_asy = mR_asy| 0;

			}

			

			/* Part 3: IID Bootstrap */

			

			decl mINB = <>;

			decl mSeINB = <>;

			decl mZ = <>;

			for(b = 0; b < B; ++b)

			{

			    

			    iidbootstrap(N, C_0~E_0, C_1~E_1, &C_0_star, &E_0_star, &C_1_star, &E_1_star);

			    

			    estimate(K, C_0_star, C_1_star, E_0_star, E_1_star, &INB_b, &vINB_b, &CI_L_b, &CI_U_b);

			    mINB = mINB | INB_b;

			    mSeINB = mSeINB | sqrt(vINB_b);

			    mZ = mZ | ((INB_b - INB) /  sqrt(vINB_b));	  

			 

			}





			// Studentized.	  Comment what follows, to line 281 to run the non-studentized version. 



			mB_esti = mB_esti | meanc(mINB); 

			mZ = sortc(mZ);

			   

			decl z_u = mZ[c_alpha -1]; // 1-a percentile goes to lower limit

			decl z_l = mZ[c_beta -1];  // a percentile goes to upper limit

	

			decl ci_l_b	= INB - z_l*sqrt(vINB);

			decl ci_u_b = INB - z_u*sqrt(vINB);

	

			if( theta> ci_l_b && theta< ci_u_b){

	

				mR_bot= mR_bot| 1;

			}

			else if( theta< ci_l_b || theta > ci_u_b){

	

				mR_bot = mR_bot| 0;

			}

			

			 

			

		// noStudentized.	Remove flags to activate this code and run the non-studentized version. 



			/*

			mB_esti = mB_esti | meanc(mINB); 

			mZ = sortc(mINB);

			   

			decl z_l = mZ[c_alpha -1]; // 1-a percentile goes to lower limit

			decl z_u = mZ[c_beta -1];  // a percentile goes to upper limit

	

			decl ci_l_b	= z_l;

			decl ci_u_b = z_u;

	

			if( theta> ci_l_b && theta< ci_u_b){

	

				mR_bot= mR_bot| 1;

			}

			else if( theta< ci_l_b || theta > ci_u_b){

	

				mR_bot = mR_bot| 0;

			}

			

			  */

			

	}		

	

	print("%c",{"N", "R", "B", "Asymptotic","Bootstrap"}, "%r",{"Coverage"}, N~R~B~meanc(mR_asy)~meanc(mR_bot)) ;
  



}



}
