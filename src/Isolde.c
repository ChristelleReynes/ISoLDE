// Version 16/03/2016
// Conditional Compilation
// Version 20/10/2015
// Enhanced Memory Allocation Error Resolution -> NBoot[0] = 0 if error

#include <R.h>
#include  <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>
#include <Rdefines.h>
#ifdef _OPENMP
    #include <omp.h>
#endif  
#include <R_ext/Rdynload.h>

void Stat9NonParF(int vline, double *vec, int *iM1, int *iM2, int *iP1, int *iP2, int nM1, int nM2, double *MSorted, double *PSorted, int *M12, int *P12, double *denom, double *diffp, double *crit);

SEXP IsoldeP1(SEXP R_NCore, SEXP R_NbReadNormFil, SEXP R_NRowReadNormFil, SEXP R_IndM, SEXP R_IndP, SEXP R_crit9glob, SEXP R_diff9glob, SEXP R_denom9glob, SEXP R_distdiffpropH0,
 SEXP R_THistodistdiffpropH0, SEXP R_Tvecprob, SEXP R_Tdistdiffcumul, SEXP R_Nboot, SEXP R_enrcrit9boot, SEXP R_moylocH0, SEXP R_pvalboot, SEXP R_WithRepeat)

{
	double *NbReadNormFil = REAL(R_NbReadNormFil), *HDDPH0 = REAL(R_THistodistdiffpropH0), *vecprob = REAL(R_Tvecprob), *enrcrit9boot = REAL(R_enrcrit9boot);
	double *distdiffcumul = REAL(R_Tdistdiffcumul), *moylocH0 = REAL(R_moylocH0), *pvalboot = REAL(R_pvalboot); //, *enrcrit9bootZ = REAL(R_enrcrit9bootZ)
	double *crit9glob = REAL(R_crit9glob), *diff9glob = REAL(R_diff9glob), *denom9glob = REAL(R_denom9glob), *distdiffpropH0 = REAL(R_distdiffpropH0);
	int *nNRNF = INTEGER(R_NRowReadNormFil), *Nboot = INTEGER(R_Nboot), *NCore = INTEGER(R_NCore), *WithRepeat = INTEGER(R_WithRepeat); // *ChoixRep1 = INTEGER(R_TChoixRep1);

	int sNRNF = length(R_NbReadNormFil); int cNRNF = floor(sNRNF / nNRNF[0]);
	int i, j, k, line, col, NbreDDPH0, Nbre2DDPH0, NColDDPH0, NumClasse, NbreVecProb, numBoot; //m, 
	double *vecBRNF, *vecDDPH0, TempDDPH0, *vecEC9B; //, *vecEC9BZ
	double AleaSelect, nu0;   //  tempc, tempd, tempe, 
	
	int nIndM1 = length(VECTOR_ELT(R_IndM, 0)); int nIndM2 = length(VECTOR_ELT(R_IndM, 1)); int nIndM = nIndM1 + nIndM2; 
	int nIndP1 = nIndM1; int nIndP2 = nIndM2; int nIndP = nIndP1 + nIndP2; 
	int *IndM, *IndM1, *IndM2, *IndP, *IndP1, *IndP2, *ChoixRepM, *ChoixRepP, *NbBootOk, *IndMPv2;	

	int AllocOk = 1;
	IndM = malloc(nIndM*sizeof(int)); if (IndM == NULL) { Rprintf( "Not enough memory to allocate buffer IndM \n"); AllocOk = 0;}
	IndM1 = malloc(nIndM1*sizeof(int)); if (IndM1 == NULL) { Rprintf( "Not enough memory to allocate buffer IndM1 \n"); AllocOk = 0;}
	IndM2 = malloc(nIndM2*sizeof(int)); if (IndM2 ==NULL) { Rprintf( "Not enough memory to allocate buffer IndM2 \n"); AllocOk = 0;}
	IndP = malloc(nIndP*sizeof(int)); if (IndP == NULL) { Rprintf( "Not enough memory to allocate buffer IndP \n"); AllocOk = 0;}
	IndP1 = malloc(nIndP1*sizeof(int)); if (IndP1 == NULL) { Rprintf( "Not enough memory to allocate buffer IndP1 \n"); AllocOk = 0;}
	IndP2 = malloc(nIndP2*sizeof(int)); if (IndP2 == NULL) { Rprintf( "Not enough memory to allocate buffer IndP2 \n"); AllocOk = 0;}
	IndMPv2 = malloc(cNRNF*Nboot[0]*sizeof(int)); if (IndMPv2 == NULL) { Rprintf( "Not enough memory to allocate buffer IndMPv2 \n"); AllocOk = 0;}
	ChoixRepM = malloc(nIndM*Nboot[0]*sizeof(int)); if (ChoixRepM == NULL) { Rprintf( "Not enough memory to allocate buffer ChoixRepM \n"); AllocOk = 0;}
	ChoixRepP = malloc(nIndM*Nboot[0]*sizeof(int)); if (ChoixRepP == NULL) { Rprintf( "Not enough memory to allocate buffer ChoixRepP \n"); AllocOk = 0;}
	NbBootOk = malloc(nNRNF[0]*sizeof(int));  if (NbBootOk == NULL) { Rprintf( "Not enough memory to allocate buffer NbBootOk \n"); AllocOk = 0;}
	if  (AllocOk == 0)
		{Nboot[0] = 0; Rprintf( "ERROR - Not enough memory to allocate buffers P1 \n"); goto AbNormalEnd;}
	
	//Set the number of threads
	int HaveCores = 1;
    #ifdef _OPENMP
        HaveCores = omp_get_num_procs(); Rprintf("Open MP enabled \n");
    #else   
         Rprintf("Open MP disabled \n");
    #endif  
	Rprintf("Have Cores : %d \n", HaveCores);
	int UseCores, UseCoresNew;
	if (NCore[0]<=0) {UseCores = 1;}
	else
		{
			if (NCore[0]>=100) {UseCores = HaveCores;}
			else {UseCores = (int) floor(HaveCores*NCore[0]/100);}
		}
	if (UseCores<1) {UseCores = 1;}
	#ifdef _OPENMP
	omp_set_num_threads(UseCores);
	#endif
	Rprintf("UseCores : %d \n", UseCores);
	//int My_Thread;	
	
	
	GetRNGstate();
	
	for (i=0; i<nIndM1; i++)
	{
		IndM1[i] = INTEGER(VECTOR_ELT(R_IndM, 0))[i]; IndM[i] = IndM1[i]; IndP1[i] = INTEGER(VECTOR_ELT(R_IndP, 0))[i]; IndP[i] = IndP1[i];
	}
	for (i=0; i<nIndM2; i++)
	{
		IndM2[i] = INTEGER(VECTOR_ELT(R_IndM, 1))[i]; IndM[i+nIndM1] = IndM2[i]; IndP2[i] = INTEGER(VECTOR_ELT(R_IndP, 1))[i]; IndP[i+nIndP1] = IndP2[i];
	}

	NbreDDPH0 = 0; Nbre2DDPH0 = 0; NColDDPH0 = (int) choose(nIndM1, 2) + (int) choose(nIndM2, 2);
	
	for (i=0; i<100; i++) {HDDPH0[i] = 0;}
	

	UseCoresNew = UseCores;
	double *TestMSortBuff, *TestPSortBuff;
	int *TestM12Buff, *TestP12Buff;
	do
	{
		AllocOk = 1; 
		TestM12Buff = malloc(UseCoresNew*nIndM*sizeof(int)); if (TestM12Buff == NULL) { Rprintf( "Not enough memory to allocate buffer TestM12Buff \n"); AllocOk = 0;}
		TestP12Buff = malloc(UseCoresNew*nIndM*sizeof(int)); if (TestP12Buff == NULL) { Rprintf( "Not enough memory to allocate buffer TestP12Buff \n"); AllocOk = 0;}
		TestMSortBuff = malloc(UseCoresNew*nIndM*sizeof(double)); if (TestMSortBuff == NULL) { Rprintf( "Not enough memory to allocate buffer TestMSortBuff \n"); AllocOk = 0;}
		TestPSortBuff = malloc(UseCoresNew*nIndM*sizeof(double)); if (TestPSortBuff == NULL) { Rprintf( "Not enough memory to allocate buffer TestPSortBuff \n"); AllocOk = 0;}
		free(TestPSortBuff); free(TestMSortBuff); free(TestP12Buff); free(TestM12Buff); 
		if (AllocOk ==0) {UseCoresNew = UseCoresNew -1; Rprintf( "WARNING - Not enough memory to allocate parallel buffers - Lowering to %d threads \n", UseCoresNew);} //NBoot[0] = 0; 
	}
	while ((AllocOk==0) && (UseCoresNew > 0));
	if (UseCoresNew ==0) {Nboot[0] = 0; Rprintf( "ERROR - Not enough memory to allocate parallel buffers \n"); goto AbNormalEnd2;}
	#ifdef _OPENMP
        omp_set_num_threads(UseCoresNew); 
    #endif

	#pragma omp parallel for \
	private(line, vecBRNF, vecDDPH0, col, j, k, TempDDPH0, NumClasse)\
 	schedule(static) reduction(+:NbreDDPH0)
	for (line = 0; line<nNRNF[0]; line++)
	{
		double *MSortBuff, *PSortBuff, intpart, tempa;
		int *M12Buff, *P12Buff;
		M12Buff = malloc(nIndM*sizeof(int)); if (M12Buff == NULL) { Rprintf( "Not enough memory to allocate buffer M12Buff \n");}
		P12Buff = malloc(nIndM*sizeof(int)); if (P12Buff == NULL) { Rprintf( "Not enough memory to allocate buffer P12Buff \n");}
		MSortBuff = malloc(nIndM*sizeof(double)); if (MSortBuff == NULL) { Rprintf( "Not enough memory to allocate buffer MSortBuff \n");}
		PSortBuff = malloc(nIndM*sizeof(double)); if (PSortBuff == NULL) { Rprintf( "Not enough memory to allocate buffer PSortBuff \n");}
		
		Stat9NonParF(1, NbReadNormFil + line*cNRNF, IndM1, IndM2, IndP1, IndP2, nIndM1, nIndM2, MSortBuff, PSortBuff, M12Buff, 
			P12Buff, denom9glob + line, diff9glob + line, crit9glob + line);
		
		vecBRNF = NbReadNormFil + line*cNRNF; vecDDPH0 = distdiffpropH0 + line*NColDDPH0; // My_Thread = omp_get_thread_num();
		col = 0;
		for (j=0; j<(nIndM1-1); j++)
		{
			for (k=(j+1); k<nIndM1; k++)
			{
				vecDDPH0[col]= 0.5 + vecBRNF[IndM[j]-1]/(vecBRNF[IndM[j]-1]+vecBRNF[IndP[j]-1])-vecBRNF[IndM[k]-1]/(vecBRNF[IndM[k]-1]+vecBRNF[IndP[k]-1]);
				if ((col) > (NColDDPH0-1)) { Rprintf( "probleme vecDDPH0 : j = %d k = %d col = %d NColDDPH0 = %d \n", j, k, col, NColDDPH0);}
				TempDDPH0 = vecDDPH0[col];
				if (TempDDPH0>=0 && TempDDPH0<=1)
				{
					if (TempDDPH0 == 0)
						NumClasse = 0;
					else
					{
						tempa = TempDDPH0/0.01;
						if (modf(tempa, &intpart) == 0)
							NumClasse = (int) intpart-1;
						else	
							NumClasse = (int) intpart;
					}		
					#pragma omp critical
					{HDDPH0[NumClasse] = HDDPH0[NumClasse] + 1; }
					NbreDDPH0 = NbreDDPH0 + 1;
				}
				col = col + 1;
			}
		}
		
		for (j=0; j<(nIndM2-1); j++)
		{
			for (k=(j+1); k<nIndM2; k++)
			{
				vecDDPH0[col]= 0.5 + vecBRNF[IndM[j+nIndM1]-1]/(vecBRNF[IndM[j+nIndM1]-1]+vecBRNF[IndP[j+nIndM1]-1])-vecBRNF[IndM[k+nIndM1]-1]/(vecBRNF[IndM[k+nIndM1]-1]+vecBRNF[IndP[k+nIndM1]-1]);
				if ((col) > (NColDDPH0-1)) { Rprintf( "probleme vecDDPH0 : j = %d k = %d col = %d NColDDPH0 = %d \n", j, k, col, NColDDPH0);}
				TempDDPH0 = vecDDPH0[col];
				if (TempDDPH0>=0 && TempDDPH0<=1)
				{
					if (TempDDPH0 == 0)
						NumClasse = 0;
					else
					{
						tempa = TempDDPH0/0.01;
						if (modf(tempa, &intpart) == 0)
							NumClasse = (int) intpart-1;
						else	
							NumClasse = (int) intpart;
					}		
					#pragma omp critical
					{HDDPH0[NumClasse] = HDDPH0[NumClasse] + 1; }
					NbreDDPH0 = NbreDDPH0 + 1;
				}	
				col = col + 1;
			}
		}
		
		NbBootOk[line] = 0; 
		free(M12Buff); free(P12Buff); free(PSortBuff); free(MSortBuff);
	}
	#ifdef _OPENMP
        omp_set_num_threads(UseCores); 
    #endif  
	
	NbreVecProb = 0; distdiffcumul[0] = 0; 
	for (NumClasse = 0; NumClasse<100; NumClasse++)
	{
		if (HDDPH0[NumClasse]>0)
		{
			Nbre2DDPH0 = Nbre2DDPH0 + HDDPH0[NumClasse]; 
			HDDPH0[NbreVecProb] = HDDPH0[NumClasse] / NbreDDPH0;
			vecprob[NbreVecProb] = 0.005 + NumClasse*0.01; 
			distdiffcumul[NbreVecProb +1] = distdiffcumul[NbreVecProb] + HDDPH0[NbreVecProb];
			NbreVecProb = NbreVecProb + 1;
		}	
	}

	#pragma omp parallel for \
	private(numBoot, i, j)
	for (numBoot=0; numBoot<Nboot[0]; numBoot++)
	{
		for (i=0; i<nIndM; i++) // Random sampling
		{
			j = (int)(nIndM * unif_rand()); 
			//j = (int)(nIndM * randnum3(&RNGidum)); 
			//j = ChoixRep1[numBoot*nIndM + i] - 1; //////////////////////////////
			ChoixRepM[i+nIndM*numBoot] = IndM[j]-1; ChoixRepP[i+nIndM*numBoot] = IndP[j]-1;// 0 to nIndM-1   Base 0
		}
		for (i=0; i<nIndM1; i++) {IndMPv2[i+cNRNF*numBoot] = i + 1; IndMPv2[i+nIndM+cNRNF*numBoot] = i + nIndM + 1;} //IndM1v2, IndP1v2
		for (i=0; i<nIndM2; i++) {IndMPv2[i+nIndM1+cNRNF*numBoot] = i + nIndM1 + 1; IndMPv2[i+nIndM+nIndM1+cNRNF*numBoot] = i + nIndM + nIndM1 + 1;} //IndM2v2, IndP2v2
	}
	
	UseCoresNew = UseCores;
	double *TestNbReadTot, *TestProbMat, *TestDataBootStrap, *Testtempa, *Testtempb;
	do
	{
		AllocOk = 1; 
		TestM12Buff = malloc(UseCoresNew*nIndM*sizeof(int)); if (TestM12Buff == NULL) { Rprintf( "Not enough memory to allocate buffer TestM12Buff \n"); AllocOk = 0;}
		TestP12Buff = malloc(UseCoresNew*nIndM*sizeof(int)); if (TestP12Buff == NULL) { Rprintf( "Not enough memory to allocate buffer TestP12Buff \n"); AllocOk = 0;}
		TestMSortBuff = malloc(UseCoresNew*nIndM*sizeof(double)); if (TestMSortBuff == NULL) { Rprintf( "Not enough memory to allocate buffer TestMSortBuff \n"); AllocOk = 0;}
		TestPSortBuff = malloc(UseCoresNew*nIndM*sizeof(double)); if (TestPSortBuff == NULL) { Rprintf( "Not enough memory to allocate buffer TestPSortBuff \n"); AllocOk = 0;}
		Testtempa = malloc(UseCoresNew*sizeof(double)); if (Testtempa == NULL) { Rprintf( "Not enough memory to allocate buffer Test \n"); AllocOk = 0;}
		Testtempb = malloc(UseCoresNew*sizeof(double)); if (Testtempb == NULL) { Rprintf( "Not enough memory to allocate buffer Test \n"); AllocOk = 0;}
		TestProbMat = malloc(UseCoresNew*nIndM*sizeof(double)); if (TestProbMat == NULL) { Rprintf( "Not enough memory to allocate buffer TestProbMat \n"); AllocOk = 0;}
		TestDataBootStrap = malloc(UseCoresNew*cNRNF*sizeof(double)); if (TestDataBootStrap == NULL) { Rprintf( "Not enough memory to allocate buffer TestDataBootStrap \n"); AllocOk = 0;}
		TestNbReadTot = malloc(UseCoresNew*nIndM*sizeof(double)); if (TestNbReadTot == NULL) { Rprintf( "Not enough memory to allocate buffer TestNbReadTot \n"); AllocOk = 0;}
		free(TestM12Buff); free(TestP12Buff); free(TestPSortBuff); free(TestMSortBuff); free(TestDataBootStrap); free(TestProbMat); free(Testtempa); free(Testtempb); free(TestNbReadTot); 
		if (AllocOk == 0) {UseCoresNew = UseCoresNew -1; Rprintf( "WARNING - Not enough memory to allocate parallel buffers 2 - Lowering to %d threads \n", UseCoresNew);} //NBoot[0] = 0; 
	}
	while ((AllocOk == 0) && (UseCoresNew > 0));
	if (UseCoresNew ==0) {Nboot[0] = 0; Rprintf( "ERROR - Not enough memory to allocate parallel buffers 2 \n"); goto AbNormalEnd2;}
	#ifdef _OPENMP
        omp_set_num_threads(UseCoresNew); 
    #endif   

	#pragma omp parallel for \
	private(line, vecEC9B, vecBRNF, i, j, k, AleaSelect, numBoot)\
	schedule(static)
	for (line = 0; line<nNRNF[0]; line++)
	{
		double *MSortBuff, *PSortBuff, *NbReadTot, *ProbMat, *DataBootStrap, *tempa, *tempb;
		int *M12Buff, *P12Buff;
		int *IndM1v2, *IndM2v2, *IndP1v2, *IndP2v2;
		int TestRepeat;
		
		M12Buff = malloc(nIndM*sizeof(int)); if (M12Buff == NULL) { Rprintf( "Not enough memory to allocate buffer M12Buff \n");}
		P12Buff = malloc(nIndM*sizeof(int)); if (P12Buff == NULL) { Rprintf( "Not enough memory to allocate buffer P12Buff \n");}
		MSortBuff = malloc(nIndM*sizeof(double)); if (MSortBuff == NULL) { Rprintf( "Not enough memory to allocate buffer MSortBuff \n");}
		PSortBuff = malloc(nIndM*sizeof(double)); if (PSortBuff == NULL) { Rprintf( "Not enough memory to allocate buffer PSortBuff \n");}
		tempa = malloc(sizeof(double)); if (tempa == NULL) { Rprintf( "Not enough memory to allocate buffer tempa \n");}
		tempb = malloc(sizeof(double)); if (tempb == NULL) { Rprintf( "Not enough memory to allocate buffer tempb \n");}
		ProbMat = malloc(nIndM*sizeof(double)); if (ProbMat == NULL) { Rprintf( "Not enough memory to allocate buffer ProbMat \n");}
		DataBootStrap = malloc(cNRNF*sizeof(double)); if (DataBootStrap == NULL) { Rprintf( "Not enough memory to allocate buffer DataBootStrap \n");}
		NbReadTot = malloc(nIndM*sizeof(double)); if (NbReadTot == NULL) { Rprintf( "Not enough memory to allocate buffer NbReadTot \n");}
	
		vecBRNF = NbReadNormFil + line*cNRNF; // My_Thread = omp_get_thread_num(); 
		for (numBoot=0; numBoot<Nboot[0]; numBoot++)
		{
			//## standard bootstrap step
			vecEC9B = enrcrit9boot + line*Nboot[0] + numBoot; 
			do
			{
				if (WithRepeat[0] == 0)
				{ // without repetition
					for (i=0; i<nIndM; i++) // random selection of reads between fathers and mothers
					{
						do 
						{
							TestRepeat = 0;
							AleaSelect = unif_rand(); 
							j = 0; k = NbreVecProb; 
							while ((k-j) !=1)
							{
								if (AleaSelect < distdiffcumul[(int) floor((j+k)/2)])
									k = floor((j+k)/2);
								else
									j = floor((j+k)/2);
							}
							ProbMat[i] = vecprob[j]; 
							for (k = i-1; k>=0; k--) 
							{
								if (ProbMat[i]==ProbMat[k]) {TestRepeat = 1; break;}
							}
						} 
						while (TestRepeat!=0);
					}
				}
				else
				{  // with repetition
					for (i=0; i<nIndM; i++) // random selection of reads between fathers and mothers
					{
						AleaSelect = unif_rand(); 
						j = 0; k = NbreVecProb; 
						while ((k-j) !=1)
						{
							if (AleaSelect < distdiffcumul[(int) floor((j+k)/2)])
								k = floor((j+k)/2);
							else
								j = floor((j+k)/2);
						}
						ProbMat[i] = vecprob[j]; 
					}
				}			
				
				for (i=0; i<nIndM1; i++)  // Rebuild of samples 1
				{
					NbReadTot[i] = vecBRNF[ChoixRepM[i+nIndM*numBoot]] + vecBRNF[ChoixRepP[i+nIndM*numBoot]];
					DataBootStrap[i] = NbReadTot[i] * ProbMat[i]; 
					DataBootStrap[i+nIndM] = NbReadTot[i] - DataBootStrap[i];
				}
				for (i=0; i<nIndM2; i++) // Rebuild of samples 2
				{
					NbReadTot[i+nIndM1] = vecBRNF[ChoixRepM[i+nIndM1+nIndM*numBoot]] + vecBRNF[ChoixRepP[i+nIndM1+nIndM*numBoot]];
					DataBootStrap[i+nIndM1] = NbReadTot[i+nIndM1] * ProbMat[i+nIndM1];
					DataBootStrap[i+nIndM+nIndM1] = NbReadTot[i+nIndM1] - DataBootStrap[i+nIndM1];
				}
				IndM1v2 = IndMPv2+cNRNF*numBoot; IndM2v2 = IndMPv2+nIndM1+cNRNF*numBoot; IndP1v2 = IndMPv2+nIndM+cNRNF*numBoot; IndP2v2 = IndMPv2+nIndM+nIndM1+cNRNF*numBoot;
				Stat9NonParF(line, DataBootStrap, IndM1v2, IndM2v2, IndP1v2, IndP2v2, nIndM1, nIndM2, MSortBuff, PSortBuff, M12Buff, P12Buff, tempa, tempb, vecEC9B);
				
			} while (tempa[0]==0);
			
			if (isfinite(vecEC9B[0])) { moylocH0[line] = moylocH0[line] + vecEC9B[0]; NbBootOk[line]++; }
		}
		
		free(M12Buff); free(P12Buff); free(PSortBuff); free(MSortBuff); free(DataBootStrap); free(ProbMat); free(tempa); free(tempb); 
	}
	#ifdef _OPENMP
        omp_set_num_threads(UseCores);
    #endif
	nu0 = moylocH0[0] / NbBootOk[0]; 
	
	for (line = 0; line<nNRNF[0]; line++)
	{
		moylocH0[line] = moylocH0[line] / NbBootOk[line]; 
		if (moylocH0[line] < nu0) {nu0 = moylocH0[line];}
	}
	
	nu0 = floor(nu0 * 10) / 10; Rprintf("nu0= %f \n", nu0);	
	
	#pragma omp parallel for \
	private(line, vecEC9B, numBoot)\
	schedule(static)
	for (line = 0; line<nNRNF[0]; line++)
	{
		vecEC9B = enrcrit9boot + line*Nboot[0]; //vecEC9BZ = enrcrit9bootZ + line*Nboot[0];
		for (numBoot=0; numBoot<Nboot[0]; numBoot++)
		{
			vecEC9B[numBoot] = vecEC9B[numBoot] - moylocH0[line] - nu0; 
			if ((vecEC9B[numBoot]>=crit9glob[line]) && (isfinite(vecEC9B[0]))) {pvalboot[line] = pvalboot[line] + 1;}
		}
		pvalboot[line] = pvalboot[line] / NbBootOk[line];
	}
	
	AbNormalEnd2:
	PutRNGstate();	
	
	AbNormalEnd:
	free(NbBootOk); free(IndMPv2); // free(RNGidum); free(RNGinext); free(RNGinextp); free(RNGiff); free(RNGma);//RNG
	//free(IndMv2); free(IndM1v2); free(IndM2v2); free(IndPv2); free(IndP1v2); free(IndP2v2);
	free(ChoixRepM); free(ChoixRepP); free(IndM); free(IndM1); free(IndM2); free(IndP); free(IndP1); free(IndP2);
	
	return(R_NilValue);
}

//=====================================================================================================================================

SEXP IsoldeP2(SEXP R_NCore, SEXP R_NbReadNormFil, SEXP R_NRowReadNormFil, SEXP R_IndM, SEXP R_IndP, SEXP R_crit9glob, SEXP R_NbreadH1, SEXP R_NRowNbreadH1, SEXP R_enrdiffH1, 
	SEXP R_enrcrit9bootH1bis, SEXP R_nbreadtot, SEXP R_Seuil, SEXP R_Nboot2, SEXP R_moylocHX, SEXP R_pvalbootH1bis)

{
	double *NbReadNormFil = REAL(R_NbReadNormFil), *NbreadH1 = REAL(R_NbreadH1), *enrdiffH1 = REAL(R_enrdiffH1), *enrcrit9bootH1bis = REAL(R_enrcrit9bootH1bis);
	double *nbreadtotT = REAL(R_nbreadtot), *crit9glob = REAL(R_crit9glob);
	double *moylocHX = REAL(R_moylocHX), *pvalbootH1bis = REAL(R_pvalbootH1bis); //, *enrcrit9bootH1bisZ = REAL(R_enrcrit9bootH1bisZ);
	double *Seuil = REAL(R_Seuil);
	//double *ChoixRep3 = REAL(R_TChoixRep3);
	int *nNRNF = INTEGER(R_NRowReadNormFil), *nNRH1 = INTEGER(R_NRowNbreadH1), *Nboot2 = INTEGER(R_Nboot2), *NCore = INTEGER(R_NCore);

	int sNRNF = length(R_NbReadNormFil); int cNRNF = floor(sNRNF / nNRNF[0]);
	int i, j, line, col, numBoot; 
	double *vecBRNF, *vecNRH1, *vecEDH1, *vecEC9BH1B, *vecNRTT; //, *vecEC9BH1BZ
	double lambda0; //   tempc, tempd, tempe, intpart, AleaSelect
	
	int nIndM1 = length(VECTOR_ELT(R_IndM, 0)); int nIndM2 = length(VECTOR_ELT(R_IndM, 1)); int nIndM = nIndM1 + nIndM2; 
	int nIndP1 = nIndM1; int nIndP2 = nIndM2; int nIndP = nIndP1 + nIndP2; 
	int *IndM1v2, *IndM2v2, *IndP1v2, *IndP2v2;
	int *IndM, *IndM1, *IndM2, *IndP, *IndP1, *IndP2, *ChoixRepM, *ChoixRepP, *NbBootOk, *IndMPv2;	
	double	*NbReadTot;
	
	int AllocOk = 1;
	IndM = malloc(nIndM*sizeof(int)); if (IndM == NULL) { Rprintf( "Not enough memory to allocate buffer IndM \n"); AllocOk = 0;}
	IndM1 = malloc(nIndM1*sizeof(int)); if (IndM1 == NULL) { Rprintf( "Not enough memory to allocate buffer IndM1 \n"); AllocOk = 0;}
	IndM2 = malloc(nIndM2*sizeof(int)); if (IndM2 ==NULL) { Rprintf( "Not enough memory to allocate buffer IndM2 \n"); AllocOk = 0;}
	IndP = malloc(nIndP*sizeof(int)); if (IndP == NULL) { Rprintf( "Not enough memory to allocate buffer IndP \n"); AllocOk = 0;}
	IndP1 = malloc(nIndP1*sizeof(int)); if (IndP1 == NULL) { Rprintf( "Not enough memory to allocate buffer IndP1 \n"); AllocOk = 0;}
	IndP2 = malloc(nIndP2*sizeof(int)); if (IndP2 == NULL) { Rprintf( "Not enough memory to allocate buffer IndP2 \n"); AllocOk = 0;}
	IndMPv2 = malloc(cNRNF*sizeof(int)); if (IndMPv2 == NULL) { Rprintf( "Not enough memory to allocate buffer IndMPv2 \n"); AllocOk = 0;}

	ChoixRepM = malloc(nIndM*sizeof(int)); if (ChoixRepM == NULL) { Rprintf( "Not enough memory to allocate buffer ChoixRepM \n"); AllocOk = 0;}
	ChoixRepP = malloc(nIndM*sizeof(int)); if (ChoixRepP == NULL) { Rprintf( "Not enough memory to allocate buffer ChoixRepP \n"); AllocOk = 0;}
	NbReadTot = malloc(nIndM*sizeof(double)); if (NbReadTot == NULL) { Rprintf( "Not enough memory to allocate buffer NbReadTot \n"); AllocOk = 0;}
	NbBootOk = malloc(nNRNF[0]*sizeof(int)); if (NbBootOk == NULL) { Rprintf( "Not enough memory to allocate buffer NbBootOk \n"); AllocOk = 0;}
	if  (AllocOk == 0)
		{Nboot2[0] = 0; Rprintf( "ERROR - Not enough memory to allocate buffers P2 \n"); goto AbNormalEnd;}
		
	//Set the number of threads
	int HaveCores = 1;
    #ifdef _OPENMP
        HaveCores = omp_get_num_procs(); Rprintf("Open MP enabled \n");
    #else   
         Rprintf("Open MP disabled \n");
    #endif  
	Rprintf("Have Cores : %d \n", HaveCores);
	int UseCores, UseCoresNew;
	if (NCore[0]<0) {UseCores = 1;}
	else
		{
			if (NCore[0]>=100) {UseCores = HaveCores;}
			else {UseCores = (int) floor(HaveCores*NCore[0]/100);}
		}
	if (UseCores<1) {UseCores = 1;}
	#ifdef _OPENMP
        omp_set_num_threads(UseCores);
    #endif  
	Rprintf("UseCores : %d \n", UseCores);
	//int My_Thread;	
	
	for (i=0; i<nIndM1; i++) {IndM1[i] = INTEGER(VECTOR_ELT(R_IndM, 0))[i]; IndM[i] = IndM1[i]; IndP1[i] = INTEGER(VECTOR_ELT(R_IndP, 0))[i]; IndP[i] = IndP1[i];}
	for (i=0; i<nIndM2; i++) {IndM2[i] = INTEGER(VECTOR_ELT(R_IndM, 1))[i]; IndM[i+nIndM1] = IndM2[i]; IndP2[i] = INTEGER(VECTOR_ELT(R_IndP, 1))[i]; IndP[i+nIndP1] = IndP2[i];	}
	GetRNGstate();
	
	#pragma omp parallel for \
	private(line, vecNRH1, vecEDH1, i, j, numBoot)\
	schedule(static)
	for (line = 0; line<nNRH1[0]; line++)
	{
		vecNRH1 = NbreadH1 + line*cNRNF; vecEDH1 = enrdiffH1 + line*nIndM1*nIndM2; col = 0;
		for (j=0; j<nIndM2; j++) {NbReadTot[j]= vecNRH1[IndP2[j]-1]/(vecNRH1[IndP2[j]-1]+vecNRH1[IndM2[j]-1]);}
		for (i=0; i<nIndM1; i++)
		{
			for (j=0; j<nIndM2; j++) { vecEDH1[col] = (vecNRH1[IndM1[i]-1]/(vecNRH1[IndM1[i]-1]+vecNRH1[IndP1[i]-1]))-NbReadTot[j]; col = col + 1; }
		}
	}
	
	#pragma omp parallel for \
	private(line, vecNRTT, vecBRNF, i)\
	schedule(static)
	for (line = 0; line<nNRNF[0]; line++) // calcul nbreadtot
	{
		vecNRTT = nbreadtotT + line * nIndM; vecBRNF = NbReadNormFil + line*cNRNF;
		for (i=0; i<nIndM; i++) { vecNRTT[i] = vecBRNF[IndM[i]-1] + vecBRNF[IndP[i]-1];} 
		NbBootOk[line] = 0; 
	}
	
	for (i=0; i<nIndM1; i++) {IndMPv2[i] = i + 1; IndMPv2[i+nIndM] = i + nIndM + 1;} //IndM1v2, IndP1v2
	for (i=0; i<nIndM2; i++) {IndMPv2[i+nIndM1] = i + nIndM1 + 1; IndMPv2[i+nIndM+nIndM1] = i + nIndM + nIndM1 + 1;} //IndM2v2, IndP2v2
	
	UseCoresNew = UseCores;
	double *TestMSortBuff, *TestPSortBuff, *TestProbMat, *TestDataBootStrap, *Testtempa, *Testtempb;
	int *TestM12Buff, *TestP12Buff;
	do
	{
		AllocOk = 1; 
		TestM12Buff = malloc(UseCoresNew*nIndM*sizeof(int)); if (TestM12Buff == NULL) { Rprintf( "Not enough memory to allocate buffer TestM12Buff \n"); AllocOk = 0;}
		TestP12Buff = malloc(UseCoresNew*nIndM*sizeof(int)); if (TestP12Buff == NULL) { Rprintf( "Not enough memory to allocate buffer TestP12Buff \n"); AllocOk = 0;}
		TestMSortBuff = malloc(UseCoresNew*nIndM*sizeof(double)); if (TestMSortBuff == NULL) { Rprintf( "Not enough memory to allocate buffer TestMSortBuff \n"); AllocOk = 0;}
		TestPSortBuff = malloc(UseCoresNew*nIndM*sizeof(double)); if (TestPSortBuff == NULL) { Rprintf( "Not enough memory to allocate buffer TestPSortBuff \n"); AllocOk = 0;}
		Testtempa = malloc(UseCoresNew*sizeof(double)); if (Testtempa == NULL) { Rprintf( "Not enough memory to allocate buffer Test \n"); AllocOk = 0;}
		Testtempb = malloc(UseCoresNew*sizeof(double)); if (Testtempb == NULL) { Rprintf( "Not enough memory to allocate buffer Test \n"); AllocOk = 0;}
		TestProbMat = malloc(UseCoresNew*nIndM*sizeof(double)); if (TestProbMat == NULL) { Rprintf( "Not enough memory to allocate buffer TestProbMat \n"); AllocOk = 0;}
		TestDataBootStrap = malloc(UseCoresNew*cNRNF*sizeof(double)); if (TestDataBootStrap == NULL) { Rprintf( "Not enough memory to allocate buffer TestDataBootStrap \n"); AllocOk = 0;}
		free(TestM12Buff); free(TestP12Buff); free(TestPSortBuff); free(TestMSortBuff); free(TestDataBootStrap); free(TestProbMat); free(Testtempa); free(Testtempb); 
		if (AllocOk ==0) {UseCoresNew = UseCoresNew -1; Rprintf( "WARNING - Not enough memory to allocate parallel buffers P2 - Lowering to %d threads \n", UseCoresNew);} //NBoot[0] = 0; 
	}
	while ((AllocOk==0) && (UseCoresNew > 0));
	if (UseCoresNew ==0) {Nboot2[0] = 0; Rprintf( "ERROR - Not enough memory to allocate parallel buffers P2 \n"); goto AbNormalEnd2;}
	#ifdef _OPENMP
        omp_set_num_threads(UseCoresNew); 
    #endif
	
	#pragma omp parallel for \
	private(line, numBoot, vecNRTT, vecEC9BH1B, i)\
	schedule(static)
	for (line = 0; line<nNRNF[0]; line++)//line/boot
	{
		double *MSortBuff, *PSortBuff, *ProbMat, *DataBootStrap, *tempa, *tempb;
		int *M12Buff, *P12Buff;
		
		M12Buff = malloc(nIndM*sizeof(int)); if (M12Buff == NULL) { Rprintf( "Not enough memory to allocate buffer M12Buff \n");}
		P12Buff = malloc(nIndM*sizeof(int)); if (P12Buff == NULL) { Rprintf( "Not enough memory to allocate buffer P12Buff \n");}
		MSortBuff = malloc(nIndM*sizeof(double)); if (MSortBuff == NULL) { Rprintf( "Not enough memory to allocate buffer MSortBuff \n");}
		PSortBuff = malloc(nIndM*sizeof(double)); if (PSortBuff == NULL) { Rprintf( "Not enough memory to allocate buffer PSortBuff \n");}
		tempa = malloc(sizeof(double)); if (tempa == NULL) { Rprintf( "Not enough memory to allocate buffer tempa \n");}
		tempb = malloc(sizeof(double)); if (tempb == NULL) { Rprintf( "Not enough memory to allocate buffer tempb \n");}
		ProbMat = malloc(nIndM*sizeof(double)); if (ProbMat == NULL) { Rprintf( "Not enough memory to allocate buffer ProbMat \n");}
		DataBootStrap = malloc(cNRNF*sizeof(double)); if (DataBootStrap == NULL) { Rprintf( "Not enough memory to allocate buffer DataBootStrap \n");}
	
		vecNRTT = nbreadtotT + line * nIndM; //My_Thread = omp_get_thread_num(); 
		for (numBoot=0; numBoot<Nboot2[0]; numBoot++) //		# standard bootstrap step
		{
			vecEC9BH1B = enrcrit9bootH1bis + line*Nboot2[0] + numBoot;
			for (i=0; i<nIndM; i++) 
			{ 
				ProbMat[i] = Seuil[0] + unif_rand()*(1-Seuil[0]);  // random read selection between fathers and mothers
				//ProbMat[i] = ChoixRep3[numBoot*nIndM*nNRNF[0] + nIndM*line + i] ;////////////////////////////////////////////
				DataBootStrap[i] = vecNRTT[i] * ProbMat[i]; 
				DataBootStrap[i+nIndM] = vecNRTT[i] - DataBootStrap[i]; 
			}	
			IndM1v2 = IndMPv2; IndM2v2 = IndMPv2+nIndM1; IndP1v2 = IndMPv2+nIndM; IndP2v2 = IndMPv2+nIndM+nIndM1;
			Stat9NonParF(0, DataBootStrap, IndM1v2, IndM2v2, IndP1v2, IndP2v2, nIndM1, nIndM2, MSortBuff, PSortBuff, M12Buff, P12Buff, tempa, tempb, vecEC9BH1B);
			if (isfinite(vecEC9BH1B[0])) {moylocHX[line] = moylocHX[line] + vecEC9BH1B[0]; NbBootOk[line]++;}
		}
		free(M12Buff); free(P12Buff); free(PSortBuff); free(MSortBuff); free(DataBootStrap); free(ProbMat); free(tempa); free(tempb);
	}
	#ifdef _OPENMP
        omp_set_num_threads(UseCoresNew); 
    #endif  
		
	lambda0 = moylocHX[0] / NbBootOk[0];
	for (line = 0; line<nNRNF[0]; line++)
	{
		moylocHX[line] = moylocHX[line] / NbBootOk[line]; 
		if (moylocHX[line] < lambda0) {lambda0 = moylocHX[line];}
	}
	lambda0 = floor(lambda0 * 10) / 10;
	Rprintf("lambda0 : %f \n", lambda0);

	#pragma omp parallel for \
	private(line, vecEC9BH1B, numBoot)\
	schedule(static)
	for (line = 0; line<nNRNF[0]; line++)
	{
		vecEC9BH1B = enrcrit9bootH1bis + line*Nboot2[0]; NbBootOk[line] = 0; //vecEC9BH1BZ = enrcrit9bootH1bisZ + line*Nboot2[0]; 
		for (numBoot=0; numBoot<Nboot2[0]; numBoot++)
		{
			//vecEC9BH1BZ[numBoot] = vecEC9BH1B[numBoot];
			vecEC9BH1B[numBoot] = vecEC9BH1B[numBoot] - moylocHX[line] + lambda0; 
			if (vecEC9BH1B[numBoot]<=crit9glob[line]) {pvalbootH1bis[line] = pvalbootH1bis[line] + 1; NbBootOk[line] = NbBootOk[line] + 1; }
		}
		pvalbootH1bis[line] = pvalbootH1bis[line] / Nboot2[0];
	}
	AbNormalEnd2:
	PutRNGstate();	
	//free(RDGidum); free(RDGidum2); free(RDGiy); free(RDGiv); 
	AbNormalEnd:
	free(NbBootOk); 
	free(NbReadTot); free(ChoixRepM); free(ChoixRepP); free(IndM); free(IndM1); free(IndM2); free(IndP); free(IndP1); free(IndP2); free(IndMPv2);
	
	return(R_NilValue);
	}

//=====================================================================================================================================
	
void Stat9NonParF(int vline, double *vec, int *iM1, int *iM2, int *iP1, int *iP2, int nM1, int nM2, double *MSorted, double *PSorted, int *M12, int *P12, double *denom, double *diffp, double *crit)
{
	int i, j, nM; //k
	double sM1 = 0, sM2 = 0, sP1 = 0, sP2 = 0;
	double a, meM, meP, sM, sP, meDiffM, meDiffP, CVM, CVP;
	nM = nM1 + nM2; //nMP = 2 * nM; vData = vline*nMP;
	
	for (i=0; i<nM1; i++) {sM1 = sM1 + vec[iM1[i]-1]; sP1 = sP1 + vec[iP1[i]-1]; M12[i] = iM1[i]-1; P12[i] = iP1[i]-1;}
	for (i=0; i<nM2; i++) {sM2 = sM2 + vec[iM2[i]-1]; sP2 = sP2 + vec[iP2[i]-1]; M12[i+nM1] = iM2[i]-1; P12[i+nM1] = iP2[i]-1;}
	diffp[0] = sM1/(sM1+sP1) - sP2/(sM2+sP2); sP = sP1 + sP2; sM = sM1 + sM2; 
	MSorted[0] = vec[M12[0]]; // tri M
	for (j=2; j<=nM; j++) //Pick out each element in turn.
	{ 
		a = vec[M12[j-1]]; i = j-1;
		while (i > 0 && MSorted[i-1] > a) {MSorted[i] = MSorted[i-1]; i--;} //     Look for the place to insert it.
		MSorted[i] = a; //Insert it.
	}
	
	PSorted[0] = vec[P12[0]];// tri P
	for (j=2; j<=nM; j++) //Pick out each element in turn.
	{ 
		a = vec[P12[j-1]]; i = j-1;
		while (i > 0 && PSorted[i-1] > a) {PSorted[i] = PSorted[i-1]; i--;} //     Look for the place to insert it.
		PSorted[i] = a; //Insert it.
	}
	if (nM % 2 == 1) //calcul median M and P
		{meM = MSorted[(int) floor(nM/2)]; meP = PSorted[(int)floor(nM/2)];}
	else 
		{meM = (MSorted[(int) floor(nM/2)-1] + MSorted[(int) floor(nM/2)])/2; meP = (PSorted[(int) floor(nM/2)-1] + PSorted[(int) floor(nM/2)])/2;}
	for (i=0; i<nM; i++) {MSorted[i] = fabs(MSorted[i] - meM); PSorted[i] = fabs(PSorted[i] - meP);}
	if (sM == 0) CVM = 0;
	else 
	{	
		for (j=2; j<=nM; j++) //tri  PSorted
		{ 
			a = MSorted[j-1]; i = j-1;
			while (i > 0 && MSorted[i-1] > a) {MSorted[i] = MSorted[i-1]; i--;} //     Look for the place to insert it.
			MSorted[i] = a; //Insert it.
		}
		if (nM % 2 == 1) {meDiffM = MSorted[(int) floor(nM/2)];} else {meDiffM = (MSorted[(int) floor(nM/2)-1] + MSorted[(int) floor(nM/2)])/2;}
		CVM = meDiffM/sM;
	}
	if (sP == 0) CVP = 0;
	else 
	{	
		for (j=2; j<=nM; j++) //tri  PSorted
		{ 
			a = PSorted[j-1]; i = j-1;
			while (i > 0 && PSorted[i-1] > a) {PSorted[i] = PSorted[i-1]; i--;} //     Look for the place to insert it.
			PSorted[i] = a; //Insert it.
		}
		if (nM % 2 == 1) {meDiffP = PSorted[(int) floor(nM/2)];} else {meDiffP = (PSorted[(int) floor(nM/2)-1] + PSorted[(int) floor(nM/2)])/2;}
		CVP = meDiffP/sP;
	}
	denom[0] = CVM + CVP; crit[0] = fabs(diffp[0])/sqrt(denom[0]);
} 

R_CallMethodDef callMethods[] = {
         {"Stat9NonParF", (DL_FUNC) &Stat9NonParF, 15},
         {NULL, NULL, 0}
};

void R_init_ISoLDE(DllInfo *info)
{
         R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

//=====================================================================================================================================

