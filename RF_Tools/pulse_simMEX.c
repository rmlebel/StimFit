#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *left[], int nrhs, const mxArray *right[]) {
    
    /*  Declare variables */
    mwSize Nz, Nt, ETL, iNz, iNt, indNt, ind, iETL, dimsA[2];
    const mwSize *dims;
    double *pMin, *pCOS, *pSIN, *pCTRF, *pSTRF, *pcpRF, *pspRF, *pMout, *pAout;
    double cpRF, spRF;
    double Mx, My, Mz;
    mxArray *Mout, *Aout;
    
    /*  Get pointers to input */
    pMin  = mxGetPr(right[0]);
    pCOS  = mxGetPr(right[1]);
    pSIN  = mxGetPr(right[2]);
    pCTRF = mxGetPr(right[3]);
    pSTRF = mxGetPr(right[4]);
    pcpRF = mxGetPr(right[5]);cpRF = pcpRF[0];
    pspRF = mxGetPr(right[6]);spRF = pspRF[0];
    
    
    /*  Determine number of elements */
    Nt = mxGetN(right[3]);
    dims = mxGetDimensions(right[0]);
    ETL = dims[1];
    Nz  = dims[2];
    /*mexPrintf("Nt: %i; ETL = %i; Nz = %i\n",Nt,ETL,Nz);*/
    
    /*  Copy input magnetization */
    Mout  = mxDuplicateArray(right[0]);
    pMout = mxGetPr(Mout);
    
    /*  Create angle array and assign pointer */
    dimsA[0] = ETL;dimsA[1] = Nz;
    Aout = mxCreateNumericArray(2,dimsA,mxDOUBLE_CLASS,mxREAL);
    pAout = mxGetPr(Aout);
    
    
    /*  Perfom magnetization simulation */
    /*  Loop through all time points */
    for (iNt=0; iNt<Nt; iNt++) {
        
        /*  Apply gradient rotation */
        /*  Loop through spatial positions */
        for (iNz=0; iNz<Nz; iNz++) {
            /*  Loop through all echo times */
            for (iETL=0;iETL<ETL; iETL++) {
                ind = 3*ETL*iNz + 3*iETL;
                Mx = pMout[ind];
                My = pMout[ind+1];
                pMout[ind]   =  pCOS[iNz]*Mx + pSIN[iNz]*My;
                pMout[ind+1] = -pSIN[iNz]*Mx + pCOS[iNz]*My;
            }
        }
        
        /*  Apply RF tip */
        /*  Loop through spatial positions */
        for (iNz=0; iNz<Nz; iNz++) {
            /*  Loop through all echo times */
            for (iETL=0;iETL<ETL; iETL++) {
                indNt = ETL*iNt + iETL;
                ind = 3*ETL*iNz + 3*iETL;
                Mx = pMout[ind];
                My = pMout[ind+1];
                Mz = pMout[ind+2];
                
                pMout[ind]   = Mx*(cpRF*cpRF + pCTRF[indNt]*spRF*spRF) +
                               My*(cpRF*spRF - pCTRF[indNt]*cpRF*spRF) -
                               Mz*(spRF*pSTRF[indNt]);
                pMout[ind+1] = My*(pCTRF[indNt]*cpRF*cpRF + spRF*spRF) +
                               Mx*(cpRF*spRF - cpRF*spRF*pCTRF[indNt]) +
                               Mz*(cpRF*pSTRF[indNt]);
                pMout[ind+2] = Mz*pCTRF[indNt] -
                               My*cpRF*pSTRF[indNt] +
                               Mx*spRF*pSTRF[indNt];
            }
        }
    } /* End RF time points */
    
    
    /*  Compute angle */
    /*  Loop through spatial positions */
    for (iNz=0; iNz<Nz; iNz++) {
        /*  Loop through all echo times */
        for (iETL=0;iETL<ETL; iETL++) {
            ind = 3*ETL*iNz + 3*iETL + 2; /* Index for z magnetization */
            indNt = ETL*iNz + iETL; /* Index into angle array */
            Mz = pMout[ind];
            Mz = (Mz <= -1.0) ? -0.9999999999999 : Mz;
            Mz = (Mz >=  1.0) ?  0.9999999999999 : Mz;
            pAout[indNt] = acos(Mz);
        }
    }
            
    
    /*  Assign outputs */
    left[0] = Mout;
    left[1] = Aout;
}