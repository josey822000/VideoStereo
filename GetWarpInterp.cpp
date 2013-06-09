#include "mex.h"
#ifdef _OPENMP
#include <omp.h>
#endif
/* using input position to warp the point to the corresponding frame */
template< class T>int Warp(T* F,T* segMap,double* D,T* dstF,T* dstSeg,double* dstDis,double* X,double* Y,int szX,int szY);
void mexFunction(
int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{
	/*  input   */
    /*  input   */
    

    
    if(nrhs<3)
        mexErrMsgTxt("input: (src),(x),(y)");
    if(nlhs<1)
        mexErrMsgTxt("output: (dst)");
    /*
    int* o = (int *)mxGetData(prhs[1]);
	int* f = (int *)mxGetData(prhs[2]);
    int namelen = (*o)*(*f);
    char filename [50];
    sprintf(filename,"%d_%d",*o,*f);
    */
    int szY = mxGetM(prhs[0]);
	int szX = mxGetN(prhs[0]);
    double* X = (double*)mxGetData(prhs[3]);
    double* Y = (double*)mxGetData(prhs[4]);
    double* D = (double*)mxGetData(prhs[2]);
    if(mxIsDouble(prhs[0]))
    {
        double* F = (double*)mxGetData(prhs[0]);
        double* segMap = (double*)mxGetData(prhs[1]);
        plhs[0] = mxCreateNumericMatrix(szY, szX, mxDOUBLE_CLASS, mxREAL);
        double* dstF = (double*)mxGetData(plhs[0]);
        plhs[1] = mxCreateNumericMatrix(szY, szX, mxDOUBLE_CLASS, mxREAL);
        double* dstSeg = (double*)mxGetData(plhs[1]);
        plhs[2] = mxCreateNumericMatrix(szY, szX, mxDOUBLE_CLASS, mxREAL);
        double* dstDis = (double*)mxGetData(plhs[2]);
        Warp(F,segMap,D,dstF,dstSeg,dstDis,X,Y,szX,szY);
    }
    else{
        unsigned int * F = (unsigned int *)mxGetData(prhs[0]);
        unsigned int * segMap = (unsigned int *)mxGetData(prhs[1]);
        
        plhs[0] = mxCreateNumericMatrix(szY, szX, mxUINT32_CLASS, mxREAL);
        unsigned int * dstF = (unsigned int *)mxGetData(plhs[0]);
        plhs[1] = mxCreateNumericMatrix(szY, szX, mxUINT32_CLASS, mxREAL);
        unsigned int * dstSeg = (unsigned int *)mxGetData(plhs[1]);
        plhs[2] = mxCreateNumericMatrix(szY, szX, mxDOUBLE_CLASS, mxREAL);
        double* dstDis = (double*)mxGetData(plhs[2]);
        Warp(F,segMap,D,dstF,dstSeg,dstDis,X,Y,szX,szY);
    }

    
    
    
	
    
    
    //dst[0] = 1;
    //dst[1] = 2;
}



template< class T>int Warp(T* F,T* segMap,double* D,T* dstF,T* dstSeg,double* dstDis,double* X,double* Y,int szX,int szY)
{
	int id = 0,idx = 0,idy = 0,cnt=0;    
    #pragma omp parallel for
    for(int j = 0;j<szX;j++)
        for(int i=0;i<szY;i++)
        {
            id = j*szY+i;
            idx = X[id];
            idy = Y[id];
            if(idx>=0 && idx<szX && idy>=0 && idy<szY && D[id]>dstDis[idx*szY+idy])
            {
                dstF[idx*szY+idy] = F[id];
                dstSeg[idx*szY+idy] = segMap[id];
                dstDis[idx*szY+idy] = D[id];
            }
        }
    double tmpD;
    #pragma omp parallel for
    for(int j = 1;j<szX-1;j++)
        for(int i=1;i<szY-1;i++)
        {
            id = j*szY+i;
            if(dstF[id]==0)
            {
                tmpD = 0;
                cnt = 0;
                if(dstF[id+1]!=0)
                {
                    dstSeg[id] = dstSeg[id+1];
                    tmpD += dstDis[id+1];
                    cnt++;
                }
                if(dstF[id-1]!=0)
                {
                    dstSeg[id] = dstSeg[id-1];
                    tmpD += dstDis[id-1];
                    cnt++;
                }
                if(dstF[id+szY]!=0)
                {
                    dstSeg[id] = dstSeg[id+szY];
                    tmpD += dstDis[id+szY];
                    cnt++;
                }
                if(dstF[id-szY]!=0)
                {
                    dstSeg[id] = dstSeg[id-szY];
                    tmpD += dstDis[id-szY];
                    cnt++;
                }
                if(cnt !=0)
                {
                    dstDis[id] = tmpD/cnt;
                }
            }
        }
    
    for(int i=1;i<szY-1;i++)
    {
        id = i;
        if(dstF[id]==0)
        {
            
            tmpD = 0;
            cnt = 0;
            if(dstF[id+1]!=0)
            {
                dstSeg[id] = dstSeg[id+1];
                tmpD += dstDis[id+1];
                cnt++;
            }
            if(dstF[id-1]!=0)
            {
                dstSeg[id] = dstSeg[id-1];
                tmpD += dstDis[id-1];
                cnt++;
            }
            if(dstF[id+szY]!=0)
            {
                dstSeg[id] = dstSeg[id+szY];
                tmpD += dstDis[id+szY];
                cnt++;
            }
            if(cnt !=0)
            {
                dstDis[id] = tmpD/cnt;
            }
        }
        
        //
        id = (szX-1)*szY+i;
        if(dstF[id]==0)
        {
            tmpD = 0;
            cnt = 0;
            if(dstF[id+1]!=0)
            {
                dstSeg[id] = dstSeg[id+1];
                tmpD += dstDis[id+1];
                cnt++;
            }
             
            if(dstF[id-1]!=0)
            {
                dstSeg[id] = dstSeg[id-1];
                tmpD += dstDis[id-1];
                cnt++;
            }
            if(dstF[id-szY]!=0)
            {
                dstSeg[id] = dstSeg[id-szY];
                tmpD += dstDis[id-szY];
                cnt++;
            }
            if(cnt !=0)
            {
                dstDis[id] = tmpD/cnt;
            }
        }
        
    }
    //
    for(int j=1;j<szX-1;j++)
    {
        id = j*szY;
        if(dstF[id]==0)
        {
            tmpD = 0;
            cnt = 0;
            if(dstF[id+1]!=0)
            {
                dstSeg[id] = dstSeg[id+1];
                tmpD += dstDis[id+1];
                cnt++;
            }
            if(dstF[id+szY]!=0)
            {
                dstSeg[id] = dstSeg[id+szY];
                tmpD += dstDis[id+szY];
                cnt++;
            }
            if(dstF[id-szY]!=0)
            {
                dstSeg[id] = dstSeg[id-szY];
                tmpD += dstDis[id-szY];
                cnt++;
            }
            if(cnt !=0)
            {
                dstDis[id] = tmpD/cnt;
            }
        }
        
        //
        id = j*szY+(szY-1);
        if(dstF[id]==0)
        {
            tmpD = 0;
            cnt = 0;
            if(dstF[id-1]!=0)
            {
                dstSeg[id] = dstSeg[id-1];
                tmpD += dstDis[id-1];
                cnt++;
            }
            if(dstF[id+szY]!=0)
            {
                dstSeg[id] = dstSeg[id+szY];
                tmpD += dstDis[id+szY];
                cnt++;
            }
            if(dstF[id-szY]!=0)
            {
                dstSeg[id] = dstSeg[id-szY];
                tmpD += dstDis[id-szY];
                cnt++;
            }
            if(cnt !=0)
            {
                dstDis[id] = tmpD/cnt;
            }
        }
        
    }
    id = 0;
    if(dstF[id]==0)
    {
        tmpD = 0;
        cnt = 0;
        if(dstF[id+1]!=0)
        {
            dstSeg[id] = dstSeg[id+1];
            tmpD += dstDis[id+1];
            cnt++;
        }
        
        if(dstF[id+szY]!=0)
        {
            dstSeg[id] = dstSeg[id+szY];
            tmpD += dstDis[id+szY];
            cnt++;
        }
        
        if(cnt !=0)
        {
            dstDis[id] = tmpD/cnt;
        }
    }
    id = szY-1;
    if(dstF[id]==0)
    {
        
        tmpD = 0;
        cnt = 0;
        if(dstF[id-1]!=0)
        {
            dstSeg[id] = dstSeg[id-1];
            tmpD += dstDis[id-1];
            cnt++;
        }
        if(dstF[id+szY]!=0)
        {
            dstSeg[id] = dstSeg[id+szY];
            tmpD += dstDis[id+szY];
            cnt++;
        }
        
        if(cnt !=0)
        {
            dstDis[id] = tmpD/cnt;
        }
    }
    id = (szX-1)*szY;
    if(dstF[id]==0)
    {
        
        tmpD = 0;
        cnt = 0;
        if(dstF[id+1]!=0)
        {
            dstSeg[id] = dstSeg[id+1];
            tmpD += dstDis[id+1];
            cnt++;
        }
        
        if(dstF[id-szY]!=0)
        {
            dstSeg[id] = dstSeg[id-szY];
            tmpD += dstDis[id-szY];
            cnt++;
        }
        if(cnt !=0)
        {
            dstDis[id] = tmpD/cnt;
        }
    }
    id = szX*szY-1;
    if(dstF[id]==0)
    {
        tmpD = 0;
        cnt = 0;
        
        if(dstF[id-1]!=0)
        {
            dstSeg[id] = dstSeg[id-1];
            tmpD += dstDis[id-1];
            cnt++;
        }
        
        if(dstF[id-szY]!=0)
        {
            dstSeg[id] = dstSeg[id-szY];
            tmpD += dstDis[id-szY];
            cnt++;
        }
        if(cnt !=0)
        {
            dstDis[id] = tmpD/cnt;
        }
    }
    return 0;
}
