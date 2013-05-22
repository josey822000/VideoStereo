#include "mex.h"

template< class T>int Warp(T* D1,T* D2,T* dstD1,T* dstD2,double* X1,double* Y1,double* X2,double* Y2,int szX,int szY);
void mexFunction(
int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{
	/*  input   */
    /*  input   */
    

    
    if(nrhs<3)
        mexErrMsgTxt("input: (D1),(D2),(x),(y)");
    if(nlhs<2)
        mexErrMsgTxt("output: (D1),(D2)");
    /*
    int* o = (int *)mxGetData(prhs[1]);
	int* f = (int *)mxGetData(prhs[2]);
    int namelen = (*o)*(*f);
    char filename [50];
    sprintf(filename,"%d_%d",*o,*f);
    */
    int szY = mxGetM(prhs[0]);
	int szX = mxGetN(prhs[0]);
    double* X1 = (double*)mxGetData(prhs[2]);
    double* Y1 = (double*)mxGetData(prhs[3]);
    double* X2 = (double*)mxGetData(prhs[4]);
    double* Y2 = (double*)mxGetData(prhs[5]);
    
    if(mxIsDouble(prhs[0]))
    {
        double* D1 = (double*)mxGetData(prhs[0]);
        double* D2 = (double*)mxGetData(prhs[1]);
        plhs[0] = mxCreateNumericMatrix(szY, szX, mxDOUBLE_CLASS, mxREAL);
        double* dstD1 = (double*)mxGetData(plhs[0]);
        plhs[1] = mxCreateNumericMatrix(szY, szX, mxDOUBLE_CLASS, mxREAL);
        double* dstD2 = (double*)mxGetData(plhs[1]);
        Warp(D1,D2,dstD1,dstD2,X1,Y1,X2,Y2,szX,szY);
    }
    else{
        unsigned char* D1 = (unsigned char*)mxGetData(prhs[0]);
        unsigned char* D2 = (unsigned char*)mxGetData(prhs[1]);
        
        plhs[0] = mxCreateNumericMatrix(szY, szX, mxUINT8_CLASS, mxREAL);
        unsigned char* dstD1 = (unsigned char*)mxGetData(plhs[0]);
        plhs[1] = mxCreateNumericMatrix(szY, szX, mxUINT8_CLASS, mxREAL);
        unsigned char* dstD2 = (unsigned char*)mxGetData(plhs[1]);
        Warp(D1,D2,dstD1,dstD2,X1,Y1,X2,Y2,szX,szY);
    }

    
    
    
	
    
    
    //dst[0] = 1;
    //dst[1] = 2;
}



template< class T>int Warp(T* D1,T* D2,T* dstD1,T* dstD2,double* X1,double* Y1,double* X2,double* Y2,int szX,int szY)
{
	int id = 0,idx = 0,idy = 0;    
    for(int j = 0;j<szX;j++)
        for(int i=0;i<szY;i++)
        {
            id = j*szY+i;
            idx = X1[id];
            idy = Y1[id];
            if(idx>=0 && idx<szX && idy>=0 && idy<szY )
                dstD1[idx*szY+idy] = D1[id];
            idx = X2[id];
            idy = Y2[id];
            if(idx>=0 && idx<szX && idy>=0 && idy<szY )
                dstD2[idx*szY+idy] = D2[id];
        }
	
    for(int i=0;i<szY;i++)
    {
        T tmp1=0,tmp2=0;
        for(int j = 0;j<szX;j++)
        {
            //
            if(dstD1[j*szY+i]==0)
                dstD1[j*szY+i] = tmp1;
            else
                tmp1 = dstD1[j*szY+i];
            //
            if(dstD2[j*szY+i]==0)
                dstD2[j*szY+i] = tmp2;
            else
                tmp2 = dstD2[j*szY+i];
        }
    }
    for(int i=0;i<szY;i++)
    {
        T tmp1=0,tmp2=0;
        for(int j = szX-1;j>=0;j--)
        {
            //
            if(dstD1[j*szY+i]==0)
                dstD1[j*szY+i] = tmp1;
            else
                tmp1 = dstD1[j*szY+i];
            //
            if(dstD2[j*szY+i]==0)
                dstD2[j*szY+i] = tmp2;
            else
                tmp2 = dstD2[j*szY+i];
        }
    }
    for(int j = 0;j<szX;j++)
    {
        T tmp1=0,tmp2=0;
        for(int i=0;i<szY;i++)
        {
            //
            if(dstD1[j*szY+i]==0)
                dstD1[j*szY+i] = tmp1;
            else
                tmp1 = dstD1[j*szY+i];
            //
            if(dstD2[j*szY+i]==0)
                dstD2[j*szY+i] = tmp2;
            else
                tmp2 = dstD2[j*szY+i];
        }
    }
    for(int j = 0;j<szX;j++)
    {
        T tmp1=0,tmp2=0;
        for(int i=szY-1;i>=0;i--)
        {
            //
            if(dstD1[j*szY+i]==0)
                dstD1[j*szY+i] = tmp1;
            else
                tmp1 = dstD1[j*szY+i];
            //
            if(dstD2[j*szY+i]==0)
                dstD2[j*szY+i] = tmp2;
            else
                tmp2 = dstD2[j*szY+i];
        }
    }
	
    return 0;
}
