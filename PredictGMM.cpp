#include "mex.h"

#include <opencv/cv.h>
#include <opencv/cxcore.h>
#include <opencv/highgui.h>
#include <opencv/ml.h>
using namespace cv;
int predict(float* src,double* dst,char *filename,int N);
void mexFunction(
int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{
	/*  input   */
    /*  input   */
    char *input_buf, *output_buf;
    int   buflen,status;


    /* Get the length of the input string. */
    buflen = (mxGetM(prhs[1]) * mxGetN(prhs[1])) + 1;

    /* Allocate memory for input and output strings. */
    char* filename = (char*)mxCalloc(buflen, sizeof(char));

    /* Copy the string data from prhs[0] into a C string 
    * input_buf. */
    status = mxGetString(prhs[1], filename, buflen);
    if(nrhs<2)
        mexErrMsgTxt("input: ([N 3] float samples),(modelName)   ");
    if(nlhs<1)
        mexErrMsgTxt("output: (double likelihood)");
    /*
    int* o = (int *)mxGetData(prhs[1]);
	int* f = (int *)mxGetData(prhs[2]);
    int namelen = (*o)*(*f);
    char filename [50];
    sprintf(filename,"%d_%d",*o,*f);
    */
    if(!mxIsSingle(prhs[0]))
        mexErrMsgTxt("input should be float");
    float* src = (float*)mxGetData(prhs[0]);
	int C = mxGetM(prhs[0]);
	int N = mxGetN(prhs[0]);
    if(C != 3)
        mexErrMsgTxt("dimension error");
    plhs[0] = mxCreateDoubleMatrix(1, N, mxREAL);
    double* dst = mxGetPr(plhs[0]);
    //plhs[1] = mxCreateDoubleMatrix(1, N, mxREAL);
    //double* dst2 = mxGetPr(plhs[1]);
    predict(src,dst,filename,N);
    //dst[0] = 1;
    //dst[1] = 2;
}



int predict(float* src,double* dst,char *filename,int N)
{
    
    EM* myGMM = new EM(5,EM::COV_MAT_SPHERICAL);
	FileStorage tmp(filename,cv::FileStorage::READ);
    Mat test_32F(N,3,CV_32FC1,src);
    test_32F = test_32F * (1/255.);
    //test_32F.convertTo(test_32F,CV_32F,1/255.);
    FileNode myFilenode = tmp["model"];
	myGMM->read(myFilenode);
	Mat sample(1,3,CV_32FC1);
    int height = test_32F.rows;
    for(int i=0;i<height;i++)
    {
        sample= test_32F.row(i);
        Vec2d id = myGMM->predict(sample);
        dst[i] = id.val[0];
            //dst[i] = myGMM->predict(sample);
    }
	
    return 0;
}
