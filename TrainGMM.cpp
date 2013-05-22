#include "mex.h"
#include <iostream>
#include <opencv/cv.h>
#include <opencv/cxcore.h>
#include <opencv/highgui.h>
#include <opencv/ml.h>
#include <string>
using namespace cv;
int train(float* src,char* filename,int N);
void mexFunction(
int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[])
{
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
        mexErrMsgTxt("input: ([3 N] samples),(o),(f) ");
    /*
    int* o = (int *)mxGetData(prhs[1]);
	int* f = (int *)mxGetData(prhs[2]);
    int namelen = (*o)*(*f);
    char filename [50];
    sprintf(filename,"%d_%d",*o,*f);
    */
    float* src = (float*)mxGetData(prhs[0]);
	int C = mxGetM(prhs[0]);
	int N = mxGetN(prhs[0]);
    if(!mxIsSingle(prhs[0]))
        mexErrMsgTxt("input should be float");
    if(C != 3)
        mexErrMsgTxt("dimension error");
    train(src,filename,N);
    /*
	float* b;
	plhs[0] = mxCreateNumericMatrix(4,2,mxSINGLE_CLASS,mxREAL);
	b = (float *)mxGetData(plhs[0]);
	b[0] = float(src[1]);
	b[1] = float(src[2]);
    b[2] = 3;
    b[3] = 4;*/
}



int train(float* src,char *filename,int N)
{
    EM* myGMM = new EM(5);
	FileStorage tmp(filename,cv::FileStorage::WRITE);
    Mat train_32F(N,3,CV_32FC1,src);
    train_32F.convertTo(train_32F,CV_32F,1/255.);
    myGMM->train(train_32F);
    tmp<<"model";
	myGMM->write(tmp);
	tmp.release();
	
    return 0;
}
