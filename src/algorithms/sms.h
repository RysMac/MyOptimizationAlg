#ifndef _SMS_H
#define _SMS_H

#if defined(WIN32) || defined(__WIN32__) || defined(_WIN32) || defined(WIN64) || defined(__WIN64__) || defined(_WIN64)
#define SMSWINDOWS
#else
#define SMSLINUX
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#ifdef SMSWINDOWS
#include "direct.h"
#endif 

#ifdef SMSLINUX
#include <unistd.h>
#include <sys/dir.h>
#include <sys/param.h>
#include <dlfcn.h>
#endif 

#ifndef smsmin
#define smsmin(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef smsmax
#define smsmax(a,b) (((a)>(b))?(a):(b)) 
#endif
#ifndef min
#define min(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef max
#define max(a,b) (((a)>(b))?(a):(b)) 
#endif

#ifdef SMSWINDOWS
#define TIMEB _timeb
#define ISNAN _isnan
#define FINITE _finite
#define DLLEXPORT __declspec(dllexport)
#define FTIME _ftime
#define GETCWD _getcwd
#endif

#ifdef SMSLINUX
#define GETCWD getcwd
#define _copysign copysign
#define FTIME ftime
#define TIMEB timeb
#define _MAX_PATH 4096
#define ISNAN isnan
#define FINITE finite
#define CALLBACK __attribute__((__stdcall__))
#define DLLEXPORT
#define _inline
#endif 

double SMSDot(double a[],double b[],int n);

double SMSSum(double a[],int n);

double SMSDotB(double a[],double b[],int m, int n, int k);

double SMSSumB(double a[],int m, int n, int k);


double SMSKDelta(int i,int j);

double SMSDeltaPart(double *a,int i, int j,int k);

void SMSMove(double a[],double b[],int n);

void SMSZero(double a[],int n);

double Time();

#define Power(x, y) (pow((double)(x), (double)(y)))
#define Sqrt(x)        (sqrt((double)(x)))
#define Cbrt(x) ( x<0 ? -pow(-(double)(x),1./3.) : pow((double)(x),1./3.) )
#define Abs(x)        (fabs((double)(x)))
#define Less(a,b,c)        (a<b && b<c)
#define LessEqual(a,b,c)        (a<=b && b<=c)
#define Greater(a,b,c)        (a>b && b>c)
#define GreaterEqual(a,b,c)        (a>=b && b>=c)


/*====================================================================*/
/*MathLink specific*/

void SMSSetLinkOption(int i,int j);
void PutReal(double a);
void PutRealList(double *a,int n);
void exit_util();
void MathLinkInitialize();
void SMSLinkNoEvaluations();
int *MathLinkCallCount;

#define CO_SparseArray 0
#define CO_PauseOnExit 1
#define CO_NoSubroutines 2
#define CO_Length 3

int MathLinkOptions[CO_Length];

void PutRealArray(double **ap,double *as,int n,int m,int opt,int type);


/*Matrix functions definitions*/
void MatrixExponentialOMMF(double m[9], int *r, double mf[9]);
void MatrixExponentialOMMFS(double m[6], int *r, double mf[6]);
void MatrixExponentialOMM2D(double m[5], int *r, double mf[5]);
void MatrixExponentialOMM2DS(double m[4], int *r, double mf[4]);
void MatrixExponentialOMDMF(double m[9], int *r, double mf[9], double dmf[45]);
void MatrixExponentialOMDMFS(double m[6], int *r, double mf[6], double dmf[21]);
void MatrixExponentialOMDM2D(double m[5], int *r, double mf[5], double dmf[15]);
void MatrixExponentialOMDM2DS(double m[4], int *r, double mf[4], double dmf[10]);
void MatrixExponentialOMDDMF(double m[9], int *r, double mf[9], double dmf[45], double ddmf[165]);
void MatrixExponentialOMDDMFS(double m[6], int *r, double mf[6], double dmf[21], double ddmf[74]);
void MatrixExponentialOMDDM2D(double m[5], int *r, double mf[5], double dmf[15], double ddmf[35]);
void MatrixExponentialOMDDM2DS(double m[4], int *r, double mf[4], double dmf[10], double ddmf[21]);
void MatrixLogarithmOMMF(double m[9], int *r, double mf[9]);
void MatrixLogarithmOMMFS(double m[6], int *r, double mf[6]);
void MatrixLogarithmOMM2D(double m[5], int *r, double mf[5]);
void MatrixLogarithmOMM2DS(double m[4], int *r, double mf[4]);
void MatrixLogarithmOMDMF(double m[9], int *r, double mf[9], double dmf[45]);
void MatrixLogarithmOMDMFS(double m[6], int *r, double mf[6], double dmf[21]);
void MatrixLogarithmOMDM2D(double m[5], int *r, double mf[5], double dmf[15]);
void MatrixLogarithmOMDM2DS(double m[4], int *r, double mf[4], double dmf[10]);
void MatrixLogarithmOMDDMF(double m[9], int *r, double mf[9], double dmf[45], double ddmf[165]);
void MatrixLogarithmOMDDMFS(double m[6], int *r, double mf[6], double dmf[21], double ddmf[74]);
void MatrixLogarithmOMDDM2D(double m[5], int *r, double mf[5], double dmf[15], double ddmf[35]);
void MatrixLogarithmOMDDM2DS(double m[4], int *r, double mf[4], double dmf[10], double ddmf[21]);
void MatrixSquareRootOMMF(double m[9], int *r, double mf[9]);
void MatrixSquareRootOMMFS(double m[6], int *r, double mf[6]);
void MatrixSquareRootOMM2D(double m[5], int *r, double mf[5]);
void MatrixSquareRootOMM2DS(double m[4], int *r, double mf[4]);
void MatrixSquareRootOMDMF(double m[9], int *r, double mf[9], double dmf[45]);
void MatrixSquareRootOMDMFS(double m[6], int *r, double mf[6], double dmf[21]);
void MatrixSquareRootOMDM2D(double m[5], int *r, double mf[5], double dmf[15]);
void MatrixSquareRootOMDM2DS(double m[4], int *r, double mf[4], double dmf[10]);
void MatrixSquareRootOMDDMF(double m[9], int *r, double mf[9], double dmf[45], double ddmf[165]);
void MatrixSquareRootOMDDMFS(double m[6], int *r, double mf[6], double dmf[21], double ddmf[74]);
void MatrixSquareRootOMDDM2D(double m[5], int *r, double mf[5], double dmf[15], double ddmf[35]);
void MatrixSquareRootOMDDM2DS(double m[4], int *r, double mf[4], double dmf[10], double ddmf[21]);
void MatrixPowerOMMF(double m[9], double *p, int *r, double mf[9]);
void MatrixPowerOMMFS(double m[6], double *p, int *r, double mf[6]);
void MatrixPowerOMM2D(double m[5], double *p, int *r, double mf[5]);
void MatrixPowerOMM2DS(double m[4], double *p, int *r, double mf[4]);
void MatrixPowerOMDMF(double m[9], double *p, int *r, double mf[9], double dmf[45]);
void MatrixPowerOMDMFS(double m[6], double *p, int *r, double mf[6], double dmf[21]);
void MatrixPowerOMDM2D(double m[5], double *p, int *r, double mf[5], double dmf[15]);
void MatrixPowerOMDM2DS(double m[4], double *p, int *r, double mf[4], double dmf[10]);
void MatrixPowerOMDDMF(double m[9], double *p, int *r, double mf[9], double dmf[45], double ddmf[165]);
void MatrixPowerOMDDMFS(double m[6], double *p, int *r, double mf[6], double dmf[21], double ddmf[74]);
void MatrixPowerOMDDM2D(double m[5], double *p, int *r, double mf[5], double dmf[15], double ddmf[35]);
void MatrixPowerOMDDM2DS(double m[4], double *p, int *r, double mf[4], double dmf[10], double ddmf[21]);
void MatrixPowerSeriesOMMF(double m[9], double *cp, int *n, int *r, double mf[9]);
void MatrixPowerSeriesOMMFS(double m[6], double *cp, int *n, int *r, double mf[6]);
void MatrixPowerSeriesOMM2D(double m[5], double *cp, int *n, int *r, double mf[5]);
void MatrixPowerSeriesOMM2DS(double m[4], double *cp, int *n, int *r, double mf[4]);
void MatrixPowerSeriesOMDMF(double m[9], double *cp, int *n, int *r, double mf[9], double dmf[45]);
void MatrixPowerSeriesOMDMFS(double m[6], double *cp, int *n, int *r, double mf[6], double dmf[21]);
void MatrixPowerSeriesOMDM2D(double m[5], double *cp, int *n, int *r, double mf[5], double dmf[15]);
void MatrixPowerSeriesOMDM2DS(double m[4], double *cp, int *n, int *r, double mf[4], double dmf[10]);
void MatrixPowerSeriesOMDDMF(double m[9], double *cp, int *n, int *r, double mf[9], double dmf[45], double ddmf[165]);
void MatrixPowerSeriesOMDDMFS(double m[6], double *cp, int *n, int *r, double mf[6], double dmf[21], double ddmf[74]);
void MatrixPowerSeriesOMDDM2D(double m[5], double *cp, int *n, int *r, double mf[5], double dmf[15], double ddmf[35]);
void MatrixPowerSeriesOMDDM2DS(double m[4], double *cp, int *n, int *r, double mf[4], double dmf[10], double ddmf[21]);

void MatrixFunctionCompressMMF(double m[3][3],double mc[9]);
void MatrixFunctionUncompressMFMF(double mf[9],double mffull[3][3]);
void MatrixFunctionUncompressDMFMF(double dmf[45],double dmffull[3][3][3][3]);
void MatrixFunctionUncompressDDMFMF(double ddmf[165],double ddmffull[3][3][3][3][3][3]);
void MatrixFunctionCompressMMFS(double m[3][3],double mc[6]);
void MatrixFunctionUncompressMFMFS(double mf[6],double mffull[3][3]);
void MatrixFunctionUncompressDMFMFS(double dmf[21],double dmffull[3][3][3][3]);
void MatrixFunctionUncompressDDMFMFS(double ddmf[74],double ddmffull[3][3][3][3][3][3]);
void MatrixFunctionCompressMM2D(double m[3][3],double mc[5]);
void MatrixFunctionUncompressMFM2D(double mf[5],double mffull[3][3]);
void MatrixFunctionUncompressDMFM2D(double dmf[15],double dmffull[3][3][3][3]);
void MatrixFunctionUncompressDDMFM2D(double ddmf[35],double ddmffull[3][3][3][3][3][3]);
void MatrixFunctionCompressMM2DS(double m[3][3],double mc[4]);
void MatrixFunctionUncompressMFM2DS(double mf[4],double mffull[3][3]);
void MatrixFunctionUncompressDMFM2DS(double dmf[10],double dmffull[3][3][3][3]);
void MatrixFunctionUncompressDDMFM2DS(double ddmf[21],double ddmffull[3][3][3][3][3][3]);



#endif

