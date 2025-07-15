#include "sms.h"
#include "stdlib.h"
#include "stdio.h"
#include "mathlink.h"
#include "time.h"
#include "sys/timeb.h"

double SMSDot(double a[],double b[],int n)
 {double s=0e0;int i;
  for(i=0;i<n;i++)s+=a[i]*b[i];
  return(s);
 }

double SMSDotB(double a[],double b[], int m,int n,int k)
 {double s=0e0;int i;
  for(i=m-1;i<n;i=i+k)s+=a[i]*b[i];
  return(s);
 }

double SMSSum(double a[],int n)
 {double s=0e0;int i;
  for(i=0;i<n;i++)s+=a[i];
  return(s);
 }

double SMSSumB(double a[],int m,int n,int k)
 {double s=0e0;int i;
  for(i=m-1;i<n;i=i+k)s+=a[i];
  return(s);
 }

double SMSKDelta(int i,int j)
{if(i==j){return 1e0;}else{return 0e0;};}


double SMSDeltaPart(double *a,int i,int j,int k)
{div_t d=div(i,j);
  if(d.rem || d.quot>k){return 0e0;}else{return a[d.quot-1];};
}

void SMSMove(double a[],double b[],int n)
{int i;
 for(i=0;i<n;i++)b[i]=a[i];
}

void SMSZero(double a[],int n)
{int i;
 for(i=0;i<n;i++)a[i]=0e0;
}

/*====================================================================*/
/*MathLink specific*/

void exit_util()
{
 char tmp;
 printf("MathLink exit - press any key to exit the program");
 if(MathLinkOptions[CO_PauseOnExit])tmp=getchar();
}


void PutReal(double a){
 if(ISNAN(a) || !FINITE(a) ){
   MLPutSymbol(stdlink,"Indeterminate");
 }else{
   MLPutReal(stdlink,a);
 };
}

void PutRealList(double *a,int n){
 int i,j;
 j=0;
 for(i=0;i<n;i++)j=j || ISNAN(a[i]) || (! FINITE(a[i]));
 if(j){
  MLPutFunction(stdlink, "List",n);
  for(i=0;i<n;i++){
	  if(ISNAN(a[i]) || !FINITE(a[i])){
		  MLPutSymbol(stdlink,"Indeterminate");
    }else{
      MLPutReal(stdlink,a[i]);
    };
  };
 }else{
   MLPutRealList(stdlink,a,n);
 };
}

/*
opt - 0 Automatic
  - 1 Full array
  - 2 Sparse array
type 1 - **ap  pointer array
 	 0 - *as   static array

mt - 1 full 0 sparse
*/
void PutRealArray(double **ap,double *as,int n,int m,int opt,int type){
int i,j,ns,mt;
double eps;
double a;
eps=1.0e-15;
ns=0;
if(opt==0 || opt==2){
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      if(type){
        a=ap[i][j];
      }else{
        a=as[i*m+j];
      };
      if(fabs(a)>eps || ISNAN(a) || !FINITE(a) ){
        ++ns;
      };
    };
  };
  if(opt==0 && ns>0.25*n*m){
    mt=1;
  }else{
    mt=0;
  };
}else{
  mt=1;
};

if(mt){
  MLPutFunction(stdlink, "List",n);
  for(i=0;i<n;i++){
    if(type){
      PutRealList(ap[i],m);
    }else{
      PutRealList(&as[i*m],m);
    };
  };
}else{
  MLPutFunction(stdlink, "SparseArray",2);
  MLPutFunction(stdlink, "List",ns);
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      if(type){
        a=ap[i][j];
      }else{
        a=as[i*m+j];
      };
      if(ISNAN(a) || !FINITE(a)){
        MLPutFunction(stdlink, "Rule",2);
        MLPutFunction(stdlink, "List",2);
        MLPutInteger(stdlink,i+1);
        MLPutInteger(stdlink,j+1);
        MLPutSymbol(stdlink,"Indeterminate");
      }else if ( fabs(a)>eps ){
        MLPutFunction(stdlink, "Rule",2);
        MLPutFunction(stdlink, "List",2);
        MLPutInteger(stdlink,i+1);
        MLPutInteger(stdlink,j+1);
        MLPutReal(stdlink,a);
      };
    };
  };
  MLPutFunction(stdlink, "List",2);
  MLPutInteger(stdlink,n);
  MLPutInteger(stdlink,m);
};
}

/*
*/
double Time()
{  
   struct TIMEB timebuffer;
   char *timeline;
   double current;
   FTIME(&timebuffer );
   timeline = &(ctime( & ( timebuffer.time ) )[11]);
   current=timebuffer.millitm/1000.;
   current=current+(double)timebuffer.time;
   return current;
}

void SMSSetLinkOption(int i,int j)
{
  MathLinkOptions[i]=j;
  MLPutSymbol(stdlink,"True");
}

void SMSLinkNoEvaluations()
{
 MLPutIntegerList(stdlink,MathLinkCallCount,MathLinkOptions[CO_NoSubroutines]);
}


#ifdef _WINDOWS
int PASCAL WinMain( HANDLE hinstCurrent, HANDLE hinstPrevious, LPSTR lpszCmdLine, int nCmdShow)
{
	char  buff[512];
	char FAR * buff_start = buff;
	char FAR * argv[32];
	char FAR * FAR * argv_end = argv + 32;
	int i;
	for(i=0;i<CO_Length;i++){
	MathLinkOptions[i]=0;
	}; 
	MathLinkInitialize();
	MathLinkCallCount=(int *)calloc(MathLinkOptions[CO_NoSubroutines],sizeof(int));
	for(i=0;i<MathLinkOptions[CO_NoSubroutines];i++){
	MathLinkCallCount[i]=0;
	};
	if( !MLInitializeIcon( hinstCurrent, nCmdShow)) return 1;
	MLScanString( argv, &argv_end, &lpszCmdLine, &buff_start);
	return MLMain( argv_end - argv, argv);
}
#else
int main(int argc,char  *argv[]){
  int i;
  for(i=0;i<CO_Length;i++){
    MathLinkOptions[i]=0;
  }; 
  MathLinkInitialize();
  MathLinkCallCount=(int *)calloc(MathLinkOptions[CO_NoSubroutines],sizeof(int));
  for(i=0;i<MathLinkOptions[CO_NoSubroutines];i++){
    MathLinkCallCount[i]=0;
  };
  atexit(exit_util);

  return MLMain(argc, argv);
};
#endif