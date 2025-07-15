/*
 * This file automatically produced by /usr/local/Wolfram/Wolfram/14.1/SystemFiles/Links/MathLink/DeveloperKit/Linux-x86-64/CompilerAdditions/mprep from:
 *	/home/mr/job-related/MyOptimizationAlg/src/algorithms/SubproblemAceGenTest.tm
 * mprep Revision 19 Copyright (c) Wolfram Research, Inc. 1990-2024
 */

#define MPREP_REVISION 19

#include "mathlink.h"

#if defined(__cplusplus)
#define MLVOIDPARAM
#if __cplusplus >= 201103L
#define MLNULL nullptr
#endif
#endif

#if !defined(MLVOIDPARAM)
#define MLVOIDPARAM void
#endif
#if !defined(MLNULL)
#define MLNULL 0
#endif

int MLAbort = 0;
int MLDone  = 0;
long MLSpecialCharacter = '\0';

MLINK stdlink = MLNULL;
MLEnvironment stdenv = MLNULL;
MLYieldFunctionObject stdyielder = (MLYieldFunctionObject)MLNULL;
MLMessageHandlerObject stdhandler = (MLMessageHandlerObject)MLNULL;

/********************************* end header *********************************/


void SMSSetLinkOption ( int tp_1, int tp_2);

static int tr_0( MLINK mlp)
{
	int	res = 0;
	int tp_1;
	int tp_2;
	if ( ! MLGetInteger( mlp, &tp_1) ) goto L0;
	if ( ! MLGetInteger( mlp, &tp_2) ) goto L1;
	if ( ! MLNewPacket(mlp) ) goto L2;

	SMSSetLinkOption(tp_1, tp_2);

	res = 1;
L2: L1: 
L0:	return res;
} /* tr_0 */


void SMSLinkNoEvaluations ( MLVOIDPARAM);

static int tr_1( MLINK mlp)
{
	int	res = 0;
	if ( ! MLNewPacket(mlp) ) goto L0;
	if( !mlp) return res; /* avoid unused parameter warning */

	SMSLinkNoEvaluations();

	res = 1;

L0:	return res;
} /* tr_1 */


void SubproblemAceGenTestMathLink ( MLVOIDPARAM);

static int tr_2( MLINK mlp)
{
	int	res = 0;
	if( !mlp) return res; /* avoid unused parameter warning */

	SubproblemAceGenTestMathLink();

	res = 1;

	return res;
} /* tr_2 */


static struct func {
	int   f_nargs;
	int   manual;
	int   (*f_func)(MLINK);
	const char  *f_name;
	} tramps_[3] = {
		{ 2, 0, tr_0, "SMSSetLinkOption" },
		{ 0, 0, tr_1, "SMSLinkNoEvaluations" },
		{ 0, 2, tr_2, "SubproblemAceGenTestMathLink" }
		};

static const char* evalstrs[] = {
	"SMSLinkNoEvaluations[]:=SMSLinkNoEvaluations[SMSSessionName];",
	(const char*)MLNULL,
	"SMSSetLinkOptions[i__Rule]:=SMSSetLinkOptions[SMSSessionName,i];",
	(const char*)MLNULL,
	"SMSSetLinkOptions[s_String,Rule[i_String,j_]]:=SMSSetLinkOption[",
	"s,{i,j}/.{{\"SparseArray\",True}->{0,2},{\"SparseArray\",False}->{0,",
	"1},{\"SparseArray\",Automatic}->{0,0},{\"PauseOnExit\",True}->{1,1},",
	"{\"PauseOnExit\",False}->{1,0},_:>(Print[\"Incorrect option: \",Rule",
	"[i,j]];Abort[])}];",
	(const char*)MLNULL,
	"SMSSetLinkOptions[s_String,i__Rule]:=Map[SMSSetLinkOptions[s,#]&",
	",{i}];",
	(const char*)MLNULL,
	(const char*)MLNULL
};
#define CARDOF_EVALSTRS 4

static int definepattern_( MLINK, char*, char*, int);

static int doevalstr_( MLINK, int);

int  MLDoCallPacket_( MLINK, struct func[], int);


int MLInstall( MLINK mlp)
{
	int _res;
	_res = MLConnect(mlp);
	if (_res) _res = definepattern_(mlp, (char *)"SMSSetLinkOption[\"SubproblemAceGenTest\",{i_Integer,j_Integer}]", (char *)"{i,j}", 0);
	if (_res) _res = definepattern_(mlp, (char *)"SMSLinkNoEvaluations[\"SubproblemAceGenTest\"]", (char *)"{}", 1);
	if (_res) _res = doevalstr_( mlp, 0);
	if (_res) _res = doevalstr_( mlp, 1);
	if (_res) _res = doevalstr_( mlp, 2);
	if (_res) _res = doevalstr_( mlp, 3);
	if (_res) _res = definepattern_(mlp, (char *)"SubproblemAceGenTest[solution_?(ArrayQ[#,1,(Head[#]==Real || Head[#]==Integer&)] && Dimensions[#]==={2} &)]", (char *)"{solution}", 2);
	if (_res) _res = MLPutSymbol( mlp, "End");
	if (_res) _res = MLFlush( mlp);
	return _res;
} /* MLInstall */


int MLDoCallPacket( MLINK mlp)
{
	return MLDoCallPacket_( mlp, tramps_, 3);
} /* MLDoCallPacket */

/******************************* begin trailer ********************************/

#ifndef EVALSTRS_AS_BYTESTRINGS
#	define EVALSTRS_AS_BYTESTRINGS 1
#endif


#if CARDOF_EVALSTRS
static int  doevalstr_( MLINK mlp, int n)
{
	long bytesleft, bytesnow;
#if !EVALSTRS_AS_BYTESTRINGS
	long charsleft, charsnow;
#endif
	char **s, **p;
	char *t;

	s = (char **)evalstrs;
	while( n-- > 0){
		if( *s == MLNULL) break;
		while( *s++ != MLNULL){}
	}
	if( *s == MLNULL) return 0;
	bytesleft = 0;
#if !EVALSTRS_AS_BYTESTRINGS
	charsleft = 0;
#endif
	p = s;
	while( *p){
		t = *p; while( *t) ++t;
		bytesnow = t - *p;
		bytesleft += bytesnow;
#if !EVALSTRS_AS_BYTESTRINGS
		charsleft += bytesnow;
		t = *p;
		charsleft -= MLCharacterOffset( &t, t + bytesnow, bytesnow);
		/* assert( t == *p + bytesnow); */
#endif
		++p;
	}


	MLPutNext( mlp, MLTKSTR);
#if EVALSTRS_AS_BYTESTRINGS
	p = s;
	while( *p){
		t = *p; while( *t) ++t;
		bytesnow = t - *p;
		bytesleft -= bytesnow;
		MLPut8BitCharacters( mlp, bytesleft, (unsigned char*)*p, bytesnow);
		++p;
	}
#else
	MLPut7BitCount( mlp, charsleft, bytesleft);
	p = s;
	while( *p){
		t = *p; while( *t) ++t;
		bytesnow = t - *p;
		bytesleft -= bytesnow;
		t = *p;
		charsnow = bytesnow - MLCharacterOffset( &t, t + bytesnow, bytesnow);
		/* assert( t == *p + bytesnow); */
		charsleft -= charsnow;
		MLPut7BitCharacters(  mlp, charsleft, *p, bytesnow, charsnow);
		++p;
	}
#endif
	return MLError( mlp) == MLEOK;
}
#endif /* CARDOF_EVALSTRS */


static int  definepattern_( MLINK mlp, char *patt, char *args, int func_n)
{
	MLPutFunction( mlp, "DefineExternal", (long)3);
	  MLPutString( mlp, patt);
	  MLPutString( mlp, args);
	  MLPutInteger( mlp, func_n);
	return !MLError(mlp);
} /* definepattern_ */


int MLDoCallPacket_( MLINK mlp, struct func functable[], int nfuncs)
{
	int len;
	int n, res = 0;
	struct func* funcp;

	if( ! MLGetInteger( mlp, &n) ||  n < 0 ||  n >= nfuncs) goto L0;
	funcp = &functable[n];

	if( funcp->f_nargs >= 0
	&& ( ! MLTestHead(mlp, "List", &len)
	     || ( !funcp->manual && (len != funcp->f_nargs))
	     || (  funcp->manual && (len <  funcp->f_nargs))
	   )
	) goto L0;

	stdlink = mlp;
	res = (*funcp->f_func)( mlp);

L0:	if( res == 0)
		res = MLClearError( mlp) && MLPutSymbol( mlp, "$Failed");
	return res && MLEndPacket( mlp) && MLNewPacket( mlp);
} /* MLDoCallPacket_ */


mlapi_packet MLAnswer( MLINK mlp)
{
	mlapi_packet pkt = 0;
	int waitResult;

	while( ! MLDone && ! MLError(mlp)
		&& (waitResult = MLWaitForLinkActivity(mlp),waitResult) &&
		waitResult == MLWAITSUCCESS && (pkt = MLNextPacket(mlp), pkt) &&
		pkt == CALLPKT)
	{
		MLAbort = 0;
		if(! MLDoCallPacket(mlp))
			pkt = 0;
	}
	MLAbort = 0;
	return pkt;
} /* MLAnswer */



/*
	Module[ { me = $ParentLink},
		$ParentLink = contents of RESUMEPKT;
		Message[ MessageName[$ParentLink, "notfe"], me];
		me]
*/

static int refuse_to_be_a_frontend( MLINK mlp)
{
	int pkt;

	MLPutFunction( mlp, "EvaluatePacket", 1);
	  MLPutFunction( mlp, "Module", 2);
	    MLPutFunction( mlp, "List", 1);
		  MLPutFunction( mlp, "Set", 2);
		    MLPutSymbol( mlp, "me");
	        MLPutSymbol( mlp, "$ParentLink");
	  MLPutFunction( mlp, "CompoundExpression", 3);
	    MLPutFunction( mlp, "Set", 2);
	      MLPutSymbol( mlp, "$ParentLink");
	      MLTransferExpression( mlp, mlp);
	    MLPutFunction( mlp, "Message", 2);
	      MLPutFunction( mlp, "MessageName", 2);
	        MLPutSymbol( mlp, "$ParentLink");
	        MLPutString( mlp, "notfe");
	      MLPutSymbol( mlp, "me");
	    MLPutSymbol( mlp, "me");
	MLEndPacket( mlp);

	while( (pkt = MLNextPacket( mlp), pkt) && pkt != SUSPENDPKT)
		MLNewPacket( mlp);
	MLNewPacket( mlp);
	return MLError( mlp) == MLEOK;
}


int MLEvaluate( MLINK mlp, char *s)
{
	if( MLAbort) return 0;
	return MLPutFunction( mlp, "EvaluatePacket", 1L)
		&& MLPutFunction( mlp, "ToExpression", 1L)
		&& MLPutString( mlp, s)
		&& MLEndPacket( mlp);
} /* MLEvaluate */


int MLEvaluateString( MLINK mlp, char *s)
{
	int pkt;
	if( MLAbort) return 0;
	if( MLEvaluate( mlp, s)){
		while( (pkt = MLAnswer( mlp), pkt) && pkt != RETURNPKT)
			MLNewPacket( mlp);
		MLNewPacket( mlp);
	}
	return MLError( mlp) == MLEOK;
} /* MLEvaluateString */


void MLDefaultHandler( MLINK mlp, int message, int n)
{
	switch (message){
	case MLTerminateMessage:
		MLDone = 1;
	case MLInterruptMessage:
	case MLAbortMessage:
		MLAbort = 1;
	default:
		return;
	}
}

static int MLMain_( char **argv, char **argv_end, char *commandline)
{
	MLINK mlp;
	int err;

	if( !stdenv)
		stdenv = MLInitialize( (MLEnvironmentParameter)MLNULL);

	if( stdenv == (MLEnvironment)MLNULL) goto R0;

	if( !stdhandler)
		stdhandler = (MLMessageHandlerObject)MLDefaultHandler;


	mlp = commandline
		? MLOpenString( stdenv, commandline, &err)
		: MLOpenArgcArgv( stdenv, (int)(argv_end - argv), argv, &err);
	if( mlp == (MLINK)MLNULL){
		MLAlert( stdenv, MLErrorString( stdenv, err));
		goto R1;
	}

	if( stdyielder) MLSetYieldFunction( mlp, stdyielder);
	if( stdhandler) MLSetMessageHandler( mlp, stdhandler);

	if( MLInstall( mlp))
		while( MLAnswer( mlp) == RESUMEPKT){
			if( ! refuse_to_be_a_frontend( mlp)) break;
		}

	MLClose( mlp);
R1:	MLDeinitialize( stdenv);
	stdenv = (MLEnvironment)MLNULL;
R0:	return !MLDone;
} /* MLMain_ */


int MLMainString( char *commandline)
{
	return MLMain_( (charpp_ct)MLNULL, (charpp_ct)MLNULL, commandline);
}

int MLMainArgv( char** argv, char** argv_end) /* note not FAR pointers */
{   
	static char * far_argv[128];
	int count = 0;
	
	while(argv < argv_end)
		far_argv[count++] = *argv++;
		 
	return MLMain_( far_argv, far_argv + count, (charp_ct)MLNULL);

}

int MLMain( int argc, char **argv)
{
 	return MLMain_( argv, argv + argc, (char *)MLNULL);
}
 
