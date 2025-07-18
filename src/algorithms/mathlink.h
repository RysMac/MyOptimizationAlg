/*************************************************************************

        Copyright 1986 through 2024 by Wolfram Research Inc.
        All rights reserved

*************************************************************************/

#ifndef _MATHLINK_H
#define _MATHLINK_H

#ifndef _MLVERS_H
#define _MLVERS_H

#ifndef _MLPLATFM_H
#define _MLPLATFM_H


#if !WINDOWS_MATHLINK && !WINDOWS_WSTP && !UNIX_MATHLINK && !UNIX_WSTP
#	define	UNIX_MATHLINK	1

#endif

#if defined(__GNUC__)
#if defined(__GNUC_PATCHLEVEL__)
#define GCC_MATHLINK_VERSION (__GNUC__ * 10000 \
						+ __GNUC_MINOR__ * 100 \
						+ __GNUC_PATCHLEVEL__)
#else
#define GCC_MATHLINK_VERSION (__GNUC__ * 10000 \
						+ __GNUC_MINOR__ * 100)
#endif
#endif

#if WINDOWS_MATHLINK
#define LITTLEENDIAN_NUMERIC_TYPES 1
#if defined(WIN64) || defined(_WIN64)
#define WIN64_MATHLINK 1
#if( _M_X64 || __x86_64 || __x86_64__ || x86_64)
#define X64_WINDOWS_MATHLINK 1
#elif (_M_ARM64 || __arm64__)
#define ARM64_WINDOWS_MATHLINK 1
#else
#error
#endif
#elif defined(WIN32) || defined(_WIN32)
#define WIN32_MATHLINK 1
#if( _M_IX86 || __i386 || __i386__ || i386)
#define I86_WINDOWS_MATHLINK 1
#else
#error
#endif
#else
#error "WINDOWS_MATHLINK is true but neither WIN32 nor WIN64 are defined."
#endif
#elif UNIX_MATHLINK
#if (__MACH || __MACH__ || MACH) && !defined(DARWIN_MATHLINK)
#define LITTLEENDIAN_NUMERIC_TYPES 1
#define DARWIN_MATHLINK 1
#if __i386 || __i386__ || i386
#define INTEL_DARWIN_MATHLINK 1
#define X86_DARWIN_MATHLINK 1
#elif __x86_64 || __x86_64__ || x86_64
#define INTEL_DARWIN_MATHLINK 1
#define X86_64_DARWIN_MATHLINK 1
#elif __arm__
#define ARM_DARWIN_MATHLINK 1
#elif __arm64__
#define ARM64_DARWIN_MATHLINK 1
#else
			#error not yet implemented
#endif

#if __clang__
#define CLANG_MATHLINK 1
#endif
#elif (__linux || __linux__ || linux) && !defined(LINUX_MATHLINK)
#define LITTLEENDIAN_NUMERIC_TYPES 1
#define LINUX_MATHLINK 1
#if __x86_64 || __x86_64__ || x86_64
#define X86_64_LINUX_MATHLINK 1
#elif __i386 || __i386__ || i386
#define I86_LINUX_MATHLINK 1
#elif __ANDROID || __ANDROID__ || ANDROID
#define ANDROID_LINUX_MATHLINK 1
#if __arm || __arm__ || arm
#define ANDROID_ARM_LINUX_MATHLINK 1
#endif
#elif __arm || __arm__ || arm
#define ARM_LINUX_MATHLINK 1
#elif __ARM_ARCH == 8
#define ARM64_LINUX_MATHLINK 1
#else
			not yet implemented
#endif
#else
#error
#endif
#endif



#if !defined(LITTLEENDIAN_NUMERIC_TYPES) && !defined(BIGENDIAN_NUMERIC_TYPES)
#error not yet implemented
#endif


#endif /* _MLPLATFM_H */


#ifndef MLVERSION
#define MLVERSION 6
#endif
#ifndef MLINTERFACE
#define MLINTERFACE 4
#endif

#if MLVERSION != 6 || MLINTERFACE != 4
#error "This header file only works with MLVERSION=6, MLINTERFACE=4"
#endif

#define MLREVISION 51
#define MLMATHVERSION 14.1.0


/*
 * MathLink adopts a simple versioning strategy that can be adapted to many
 * compile-time and run-time environments.  In particular, it is amenable to
 * the various shared library facilities in use.  (Although certain of these
 * facilities provide more sophisticated mechanisms than are required by the
 * following simple strategy.)
 *
 * MathLink evolves by improving its implementation and by improving its
 * interface.  The values of MLREVISION or MLINTERFACE defined here are
 * incremented whenever an improvement is made and released.
 *
 * MLREVISION is the current revision number. It is incremented every time
 * a change is made to the source and MathLink is rebuilt and distributed
 * on any platform.  (Bug fixes, optimizations, or other improvements
 * transparent to the interface increment only this number.)
 *
 * MLINTERFACE is a name for a documented interface to MathLink.  This
 * number is incremented whenever a named constant or function is added,
 * removed, or its behavior is changed in a way that could break existing
 * correct* client programs.  It is expected that the interface to MathLink
 * is improved over time so that implemenations of higher numbered
 * interfaces are more complete or more convenient to use for writing
 * effective client programs.  In particular, a specific interface provides
 * all the useful functionality of an earlier interface.
 *
 *     *(It is possible that an incorrect MathLink program still works
 *     because it relies on some undocumented detail of a particular
 *     revision.  It may not always be possible to change the interface
 *     number when such a detail changes.  For example, one program may
 *     be relying on a bug in MathLink that a great many other programs
 *     need fixed.  In this case, we would likely choose to potentially
 *     break the incorrect program in order to fix the correct programs
 *     by incrementing the revision number leaving the interface number
 *     unchanged.  It is possible to bind to a particular revision of a
 *     MathLink interface if that is important for some programs.  One
 *     could use a statically linked version of the library, make use of
 *     the search algorithm used by the runtime loader, or dynamically
 *     load the MathLink library manually.)
 *
 *
 * If a distributed MathLink implmentation were labeled with its revision
 * and interface numbers in dotted notation so that, say, ML.1.6 means the
 * sixth revision of interface one, then the following may represent the
 * distribution history of MathLink.
 *
 *     first distribution
 *         ML.1.5   (Perhaps earlier revisions were never
 *                   distributed for this platform.)
 *
 *     second distribution
 *         ML.1.6   (Bug fixes or other improvements were
 *                   made that don't affect the interface.)
 *
 *     third distribution
 *         ML.2.7   (Perhaps some new functions were added.)
 *
 *         ML.1.7   (And improvements were made that don't
 *                   affect the old interface.)
 *
 *     fourth distribution
 *         ML.3.8   (Perhaps the return values of an existing
 *                   function changed.)
 *         ML.2.8   (Revision 8 also adds improvements transparent
 *                   to interface 2.)
 *         ML.1.7   (Clients of interface 1 see no improvements
 *                   in this eighth revision.)
 *
 * Note that the distribution history may not be the same on different
 * platforms.  But revision numbers represent a named body of source code
 * across all platforms.
 *
 * The mechanism for deploying this strategy differs between platforms
 * because of differing platform-specific facilities and conventions.
 * The interface and revision numbers may form part of the filename of
 * the MathLink library, or they may not.  This information is always
 * available in some conventional form so that it is easy and natural for
 * client programs to bind with and use the best available implementation
 * of a particular MathLink interface.  The details are described in the
 * MathLink Developer's Guide for each platform.
 */


#ifndef MLCREATIONID
#define MLCREATIONID 114411
#endif

#define MLAPI1REVISION 1   /* the first revision to support interface 1 */
#define MLAPI2REVISION 6   /* the first revision to support interface 2 */
#define MLAPI3REVISION 16  /* the first revision to support interface 3 */
#define MLAPI4REVISION 25  /* the first revision to support interface 4 */


#define MLAPIREVISION MLAPI4REVISION

		/*
		 * Interface 4 adds the following exported functions:
         *
         * MLCreateLinkWithExternalProtocol
         * MLDoNotHandleSignalParameter
         * MLEnableLinkLock
         * MLFilterArgv
         * MLGetAvailableLinkProtocolNames
         * MLGetInteger8
         * MLGetInteger8Array
         * MLGetInteger8ArrayData
         * MLGetInteger8List
         * MLGetLinksFromEnvironment
         * MLGetNumberAsByteString
         * MLGetNumberAsString
         * MLGetNumberAsUCS2String
         * MLGetNumberAsUTF16String
         * MLGetNumberAsUTF32String
         * MLGetNumberAsUTF8String
         * MLHandleSignal
         * MLIsLinkLoopback
         * MLLinkEnvironment
         * MLLogFileNameForLink
         * MLLogStreamToFile
         * MLLowLevelDeviceName
         * MLPutInteger8
         * MLPutInteger8Array
         * MLPutInteger8ArrayData
         * MLPutInteger8List
         * MLPutRealNumberAsUCS2String
         * MLPutRealNumberAsUTF16String
         * MLPutRealNumberAsUTF32String
         * MLPutRealNumberAsUTF8String
         * MLReleaseInteger8Array
         * MLReleaseInteger8List
         * MLReleaseLinkProtocolNames
         * MLReleaseLinksFromEnvironment
         * MLReleaseLogFileNameForLink
         * MLReleaseLowLevelDeviceName
         * MLReleaseParameters
         * MLSetThreadSafeLinksParameter
         * MLStopHandlingSignal
         * MLStopLoggingStream
         * MLValid
         * MLWaitForLinkActivity
         * MLWaitForLinkActivityWithCallback
         *
         * Interface 4 removes the following API functions.
         *
         * MLCreate0
         * MLDestroy
         * MLDoNotHandleSignalParameter0
         * MLFilterArgv0
         * MLGetByteString0
         * MLGetLinkedEnvIDString0
         * MLGetString0
         * MLGetUCS2String0
         * MLGetUTF16String0
         * MLGetUTF32String0
         * MLGetUTF8String0
         * MLGetUnicodeString0
         * MLHandleSignal0
         * MLLinkSelect
         * MLLoopbackOpen0
         * MLMake
         * MLPutRealByteString0
         * MLPutRealUCS2String0
         * MLPutRealUTF16String0
         * MLPutRealUTF32String0
         * MLPutRealUTF8String0
         * MLPutRealUnicodeString0
         * MLSetEnvIDString0
         * MLSetMessageHandler0
         * MLSetSignalHandler0
         * MLSetYieldFunction0
         * MLStopHandlingSignal0
         * MLUnsetSignalHandler0
         * MLValid0
         * MLVersionNumber0
         *
		 */


                /*
                 * Interface 3 adds the following exported functions:
                 *      MLClearAllSymbolReplacements
                 *      MLClearSymbolReplacement
                 *      MLConvertUCS2String
                 *      MLConvertUCS2StringNL
                 *      MLConvertUTF8String
                 *      MLConvertUTF8StringNL
                 *      MLConvertUTF16String
                 *      MLConvertUTF16StringNL
                 *      MLConvertUTF32String
                 *      MLConvertUTF32StringNL
                 *      MLEnvironmentData
                 *      MLGetDomainNameList
                 *      MLGetInteger16
                 *      MLGetInteger16Array
                 *      MLGetInteger16ArrayData
                 *      MLGetInteger16List
                 *      MLGetInteger32
                 *      MLGetInteger32Array
                 *      MLGetInteger32ArrayData
                 *      MLGetInteger32List
                 *      MLGetInteger64
                 *      MLGetInteger64Array
                 *      MLGetInteger64ArrayData
                 *      MLGetInteger64List
                 *      MLGetLinkedEnvIDString
                 *      MLGetMessageHandler
                 *      MLGetNetworkAddressList
                 *      MLGetReal128
                 *      MLGetReal128Array
                 *      MLGetReal128ArrayData
                 *      MLGetReal128List
                 *      MLGetReal32
                 *      MLGetReal32Array
                 *      MLGetReal32ArrayData
                 *      MLGetReal32List
                 *      MLGetReal64
                 *      MLGetReal64Array
                 *      MLGetReal64ArrayData
                 *      MLGetReal64List
                 *      MLGetUCS2Characters
                 *      MLGetUCS2String
                 *      MLGetUCS2Symbol
                 *      MLGetUTF16Characters
                 *      MLGetUTF16String
                 *      MLGetUTF16Symbol
                 *      MLGetUTF32Characters
                 *      MLGetUTF32String
                 *      MLGetUTF32Symbol
                 *      MLGetUTF8Characters
                 *      MLGetUTF8String
                 *      MLGetUTF8Symbol
                 *      MLGetYieldFunction
                 *      MLLinkName
                 *      MLOldConvertUCS2String
                 *      MLOpenArgcArgv
                 *      MLPutInteger16
                 *      MLPutInteger16Array
                 *      MLPutInteger16ArrayData
                 *      MLPutInteger16List
                 *      MLPutInteger32
                 *      MLPutInteger32Array
                 *      MLPutInteger32ArrayData
                 *      MLPutInteger32List
                 *      MLPutInteger64
                 *      MLPutInteger64Array
                 *      MLPutInteger64ArrayData
                 *      MLPutInteger64List
                 *      MLPutMessageWithArg
                 *      MLPutReal128
                 *      MLPutReal128Array
                 *      MLPutReal128ArrayData
                 *      MLPutReal128List
                 *      MLPutReal32
                 *      MLPutReal32Array
                 *      MLPutReal32ArrayData
                 *      MLPutReal32List
                 *      MLPutReal64
                 *      MLPutReal64Array
                 *      MLPutReal64ArrayData
                 *      MLPutReal64List
                 *      MLPutUCS2Characters
                 *      MLPutUCS2String
                 *      MLPutUCS2Symbol
                 *      MLPutUTF16Characters
                 *      MLPutUTF16String
                 *      MLPutUTF16Symbol
                 *      MLPutUTF32Characters
                 *      MLPutUTF32String
                 *      MLPutUTF32Symbol
                 *      MLPutUTF8Characters
                 *      MLPutUTF8String
                 *      MLPutUTF8Symbol
                 *      MLReadyParallel
                 *      MLReleaseBinaryNumberArray
                 *      MLReleaseByteArray
                 *      MLReleaseByteString
                 *      MLReleaseByteSymbol
                 *      MLReleaseDomainNameList
                 *      MLReleaseInteger16Array
                 *      MLReleaseInteger16List
                 *      MLReleaseInteger32Array
                 *      MLReleaseInteger32List
                 *      MLReleaseInteger64Array
                 *      MLReleaseInteger64List
                 *      MLReleaseNetworkAddressList
                 *      MLReleaseReal128Array
                 *      MLReleaseReal128List
                 *      MLReleaseReal32Array
                 *      MLReleaseReal32List
                 *      MLReleaseReal64Array
                 *      MLReleaseReal64List
                 *      MLReleaseString
                 *      MLReleaseSymbol
                 *      MLReleaseUCS2String
                 *      MLReleaseUCS2Symbol
                 *      MLReleaseUTF16String
                 *      MLReleaseUTF16Symbol
                 *      MLReleaseUTF32String
                 *      MLReleaseUTF32Symbol
                 *      MLReleaseUTF8String
                 *      MLReleaseUTF8Symbol
                 *      MLSetEncodingParameter
                 *      MLSetEnvIDString
                 *      MLSetEnvironmentData
                 *      MLSetSignalHandler
                 *      MLSetSignalHandlerFromFunction
                 *      MLSetSymbolReplacement
                 *      MLTestHead
                 *      MLUnsetSignalHandler
                 *      MLVersionNumbers
                 *
                 *      Interface 3 removes all the special MathLink types such as kcharp_ct, long_st, etc.
                 *      The API functions now all take standard C data types.
                 */

       		/*
                 * Interface 2 adds the following exported functions:
                 *      MLGetBinaryNumberArray0
                 *      MLTransfer0
                 *      MLNextCharacter0
                 * And, for WINDOWS_MATHLINK, some constants in "mlntypes.h"
                 * were changed in a way that causes MLGetRawType to return
                 * different values.
                 *
                 */


/* It may be possible for an implementation of one MathLink interface to
 * fully support an earlier interface.  MLNewParameters() may succeed when
 * passed an interface number less than the value of MLAPIREVISION when the
 * library was built.  This would happen, if the newer interface is a proper
 * superset of the older interface, or if the implementation can adjust its
 * behavior at runtime to conform to the older requested interface.
 */

#if WINDOWS_MATHLINK
#define MLOLDDEFINITION MLAPI3REVISION
#else
#define MLOLDDEFINITION MLAPI1REVISION
#endif


#endif /* _MLVERS_H */



#ifndef ML_EXTERN_C

#if defined(__cplusplus)
#define ML_C "C"
#define ML_EXTERN_C extern "C" {
#define ML_END_EXTERN_C }
#define WS_C "C"
#define WS_EXTERN_C extern "C" {
#define WS_END_EXTERN_C }
#else
#define ML_C
#define ML_EXTERN_C
#define ML_END_EXTERN_C
#define WS_C
#define WS_EXTERN_C
#define WS_END_EXTERN_C
#endif

#endif






#ifndef _MLDEVICE_H
#define _MLDEVICE_H




#ifndef _MLTYPES_H
#define _MLTYPES_H




#ifndef _MLBASICTYPES_H
#define _MLBASICTYPES_H




#ifndef _MLINT64_H
#define _MLINT64_H


#if WIN64_MATHLINK || X86_64_LINUX_MATHLINK || ARM64_LINUX_MATHLINK || ARM64_DARWIN_MATHLINK || X86_64_DARWIN_MATHLINK
#define ML64BIT_MATHLINK 1
#endif

#endif /* MLINT64_H */



#if ML64BIT_MATHLINK
#define ML_SMALLEST_SIGNED_64BIT        -9223372036854775807L - 1
#define ML_LARGEST_SIGNED_64BIT         9223372036854775807L
#define ML_LARGEST_UNSIGNED_64BIT       18446744073709551615UL
#else
#define ML_SMALLEST_SIGNED_64BIT        -9223372036854775807LL - 1
#define ML_LARGEST_SIGNED_64BIT         9223372036854775807LL
#define ML_LARGEST_UNSIGNED_64BIT       18446744073709551615ULL
#endif /* ML64BIT_MATHLINK */

#ifndef NO_INT64_STRUCT
#define NO_INT64_STRUCT

#if ML64BIT_MATHLINK
typedef unsigned int wint;
#else
typedef unsigned long wint;
#endif /* ML64BIT_MATHLINK */

#if LITTLEENDIAN_NUMERIC_TYPES
typedef struct _wint{
	wint low, hi;
} wint64;
#else
typedef struct _wint{
	wint hi, low;
} wint64;
#endif /* LITTLEENDIAN_NUMERIC_TYPES */


#endif /* NO_INT64_STRUCT */


#if WINDOWS_MATHLINK

typedef long mllong32;
typedef long wslong32;
typedef unsigned long mlulong32;
typedef unsigned long wsulong32;


typedef __int64 mlint64;
typedef __int64 wsint64;
typedef unsigned __int64 mluint64;
typedef unsigned __int64 wsuint64;

#if WIN64_MATHLINK
typedef mlint64 mlbigint;
typedef wsint64 wsbigint;
typedef mluint64 mlbiguint;
typedef wsuint64 wsbiguint;
#else
typedef mllong32 mlbigint;
typedef wslong32 wsbigint;
typedef mlulong32 mlbiguint;
typedef wsulong32 wsbiguint;
#endif


#elif X86_64_LINUX_MATHLINK || ARM64_LINUX_MATHLINK || ARM64_DARWIN_MATHLINK || X86_64_DARWIN_MATHLINK

typedef int mllong32;
typedef int wslong32;
typedef unsigned int mlulong32;
typedef unsigned int wsulong32;


typedef long mlint64;
typedef long wsint64;
typedef unsigned long mluint64;
typedef unsigned long wsuint64;

typedef mlint64 mlbigint;
typedef mluint64 mlbiguint;

#else /* All other 32 bit platforms */

typedef long mllong32;
typedef long wslong32;
typedef unsigned long mlulong32;
typedef unsigned long wsulong32;


typedef long long mlint64;
typedef long long wsint64;
typedef unsigned long long mluint64;
typedef unsigned long long wsuint64;

typedef mllong32 mlbigint;
typedef mlulong32 mlbiguint;

#endif /* WINDOWS_MATHLINK */









#endif /* _MLBASICTYPES_H */



#if WINDOWS_MATHLINK
#define MLCB __stdcall
#define WSCB __stdcall
#define MLAPI __stdcall
#define WSAPI __stdcall

#else
#define MLCB
#define WSCB
#define MLAPI
#define WSAPI
#endif

#if defined(LINUX_MATHLINK)
#define MLATTR __attribute__ ((visibility("default")))
#define WSATTR __attribute__ ((visibility("default")))
#else
#define MLATTR
#define WSATTR
#endif

#define MLAPI_ MLAPI


#ifndef MLDEFN
#define MLDEFN( rtype, name, params) extern MLATTR rtype MLAPI name params
#endif
#ifndef MLDECL
#define MLDECL( rtype, name, params) extern rtype MLAPI name params
#endif

#ifndef WSDEFN
#define WSDEFN( rtype, name, params) extern WSATTR rtype WSAPI name params
#endif
#ifndef WSDECL
#define WSDECL( rtype, name, params) extern rtype WSAPI name params
#endif

#ifndef ML_DEFN
#define ML_DEFN( rtype, name, params) extern rtype MLAPI_ name params
#endif
#ifndef ML_DECL
#define ML_DECL( rtype, name, params) extern ML_C rtype MLAPI_ name params
#endif



#ifndef MLCBPROC
#define MLCBPROC( rtype, name, params) typedef rtype (MLCB * name) params
#endif
#ifndef MLCBDECL
#define MLCBDECL( rtype, name, params) extern rtype MLCB name params
#endif
#ifndef MLCBDEFN
#define MLCBDEFN( rtype, name, params) extern rtype MLCB name params
#endif

#ifndef WSCBPROC
#define WSCBPROC( rtype, name, params) typedef rtype (WSCB * name) params
#endif
#ifndef WSCBDECL
#define WSCBDECL( rtype, name, params) extern rtype WSCB name params
#endif
#ifndef WSCBDEFN
#define WSCBDEFN( rtype, name, params) extern rtype WSCB name params
#endif




/* move into mlalert.h */
#ifndef MLDPROC
#define MLDPROC MLCBPROC
#endif
#ifndef WSDPROC
#define WSDPROC WSCBPROC
#endif
#ifndef MLDDECL
#define MLDDECL MLCBDECL
#endif
#ifndef WSDDECL
#define WSDDECL WSCBDECL
#endif
#ifndef MLDDEFN
#define MLDDEFN MLCBDEFN
#endif
#ifndef WSDDEFN
#define WSDDEFN WSCBDEFN
#endif




/* move into ml3state.h or mlstrenv.h */
#ifndef MLTPROC
#define MLTPROC MLCBPROC
#endif
#ifndef WSTPROC
#define WSTPROC WSCBPROC
#endif
#ifndef MLTDECL
#define MLTDECL MLCBDECL
#endif
#ifndef WSTDECL
#define WSTDECL WSCBDECL
#endif
#ifndef MLTDEFN
#define MLTDEFN MLCBDEFN
#endif
#ifndef WSTDEFN
#define WSTDEFN WSCBDEFN
#endif


/* move into mlnumenv.h */
#ifndef MLNPROC
#define MLNPROC MLCBPROC
#endif
#ifndef WSNPROC
#define WSNPROC WSCBPROC
#endif
#ifndef MLNDECL
#define MLNDECL MLCBDECL
#endif
#ifndef WSNDECL
#define WSNDECL WSCBDECL
#endif
#ifndef MLNDEFN
#define MLNDEFN MLCBDEFN
#endif
#ifndef WSNDEFN
#define WSNDEFN WSCBDEFN
#endif


/* move into mlalloc.h */
#ifndef MLAPROC
#define MLAPROC MLCBPROC
#endif
#ifndef MLADECL
#define MLADECL MLCBDECL
#endif
#ifndef MLADEFN
#define MLADEFN MLCBDEFN
#endif
#ifndef MLFPROC
#define MLFPROC MLCBPROC
#endif
#ifndef MLFDECL
#define MLFDECL MLCBDECL
#endif
#ifndef MLFDEFN
#define MLFDEFN MLCBDEFN
#endif
#ifndef MLRAPROC
#define MLRAPROC MLCBPROC
#endif
#ifndef MLRADECL
#define MLRADECL MLCBDECL
#endif
#ifndef MLRADEFN
#define MLRADEFN MLCBDEFN
#endif


/* move into mlstddev.h */
#ifndef MLYPROC
#define MLYPROC MLCBPROC
#endif
#ifndef WSYPROC
#define WSYPROC WSCBPROC
#endif
#ifndef MLYDECL
#define MLYDECL MLCBDECL
#endif
#ifndef WSYDECL
#define WSYDECL WSCBDECL
#endif
#ifndef MLYDEFN
#define MLYDEFN MLCBDEFN
#endif
#ifndef WSYDEFN
#define WSYDEFN WSCBDEFN
#endif
#ifndef MLMPROC
#define MLMPROC MLCBPROC
#endif
#ifndef WSMPROC
#define WSMPROC WSCBPROC
#endif
#ifndef MLMDECL
#define MLMDECL MLCBDECL
#endif
#ifndef WSMDECL
#define WSMDECL WSCBDECL
#endif
#ifndef MLMDEFN
#define MLMDEFN MLCBDEFN
#endif
#ifndef WSMDEFN
#define WSMDEFN WSCBDEFN
#endif


/* move into mlmake.h */
#ifndef MLUPROC
#define MLUPROC MLCBPROC
#endif
#ifndef MLUDECL
#define MLUDECL MLCBDECL
#endif
#ifndef MLUDEFN
#define MLUDEFN MLCBDEFN
#endif


/* move into mlmake.h */
#ifndef MLBPROC
#define MLBPROC MLCBPROC
#endif
#ifndef MLBDECL
#define MLBDECL MLCBDECL
#endif
#ifndef MLBDEFN
#define MLBDEFN MLCBDEFN
#endif

#ifndef MLDMPROC
#define MLDMPROC MLCBPROC
#endif
#ifndef MLDMDECL
#define MLDMDECL MLCBDECL
#endif
#ifndef MLDMDEFN
#define MLDMDEFN MLCBDEFN
#endif

#ifndef MLWPROC
#define MLWPROC MLCBPROC
#endif

#ifndef __uint_ct__
#define __uint_ct__ unsigned int
#endif
#ifndef __int_ct__
#define __int_ct__ int
#endif


typedef unsigned char        uchar_ct;
typedef uchar_ct       * ucharp_ct;
typedef ucharp_ct      * ucharpp_ct;
typedef ucharpp_ct     * ucharppp_ct;
typedef unsigned short       ushort_ct;
typedef ushort_ct      * ushortp_ct;
typedef ushortp_ct     * ushortpp_ct;
typedef ushortpp_ct    * ushortppp_ct;
typedef __uint_ct__          uint_ct;
typedef __uint_ct__    * uintp_ct;
typedef uintp_ct       * uintpp_ct;
typedef __int_ct__           int_ct;
typedef void           * voidp_ct;
typedef voidp_ct       * voidpp_ct;
typedef char           * charp_ct;
typedef charp_ct       * charpp_ct;
typedef charpp_ct      * charppp_ct;
typedef long                 long_ct;
typedef long_ct        * longp_ct;
typedef longp_ct       * longpp_ct;
typedef unsigned long        ulong_ct;
typedef ulong_ct       * ulongp_ct;




typedef const unsigned short * kushortp_ct;
typedef const unsigned short * * kushortpp_ct;
typedef const unsigned int * kuintp_ct;
typedef const unsigned int * * kuintpp_ct;
typedef const unsigned char * kucharp_ct;
typedef const unsigned char * * kucharpp_ct;
typedef const char * kcharp_ct;
typedef const char * * kcharpp_ct;
typedef const void * kvoidp_ct;


typedef void * MLPointer;

#ifndef __MLENVPARAM__
	typedef void * MLENVPARAM;
	typedef MLENVPARAM MLEnvironmentParameter;
#define __MLENVPARAM__
#endif

#ifndef __MLENV__
	typedef struct ml_environment *MLENV;
	typedef MLENV MLEnvironment;
#define __MLENV__ /* Don't rename...WolframLibrary.h depends on this name */
#endif

#ifndef __MLINK__
	typedef struct MLink *MLINK;
#define __MLINK__ /* Don't rename...WolframLibrary.h depends on this name */
#endif

#ifndef __MLMARK__
	typedef struct MLinkMark *MLMARK;
	typedef MLMARK MLINKMark;
#define __MLMARK__
#endif

#ifndef __mlapi_token__
#define __mlapi_token__ int_ct
#endif
typedef __mlapi_token__   mlapi_token;


typedef unsigned long      mlapi__token;
typedef mlapi__token * mlapi__tokenp;

#ifndef __mlapi_packet__
#define __mlapi_packet__ int_ct
#endif
typedef __mlapi_packet__  mlapi_packet;


typedef long mlapi_error;
typedef long mlapi__error;

typedef long  long_st;
typedef longp_ct longp_st;
typedef longp_st* longpp_st;

typedef long long_et;


#ifndef __mlapi_result__
#define __mlapi_result__ int_ct
#endif
typedef __mlapi_result__ mlapi_result;


#define MLSUCCESS (1) /*bugcheck:  this stuff doesnt belong where it can be seen at MLAPI_ layer */
#define MLFAILURE (0)

ML_EXTERN_C

#if WINDOWS_MATHLINK
typedef int (__stdcall *__MLProcPtr__)();
#else
typedef long (* __MLProcPtr__)(void);
#endif

ML_END_EXTERN_C

#endif /* _MLTYPES_H */



ML_EXTERN_C



typedef void * dev_voidp;
typedef dev_voidp dev_type;
typedef dev_type * dev_typep;
typedef long devproc_error;
typedef unsigned long devproc_selector;


#define MLDEV_WRITE_WINDOW  0
#define MLDEV_WRITE         1
#define MLDEV_HAS_DATA      2
#define MLDEV_READ          3
#define MLDEV_READ_COMPLETE 4
#define MLDEV_ACKNOWLEDGE   5

#define T_DEV_WRITE_WINDOW  MLDEV_WRITE_WINDOW
#define T_DEV_WRITE         MLDEV_WRITE
#define T_DEV_HAS_DATA      MLDEV_HAS_DATA
#define T_DEV_READ          MLDEV_READ
#define T_DEV_READ_COMPLETE MLDEV_READ_COMPLETE


typedef struct read_buf {
	unsigned short length;
	unsigned char* ptr;
} read_buf;

typedef read_buf * read_bufp;
typedef read_bufp * read_bufpp;

MLDMPROC( devproc_error, MLDeviceProcPtr, ( dev_type dev, devproc_selector selector, dev_voidp p1, dev_voidp p2));
MLDMDECL( devproc_error, MLDeviceMain, ( dev_type dev, devproc_selector selector, dev_voidp p1, dev_voidp p2));

typedef MLDeviceProcPtr MLDeviceUPP;
#define CallMLDeviceProc(userRoutine, thing, selector, p1, p2) (*(userRoutine))((thing), (selector), (dev_voidp)(p1), (dev_voidp)(p2))
#define NewMLDeviceProc(userRoutine) (userRoutine)

typedef MLDeviceUPP dev_main_type;
typedef dev_main_type * dev_main_typep;

ML_END_EXTERN_C


#endif /* _MLDEVICE_H */



#ifndef _MLAPI_H
#define _MLAPI_H









ML_EXTERN_C

#ifndef _MLALLOC_H
#define _MLALLOC_H





#if WIN64_MATHLINK
MLAPROC( void*, MLAllocatorProcPtr, (unsigned __int64));
#else
MLAPROC( void*, MLAllocatorProcPtr, (unsigned long));
#endif



typedef MLAllocatorProcPtr MLAllocatorUPP;
#define CallMLAllocatorProc(userRoutine, size) (*(userRoutine))((size))
#define NewMLAllocatorProc(userRoutine) (userRoutine)




MLFPROC( void, MLDeallocatorProcPtr, (void*));

typedef MLDeallocatorProcPtr MLDeallocatorUPP;
#define CallMLDeallocatorProc(userRoutine, p) (*(userRoutine))((p))
#define NewMLDeallocatorProc(userRoutine) (userRoutine)



#endif /* _MLALLOC_H */


/* explicitly not protected by _MLALLOC_H in case MLDECL is redefined for multiple inclusion */


/* just some type-safe casts */
MLDECL( __MLProcPtr__, MLAllocatorCast,   ( MLAllocatorProcPtr f));
MLDECL( __MLProcPtr__, MLDeallocatorCast, ( MLDeallocatorProcPtr f));

ML_END_EXTERN_C


typedef MLAllocatorUPP MLAllocator;
typedef MLAllocator * MLAllocatorp;
#define MLCallAllocator CallMLAllocatorProc
#define MLNewAllocator NewMLAllocatorProc

typedef MLDeallocatorUPP MLDeallocator;
typedef MLDeallocator * MLDeallocatorp;
#define MLCallDeallocator CallMLDeallocatorProc
#define MLNewDeallocator NewMLDeallocatorProc

#define MLallocator MLAllocator
#define MLdeallocator MLDeallocator

#endif /* _MLAPI_H */



#ifndef _MLNTYPES_H
#define _MLNTYPES_H




#ifndef _MLNUMENV_H
#define _MLNUMENV_H





/* mlne__s2 must convert empty strings to zero */



ML_EXTERN_C


#define REALBIT 4
#define REAL_MASK (1 << REALBIT)
#define XDRBIT 5
#define XDR_MASK (1 << XDRBIT)
#define BINARYBIT 7
#define BINARY_MASK (1 << BINARYBIT)
#define SIZEVARIANTBIT 6
#define SIZEVARIANT_MASK (1 << SIZEVARIANTBIT)


#define MLTK_INVALID                                          155


#define MLNE__IMPLIED_SIZE( tok, num_dispatch) ((tok) & XDR_MASK || !((tok) & SIZEVARIANT_MASK) \
		? (tok) & 0x08 ? (tok) & (0x0E + 2) : (1 << ((tok)>>1 & 0x03)) \
		: call_num_dispatch( (num_dispatch), MLNE__SIZESELECTOR((tok)), 0,0,0))

/* Range[-128, 127] */
/* 160 -> ((unsigned char)'\240') */
#define	MLTK_8BIT_SIGNED_2sCOMPLEMENT_INTEGER                 160
/* Range[0, 255] */
/* 161 -> ((unsigned char)'\241') */
#define	MLTK_8BIT_UNSIGNED_2sCOMPLEMENT_INTEGER               161
#define MLTK_8BIT_UNSIGNED_INTEGER MLTK_8BIT_UNSIGNED_2sCOMPLEMENT_INTEGER

/* Range[-32768, 32767] */
/* 162 -> ((unsigned char)'\242') */
#define	MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER      162
/* Range[0, 65535] */
/* 163 -> ((unsigned char)'\243') */
#define	MLTK_16BIT_UNSIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER    163
#define	MLTK_16BIT_UNSIGNED_BIGENDIAN_INTEGER MLTK_16BIT_UNSIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
/* Range[-2147483648, 2147483647] */
/* 164 -> ((unsigned char)'\244') */
#define	MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER      164
/* Range[0, 4294967295] */
/* 165 -> ((unsigned char)'\245') */
#define	MLTK_32BIT_UNSIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER    165
#define	MLTK_32BIT_UNSIGNED_BIGENDIAN_INTEGER MLTK_32BIT_UNSIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
/* Range[-9223372036854775808, 9223372036854775807] */
/* 166 -> ((unsigned char)'\246') */
#define	MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER      166
/* Range[0, 18446744073709551615] */
/* 167 -> ((unsigned char)'\247') */
#define	MLTK_64BIT_UNSIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER    167
#define	MLTK_64BIT_UNSIGNED_BIGENDIAN_INTEGER MLTK_64BIT_UNSIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER


/* Range[-32768, 32767] */
/* 226 -> ((unsigned char)'\342') */
#define	MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER   226
/* Range[0, 65535] */
/* 227 -> ((unsigned char)'\343') */
#define	MLTK_16BIT_UNSIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER 227
#define	MLTK_16BIT_UNSIGNED_LITTLEENDIAN_INTEGER MLTK_16BIT_UNSIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
/* Range[-2147483648, 2147483647] */
/* 228 -> ((unsigned char)'\344') */
#define	MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER   228
/* Range[0, 4294967295] */
/* 229 -> ((unsigned char)'\345') */
#define	MLTK_32BIT_UNSIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER 229
#define	MLTK_32BIT_UNSIGNED_LITTLEENDIAN_INTEGER MLTK_32BIT_UNSIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
/* Range[-9223372036854775808, 9223372036854775807] */
/* 230 -> ((unsigned char)'\346') */
#define	MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER   230
/* Range[0, 18446744073709551615] */
/* 231 -> ((unsigned char)'\347') */
#define	MLTK_64BIT_UNSIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER 231
#define	MLTK_64BIT_UNSIGNED_LITTLEENDIAN_INTEGER MLTK_64BIT_UNSIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER

/* Interval[{-3.402823e+38, 3.402823e+38}] */
/* 180 -> ((unsigned char)'\264')    10110100   */
#define	MLTK_BIGENDIAN_IEEE754_SINGLE	                      180
/* Interval[{-1.79769313486232e+308, 1.79769313486232e+308}] */
/* 182 -> ((unsigned char)'\266')    10110110   */
#define	MLTK_BIGENDIAN_IEEE754_DOUBLE	                      182

/* 184 -> ((unsigned char)'\270')    10111000   */
#define MLTK_BIGENDIAN_128BIT_DOUBLE                          184

/* Interval[{-3.402823e+38, 3.402823e+38}] */
/* 244 -> ((unsigned char)'\364')    11110100   */
#define	MLTK_LITTLEENDIAN_IEEE754_SINGLE	                  244
/* Interval[{-1.79769313486232e+308, 1.79769313486232e+308}] */
/* 246 -> ((unsigned char)'\366')    11110110   */
#define	MLTK_LITTLEENDIAN_IEEE754_DOUBLE	                  246

/* 248 -> ((unsigned char)'\370')    11111000   */
#define MLTK_LITTLEENDIAN_128BIT_DOUBLE                       248


/* Note, if the future brings...
 * #define MLTK_128BIT_UNSIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER   ((unsigned char)'\257')
 * with  Range[0, 340282366920938463463374607431768211456 (*approximately 3.40282e+38*)]
 * the dynamic range is still a monotonically increasing function of the token value.
 * An implementation might choose to set the high varient bit to mainain this property
 * and dispatch more efficiently by avoiding overflow checks
 */

#define MLNE__SELECTOR( dtok, stok) \
	(((dtok) << 8) | (stok)) /* maybe should mask of high word and cast stok */
#define MLNE__SIZESELECTOR( tok) MLNE__SELECTOR( 0, tok)
#define MLNE__INITSELECTOR (0)
#define MLNE__TOSTRINGSELECTOR( tok) MLNE__SELECTOR( MLNE__IS_REAL(tok) ? MLTKREAL : MLTKINT, tok)
#define MLNE__FROMSTRINGSELECTOR( dtok, stok) MLNE__SELECTOR( dtok, stok)

#define MLNE__STOK( selector) ( (selector) & 0x000000FF)
#define MLNE__DTOK( selector) ( ((selector) & 0x0000FF00)>>8)

#define MLNE__IS_BINARY( tok) ((tok) & BINARY_MASK)
#define MLNE__IS_REAL( tok) ((tok) & REAL_MASK)
#define MLNE__TEXT_TOKEN( tok) (MLNE__IS_REAL( tok) ? MLTKREAL : MLTKINT)




ML_END_EXTERN_C


#endif /* _MLNUMENV_H */




/****************  Special Token types: ****************/

/* MLTK_CSHORT_P         193
   MLTK_CINT_P           194
   MLTK_CLONG_P          195
   MLTK_CFLOAT_P         209
   MLTK_CDOUBLE_P        210
   MLTK_CLONGDOUBLE_P    211 */

#define MLTK_CSHORT_P       (( BINARY_MASK | SIZEVARIANT_MASK | 1))
#define MLTK_CINT_P         (( BINARY_MASK | SIZEVARIANT_MASK | 2))
#define MLTK_CLONG_P        (( BINARY_MASK | SIZEVARIANT_MASK | 3))
#define MLTK_CFLOAT_P       (( BINARY_MASK | SIZEVARIANT_MASK | REAL_MASK | 1))
#define MLTK_CDOUBLE_P      (( BINARY_MASK | SIZEVARIANT_MASK | REAL_MASK | 2))
#define MLTK_CLONGDOUBLE_P  (( BINARY_MASK | SIZEVARIANT_MASK | REAL_MASK | 3))


#define MLTK_64BIT_LITTLEENDIAN_STRUCTURE 196
#define MLTK_64BIT_BIGENDIAN_STRUCTURE    197

/* 158 -> ((unsigned char)'\236') - used in Solaris numerics definitions */
#define MLTK_128BIT_EXTENDED 158
#define MLTK_128BIT_LONGDOUBLE 158


/* Interval[{-1.189731495357231765e+4932, 1.189731495357231765e+4932}] */
/* 218 -> ((unsigned char)'\332') */
#define MLTK_96BIT_HIGHPADDED_INTEL_80BIT_EXTENDED 218

/* Interval[{-1.189731495357231765e+4932, 1.189731495357231765e+4932}] */
/* ((unsigned char)'\330') */
#define MLTK_INTEL_80BIT_EXTENDED 216

/********************  MASTIFF  ****************************/
#define MLMASTIFF_NUMERICS_ID    "mastiff"
#define MLMASTIFF_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLMASTIFF_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLMASTIFF_CLONG          MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLMASTIFF_CINT64         MLTK_64BIT_BIGENDIAN_STRUCTURE
#define MLMASTIFF_CSIZE_T        MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLMASTIFF_CFLOAT         MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLMASTIFF_CDOUBLE        MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLMASTIFF_CLONGDOUBLE    MLTK_128BIT_EXTENDED
#define MLMASTIFF_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLMASTIFF_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLMASTIFF_MLLONG         MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLMASTIFF_MLSIZE_T       MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLMASTIFF_MLINT64        MLTK_64BIT_BIGENDIAN_STRUCTURE
#define MLMASTIFF_MLFLOAT        MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLMASTIFF_MLDOUBLE       MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLMASTIFF_MLLONGDOUBLE   MLTK_128BIT_EXTENDED

/********************  JAPANESECHIN  ****************************/
#define MLJAPANESECHIN_NUMERICS_ID    "japanesechin"
#define MLJAPANESECHIN_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLJAPANESECHIN_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLJAPANESECHIN_CLONG          MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLJAPANESECHIN_CINT64         MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLJAPANESECHIN_CSIZE_T        MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLJAPANESECHIN_CFLOAT         MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLJAPANESECHIN_CDOUBLE        MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLJAPANESECHIN_CLONGDOUBLE    MLTK_128BIT_EXTENDED
#define MLJAPANESECHIN_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLJAPANESECHIN_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLJAPANESECHIN_MLLONG         MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLJAPANESECHIN_MLINT64        MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLJAPANESECHIN_MLSIZE_T       MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLJAPANESECHIN_MLFLOAT        MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLJAPANESECHIN_MLDOUBLE       MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLJAPANESECHIN_MLLONGDOUBLE   MLTK_128BIT_EXTENDED

/********************  BORZOI  ****************************/
/* The borzoi numerics environment specifically does not have MLBORZOI_CLONGDOUBLE or
MLBORZOI_MLLONGDOUBLE */

#define MLBORZOI_NUMERICS_ID          "borzoi"
#define MLBORZOI_CSHORT               MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBORZOI_CINT                 MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBORZOI_CLONG                MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBORZOI_CSIZE_T              MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBORZOI_CINT64               MLTK_64BIT_BIGENDIAN_STRUCTURE
#define MLBORZOI_CFLOAT               MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLBORZOI_CDOUBLE              MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLBORZOI_MLSHORT              MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBORZOI_MLINT                MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBORZOI_MLLONG               MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBORZOI_MLSIZE_T             MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBORZOI_MLINT64              MLTK_64BIT_BIGENDIAN_STRUCTURE
#define MLBORZOI_MLFLOAT              MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLBORZOI_MLDOUBLE             MLTK_BIGENDIAN_IEEE754_DOUBLE

/********************  BRIARD  ****************************/
/* The briard numerics environment purposefully does not have MLBRIARD_CLONGDOUBLE or
MLBRIARD_MLLONGDOUBLE */

#define MLBRIARD_NUMERICS_ID          "briard"
#define MLBRIARD_CSHORT               MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBRIARD_CINT                 MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBRIARD_CLONG                MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBRIARD_CINT64               MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBRIARD_CSIZE_T              MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBRIARD_CFLOAT               MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLBRIARD_CDOUBLE              MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLBRIARD_MLSHORT              MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBRIARD_MLINT                MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBRIARD_MLLONG               MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBRIARD_MLINT64              MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBRIARD_MLSIZE_T             MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBRIARD_MLFLOAT              MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLBRIARD_MLDOUBLE             MLTK_BIGENDIAN_IEEE754_DOUBLE

/********************  KEESHOND  ****************************/
#define MLKEESHOND_NUMERICS_ID    "keeshond"
#define MLKEESHOND_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLKEESHOND_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLKEESHOND_CLONG          MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLKEESHOND_CSIZE_T        MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLKEESHOND_CINT64         MLTK_64BIT_BIGENDIAN_STRUCTURE
#define MLKEESHOND_CFLOAT         MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLKEESHOND_CDOUBLE        MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLKEESHOND_CLONGDOUBLE    MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLKEESHOND_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLKEESHOND_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLKEESHOND_MLLONG         MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLKEESHOND_MLSIZE_T       MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLKEESHOND_MLINT64        MLTK_64BIT_BIGENDIAN_STRUCTURE
#define MLKEESHOND_MLFLOAT        MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLKEESHOND_MLDOUBLE       MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLKEESHOND_MLLONGDOUBLE   MLTK_BIGENDIAN_IEEE754_DOUBLE

/********************  KOMONDOR  ****************************/
#define MLKOMONDOR_NUMERICS_ID    "komondor"
#define MLKOMONDOR_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLKOMONDOR_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLKOMONDOR_CLONG          MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLKOMONDOR_CINT64         MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLKOMONDOR_CSIZE_T        MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLKOMONDOR_CFLOAT         MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLKOMONDOR_CDOUBLE        MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLKOMONDOR_CLONGDOUBLE    MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLKOMONDOR_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLKOMONDOR_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLKOMONDOR_MLLONG         MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLKOMONDOR_MLSIZE_T       MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLKOMONDOR_MLINT64        MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLKOMONDOR_MLFLOAT        MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLKOMONDOR_MLDOUBLE       MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLKOMONDOR_MLLONGDOUBLE   MLTK_BIGENDIAN_IEEE754_DOUBLE

/********************  NORWEGIANELKHOUND  ****************************/
#define MLNORWEGIANELKHOUND_NUMERICS_ID    "norwegianelkhound"
#define MLNORWEGIANELKHOUND_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLNORWEGIANELKHOUND_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLNORWEGIANELKHOUND_CLONG          MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLNORWEGIANELKHOUND_CINT64         MLTK_64BIT_LITTLEENDIAN_STRUCTURE
#define MLNORWEGIANELKHOUND_CSIZE_T        MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLNORWEGIANELKHOUND_CFLOAT         MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLNORWEGIANELKHOUND_CDOUBLE        MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLNORWEGIANELKHOUND_CLONGDOUBLE    MLTK_96BIT_HIGHPADDED_INTEL_80BIT_EXTENDED
#define MLNORWEGIANELKHOUND_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLNORWEGIANELKHOUND_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLNORWEGIANELKHOUND_MLLONG         MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLNORWEGIANELKHOUND_MLSIZE_T       MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLNORWEGIANELKHOUND_MLINT64        MLTK_64BIT_LITTLEENDIAN_STRUCTURE
#define MLNORWEGIANELKHOUND_MLFLOAT        MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLNORWEGIANELKHOUND_MLDOUBLE       MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLNORWEGIANELKHOUND_MLLONGDOUBLE   MLTK_96BIT_HIGHPADDED_INTEL_80BIT_EXTENDED

/********************  NORWICHTERRIOR  ****************************/
#define MLNORWICHTERRIOR_NUMERICS_ID    "norwichterrior"
#define MLNORWICHTERRIOR_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLNORWICHTERRIOR_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLNORWICHTERRIOR_CLONG          MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLNORWICHTERRIOR_CINT64         MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLNORWICHTERRIOR_CSIZE_T        MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLNORWICHTERRIOR_CFLOAT         MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLNORWICHTERRIOR_CDOUBLE        MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLNORWICHTERRIOR_CLONGDOUBLE    MLTK_LITTLEENDIAN_128BIT_DOUBLE
#define MLNORWICHTERRIOR_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLNORWICHTERRIOR_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLNORWICHTERRIOR_MLLONG         MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLNORWICHTERRIOR_MLSIZE_T       MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLNORWICHTERRIOR_MLINT64        MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLNORWICHTERRIOR_MLFLOAT        MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLNORWICHTERRIOR_MLDOUBLE       MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLNORWICHTERRIOR_MLLONGDOUBLE   MLTK_LITTLEENDIAN_128BIT_DOUBLE

/********************  SAINTBERNARD  ****************************/
#define MLSAINTBERNARD_NUMERICS_ID    "saintbernarnd"
#define MLSAINTBERNARD_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLSAINTBERNARD_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLSAINTBERNARD_CLONG          MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLSAINTBERNARD_CINT64         MLTK_64BIT_BIGENDIAN_STRUCTURE
#define MLSAINTBERNARD_CSIZE_T        MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLSAINTBERNARD_CFLOAT         MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLSAINTBERNARD_CDOUBLE        MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLSAINTBERNARD_CLONGDOUBLE    MLTK_BIGENDIAN_128BIT_DOUBLE
#define MLSAINTBERNARD_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLSAINTBERNARD_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLSAINTBERNARD_MLLONG         MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLSAINTBERNARD_MLINT64        MLTK_64BIT_BIGENDIAN_STRUCTURE
#define MLSAINTBERNARD_MLSIZE_T       MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLSAINTBERNARD_MLFLOAT        MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLSAINTBERNARD_MLDOUBLE       MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLSAINTBERNARD_MLLONGDOUBLE   MLTK_BIGENDIAN_128BIT_DOUBLE

/********************  BERNESEMOUNTAINDOG  ****************************/
#define MLBERNESEMOUNTAINDOG_NUMERICS_ID    "bernesemountaindog"
#define MLBERNESEMOUNTAINDOG_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBERNESEMOUNTAINDOG_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBERNESEMOUNTAINDOG_CLONG          MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBERNESEMOUNTAINDOG_CINT64         MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBERNESEMOUNTAINDOG_CSIZE_T        MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBERNESEMOUNTAINDOG_CFLOAT         MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLBERNESEMOUNTAINDOG_CDOUBLE        MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLBERNESEMOUNTAINDOG_CLONGDOUBLE    MLTK_BIGENDIAN_128BIT_DOUBLE
#define MLBERNESEMOUNTAINDOG_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBERNESEMOUNTAINDOG_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBERNESEMOUNTAINDOG_MLLONG         MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBERNESEMOUNTAINDOG_MLINT64        MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBERNESEMOUNTAINDOG_MLSIZE_T       MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBERNESEMOUNTAINDOG_MLFLOAT        MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLBERNESEMOUNTAINDOG_MLDOUBLE       MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLBERNESEMOUNTAINDOG_MLLONGDOUBLE   MLTK_BIGENDIAN_128BIT_DOUBLE

/********************  SETTER  ****************************/
#define MLSETTER_NUMERICS_ID    "setter"
#define MLSETTER_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLSETTER_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLSETTER_CLONG          MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLSETTER_CINT64         MLTK_64BIT_LITTLEENDIAN_STRUCTURE
#define MLSETTER_CSIZE_T        MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLSETTER_CFLOAT         MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLSETTER_CDOUBLE        MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLSETTER_CLONGDOUBLE    MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLSETTER_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLSETTER_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLSETTER_MLLONG         MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLSETTER_MLINT64        MLTK_64BIT_LITTLEENDIAN_STRUCTURE
#define MLSETTER_MLSIZE_T       MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLSETTER_MLFLOAT        MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLSETTER_MLDOUBLE       MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLSETTER_MLLONGDOUBLE   MLTK_INTEL_80BIT_EXTENDED

/********************  FRENCH_BULLDOG  ****************************/
#define MLFRENCH_BULLDOG_NUMERICS_ID    "french_bulldog"
#define MLFRENCH_BULLDOG_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLFRENCH_BULLDOG_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLFRENCH_BULLDOG_CLONG          MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLFRENCH_BULLDOG_CINT64         MLTK_64BIT_LITTLEENDIAN_STRUCTURE
#define MLFRENCH_BULLDOG_CSIZE_T        MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLFRENCH_BULLDOG_CFLOAT         MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLFRENCH_BULLDOG_CDOUBLE        MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLFRENCH_BULLDOG_CLONGDOUBLE    MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLFRENCH_BULLDOG_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLFRENCH_BULLDOG_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLFRENCH_BULLDOG_MLLONG         MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLFRENCH_BULLDOG_MLINT64        MLTK_64BIT_LITTLEENDIAN_STRUCTURE
#define MLFRENCH_BULLDOG_MLSIZE_T       MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLFRENCH_BULLDOG_MLFLOAT        MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLFRENCH_BULLDOG_MLDOUBLE       MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLFRENCH_BULLDOG_MLLONGDOUBLE   MLTK_LITTLEENDIAN_IEEE754_DOUBLE

/********************  BICHON_FRISE  ****************************/
#define MLBICHON_FRISE_NUMERICS_ID    "bichon_frise"
#define MLBICHON_FRISE_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBICHON_FRISE_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBICHON_FRISE_CLONG          MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBICHON_FRISE_CINT64         MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBICHON_FRISE_CSIZE_T        MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBICHON_FRISE_CFLOAT         MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLBICHON_FRISE_CDOUBLE        MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLBICHON_FRISE_CLONGDOUBLE    MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLBICHON_FRISE_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBICHON_FRISE_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBICHON_FRISE_MLLONG         MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBICHON_FRISE_MLINT64        MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBICHON_FRISE_MLSIZE_T       MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBICHON_FRISE_MLFLOAT        MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLBICHON_FRISE_MLDOUBLE       MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLBICHON_FRISE_MLLONGDOUBLE   MLTK_LITTLEENDIAN_IEEE754_DOUBLE

/********************  HELEN  ****************************/
#define MLHELEN_NUMERICS_ID    "helen"
#define MLHELEN_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLHELEN_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLHELEN_CLONG          MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLHELEN_CINT64         MLTK_64BIT_LITTLEENDIAN_STRUCTURE
#define MLHELEN_CSIZE_T        MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLHELEN_CFLOAT         MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLHELEN_CDOUBLE        MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLHELEN_CLONGDOUBLE    MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLHELEN_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLHELEN_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLHELEN_MLLONG         MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLHELEN_MLINT64        MLTK_64BIT_LITTLEENDIAN_STRUCTURE
#define MLHELEN_MLSIZE_T       MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLHELEN_MLFLOAT        MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLHELEN_MLDOUBLE       MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLHELEN_MLLONGDOUBLE   MLTK_LITTLEENDIAN_IEEE754_DOUBLE

/********************  BEAGLE  ****************************/
#define MLBEAGLE_NUMERICS_ID    "beagle"
#define MLBEAGLE_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBEAGLE_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBEAGLE_CLONG          MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBEAGLE_CINT64         MLTK_64BIT_BIGENDIAN_STRUCTURE
#define MLBEAGLE_CSIZE_T        MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBEAGLE_CFLOAT         MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLBEAGLE_CDOUBLE        MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLBEAGLE_CLONGDOUBLE    MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLBEAGLE_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBEAGLE_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBEAGLE_MLLONG         MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBEAGLE_MLINT64        MLTK_64BIT_BIGENDIAN_STRUCTURE
#define MLBEAGLE_MLSIZE_T       MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBEAGLE_MLFLOAT        MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLBEAGLE_MLDOUBLE       MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLBEAGLE_MLLONGDOUBLE   MLTK_BIGENDIAN_IEEE754_DOUBLE

/********************  BULLTERRIER  ****************************/
#define MLBULLTERRIER_NUMERICS_ID    "bullterrier"
#define MLBULLTERRIER_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBULLTERRIER_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBULLTERRIER_CLONG          MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBULLTERRIER_CINT64         MLTK_64BIT_BIGENDIAN_STRUCTURE
#define MLBULLTERRIER_CSIZE_T        MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBULLTERRIER_CFLOAT         MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLBULLTERRIER_CDOUBLE        MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLBULLTERRIER_CLONGDOUBLE    MLTK_BIGENDIAN_128BIT_DOUBLE
#define MLBULLTERRIER_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBULLTERRIER_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBULLTERRIER_MLLONG         MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBULLTERRIER_MLINT64        MLTK_64BIT_BIGENDIAN_STRUCTURE
#define MLBULLTERRIER_MLSIZE_T       MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBULLTERRIER_MLFLOAT        MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLBULLTERRIER_MLDOUBLE       MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLBULLTERRIER_MLLONGDOUBLE   MLTK_BIGENDIAN_128BIT_DOUBLE

/********************  BORDERTERRIER  ****************************/
#define MLBORDERTERRIER_NUMERICS_ID    "borderterrier"
#define MLBORDERTERRIER_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBORDERTERRIER_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBORDERTERRIER_CLONG          MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBORDERTERRIER_CINT64         MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBORDERTERRIER_CSIZE_T        MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBORDERTERRIER_CFLOAT         MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLBORDERTERRIER_CDOUBLE        MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLBORDERTERRIER_CLONGDOUBLE    MLTK_BIGENDIAN_128BIT_DOUBLE
#define MLBORDERTERRIER_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBORDERTERRIER_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBORDERTERRIER_MLLONG         MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBORDERTERRIER_MLINT64        MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBORDERTERRIER_MLSIZE_T       MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLBORDERTERRIER_MLFLOAT        MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLBORDERTERRIER_MLDOUBLE       MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLBORDERTERRIER_MLLONGDOUBLE   MLTK_BIGENDIAN_128BIT_DOUBLE

/********************  BASENJI  ****************************/
#define MLBASENJI_NUMERICS_ID    "basenji"
#define MLBASENJI_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBASENJI_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBASENJI_CLONG          MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBASENJI_CINT64         MLTK_64BIT_LITTLEENDIAN_STRUCTURE
#define MLBASENJI_CSIZE_T        MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBASENJI_CFLOAT         MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLBASENJI_CDOUBLE        MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLBASENJI_CLONGDOUBLE    MLTK_LITTLEENDIAN_128BIT_DOUBLE
#define MLBASENJI_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBASENJI_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBASENJI_MLLONG         MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBASENJI_MLINT64        MLTK_64BIT_LITTLEENDIAN_STRUCTURE
#define MLBASENJI_MLSIZE_T       MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBASENJI_MLFLOAT        MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLBASENJI_MLDOUBLE       MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLBASENJI_MLLONGDOUBLE   MLTK_LITTLEENDIAN_128BIT_DOUBLE

/********************  SHARPEI  ****************************/
#define MLSHARPEI_NUMERICS_ID    "sharpei"
#define MLSHARPEI_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLSHARPEI_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLSHARPEI_CLONG          MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLSHARPEI_CINT64         MLTK_64BIT_LITTLEENDIAN_STRUCTURE
#define MLSHARPEI_CSIZE_T        MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLSHARPEI_CFLOAT         MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLSHARPEI_CDOUBLE        MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLSHARPEI_CLONGDOUBLE    MLTK_LITTLEENDIAN_128BIT_DOUBLE
#define MLSHARPEI_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLSHARPEI_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLSHARPEI_MLLONG         MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLSHARPEI_MLINT64        MLTK_64BIT_LITTLEENDIAN_STRUCTURE
#define MLSHARPEI_MLSIZE_T       MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLSHARPEI_MLFLOAT        MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLSHARPEI_MLDOUBLE       MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLSHARPEI_MLLONGDOUBLE   MLTK_LITTLEENDIAN_128BIT_DOUBLE

/********************  TIBETANMASTIFF  ****************************/
#define MLTIBETANMASTIFF_NUMERICS_ID    "tibetanmastiff"
#define MLTIBETANMASTIFF_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLTIBETANMASTIFF_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLTIBETANMASTIFF_CLONG          MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLTIBETANMASTIFF_CINT64         MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLTIBETANMASTIFF_CSIZE_T        MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLTIBETANMASTIFF_CFLOAT         MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLTIBETANMASTIFF_CDOUBLE        MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLTIBETANMASTIFF_CLONGDOUBLE    MLTK_LITTLEENDIAN_128BIT_DOUBLE
#define MLTIBETANMASTIFF_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLTIBETANMASTIFF_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLTIBETANMASTIFF_MLLONG         MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLTIBETANMASTIFF_MLINT64        MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLTIBETANMASTIFF_MLSIZE_T       MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLTIBETANMASTIFF_MLFLOAT        MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLTIBETANMASTIFF_MLDOUBLE       MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLTIBETANMASTIFF_MLLONGDOUBLE   MLTK_LITTLEENDIAN_128BIT_DOUBLE

/********************  GREATDANE  ****************************/
#define MLGREATDANE_NUMERICS_ID    "greatdane"
#define MLGREATDANE_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLGREATDANE_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLGREATDANE_CLONG          MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLGREATDANE_CINT64         MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLGREATDANE_CSIZE_T        MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLGREATDANE_CFLOAT         MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLGREATDANE_CDOUBLE        MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLGREATDANE_CLONGDOUBLE    MLTK_LITTLEENDIAN_128BIT_DOUBLE
#define MLGREATDANE_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLGREATDANE_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLGREATDANE_MLLONG         MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLGREATDANE_MLINT64        MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLGREATDANE_MLSIZE_T       MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLGREATDANE_MLFLOAT        MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLGREATDANE_MLDOUBLE       MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLGREATDANE_MLLONGDOUBLE   MLTK_LITTLEENDIAN_128BIT_DOUBLE

/********************  REDDOG  ****************************/
#define MLREDDOG_NUMERICS_ID    "reddog"
#define MLREDDOG_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLREDDOG_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLREDDOG_CLONG          MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLREDDOG_CINT64         MLTK_64BIT_LITTLEENDIAN_STRUCTURE
#define MLREDDOG_CSIZE_T        MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLREDDOG_CFLOAT         MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLREDDOG_CDOUBLE        MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLREDDOG_CLONGDOUBLE    MLTK_96BIT_HIGHPADDED_INTEL_80BIT_EXTENDED
#define MLREDDOG_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLREDDOG_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLREDDOG_MLLONG         MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLREDDOG_MLINT64        MLTK_64BIT_LITTLEENDIAN_STRUCTURE
#define MLREDDOG_MLSIZE_T       MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLREDDOG_MLFLOAT        MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLREDDOG_MLDOUBLE       MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLREDDOG_MLLONGDOUBLE   MLTK_96BIT_HIGHPADDED_INTEL_80BIT_EXTENDED

/********************  AUSTRALIANCATTLEDOG  ****************************/
#define MLAUSTRALIANCATTLEDOG_NUMERICS_ID    "australiancattledog"
#define MLAUSTRALIANCATTLEDOG_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLAUSTRALIANCATTLEDOG_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLAUSTRALIANCATTLEDOG_CLONG          MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLAUSTRALIANCATTLEDOG_CINT64         MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLAUSTRALIANCATTLEDOG_CSIZE_T        MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLAUSTRALIANCATTLEDOG_CFLOAT         MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLAUSTRALIANCATTLEDOG_CDOUBLE        MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLAUSTRALIANCATTLEDOG_CLONGDOUBLE    MLTK_LITTLEENDIAN_128BIT_DOUBLE
#define MLAUSTRALIANCATTLEDOG_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLAUSTRALIANCATTLEDOG_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLAUSTRALIANCATTLEDOG_MLLONG         MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLAUSTRALIANCATTLEDOG_MLINT64        MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLAUSTRALIANCATTLEDOG_MLSIZE_T       MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLAUSTRALIANCATTLEDOG_MLFLOAT        MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLAUSTRALIANCATTLEDOG_MLDOUBLE       MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLAUSTRALIANCATTLEDOG_MLLONGDOUBLE   MLTK_LITTLEENDIAN_128BIT_DOUBLE

/********************  BOXER  ****************************/
#define MLBOXER_NUMERICS_ID    "boxer"
#define MLBOXER_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBOXER_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBOXER_CLONG          MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBOXER_CINT64         MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBOXER_CSIZE_T        MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBOXER_CFLOAT         MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLBOXER_CDOUBLE        MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLBOXER_CLONGDOUBLE    MLTK_LITTLEENDIAN_128BIT_DOUBLE
#define MLBOXER_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBOXER_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBOXER_MLLONG         MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBOXER_MLINT64        MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBOXER_MLSIZE_T       MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBOXER_MLFLOAT        MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLBOXER_MLDOUBLE       MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLBOXER_MLLONGDOUBLE   MLTK_LITTLEENDIAN_128BIT_DOUBLE

/********************  AKITAINU  ****************************/
#define MLAKITAINU_NUMERICS_ID    "akitainu"
#define MLAKITAINU_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLAKITAINU_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLAKITAINU_CLONG          MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLAKITAINU_CINT64         MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLAKITAINU_CSIZE_T        MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLAKITAINU_CFLOAT         MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLAKITAINU_CDOUBLE        MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLAKITAINU_CLONGDOUBLE    MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLAKITAINU_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLAKITAINU_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLAKITAINU_MLLONG         MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLAKITAINU_MLINT64        MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLAKITAINU_MLSIZE_T       MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLAKITAINU_MLFLOAT        MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLAKITAINU_MLDOUBLE       MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLAKITAINU_MLLONGDOUBLE   MLTK_LITTLEENDIAN_IEEE754_DOUBLE

/********************  CHIHUAHUA  ****************************/
#define MLCHIHUAHUA_NUMERICS_ID    "chihuahua"
#define MLCHIHUAHUA_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLCHIHUAHUA_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLCHIHUAHUA_CLONG          MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLCHIHUAHUA_CINT64         MLTK_64BIT_LITTLEENDIAN_STRUCTURE
#define MLCHIHUAHUA_CSIZE_T        MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLCHIHUAHUA_CFLOAT         MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLCHIHUAHUA_CDOUBLE        MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLCHIHUAHUA_CLONGDOUBLE    MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLCHIHUAHUA_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLCHIHUAHUA_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLCHIHUAHUA_MLLONG         MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLCHIHUAHUA_MLINT64        MLTK_64BIT_LITTLEENDIAN_STRUCTURE
#define MLCHIHUAHUA_MLSIZE_T       MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLCHIHUAHUA_MLFLOAT        MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLCHIHUAHUA_MLDOUBLE       MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLCHIHUAHUA_MLLONGDOUBLE   MLTK_LITTLEENDIAN_IEEE754_DOUBLE

/********************  ROTTWEILER  ****************************/
#define MLROTTWEILER_NUMERICS_ID    "rottweiler"
#define MLROTTWEILER_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLROTTWEILER_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLROTTWEILER_CLONG          MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLROTTWEILER_CINT64         MLTK_64BIT_BIGENDIAN_STRUCTURE
#define MLROTTWEILER_CSIZE_T        MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLROTTWEILER_CFLOAT         MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLROTTWEILER_CDOUBLE        MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLROTTWEILER_CLONGDOUBLE    MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLROTTWEILER_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLROTTWEILER_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLROTTWEILER_MLLONG         MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLROTTWEILER_MLINT64        MLTK_64BIT_BIGENDIAN_STRUCTURE
#define MLROTTWEILER_MLSIZE_T       MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLROTTWEILER_MLFLOAT        MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLROTTWEILER_MLDOUBLE       MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLROTTWEILER_MLLONGDOUBLE   MLTK_BIGENDIAN_IEEE754_DOUBLE

/********************  PHARAOHHOUND  ****************************/
#define MLPHARAOHHOUND_NUMERICS_ID    "pharaohhound"
#define MLPHARAOHHOUND_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLPHARAOHHOUND_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLPHARAOHHOUND_CLONG          MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLPHARAOHHOUND_CINT64         MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLPHARAOHHOUND_CSIZE_T        MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLPHARAOHHOUND_CFLOAT         MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLPHARAOHHOUND_CDOUBLE        MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLPHARAOHHOUND_CLONGDOUBLE    MLTK_BIGENDIAN_128BIT_DOUBLE
#define MLPHARAOHHOUND_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLPHARAOHHOUND_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLPHARAOHHOUND_MLLONG         MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLPHARAOHHOUND_MLINT64        MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLPHARAOHHOUND_MLSIZE_T       MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLPHARAOHHOUND_MLFLOAT        MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLPHARAOHHOUND_MLDOUBLE       MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLPHARAOHHOUND_MLLONGDOUBLE   MLTK_BIGENDIAN_128BIT_DOUBLE

/********************  TROUT  ****************************/
#define MLTROUT_NUMERICS_ID    "trout"
#define MLTROUT_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLTROUT_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLTROUT_CLONG          MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLTROUT_CINT64         MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLTROUT_CSIZE_T        MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLTROUT_CFLOAT         MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLTROUT_CDOUBLE        MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLTROUT_CLONGDOUBLE    MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLTROUT_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLTROUT_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLTROUT_MLLONG         MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLTROUT_MLINT64        MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLTROUT_MLSIZE_T       MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLTROUT_MLFLOAT        MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLTROUT_MLDOUBLE       MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLTROUT_MLLONGDOUBLE   MLTK_BIGENDIAN_IEEE754_DOUBLE

/********************  PUG  ****************************/
#define MLPUG_NUMERICS_ID    "pug"
#define MLPUG_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLPUG_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLPUG_CLONG          MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLPUG_CINT64         MLTK_64BIT_BIGENDIAN_STRUCTURE
#define MLPUG_CSIZE_T        MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLPUG_CFLOAT         MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLPUG_CDOUBLE        MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLPUG_CLONGDOUBLE    MLTK_BIGENDIAN_128BIT_DOUBLE
#define MLPUG_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLPUG_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLPUG_MLLONG         MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLPUG_MLINT64        MLTK_64BIT_BIGENDIAN_STRUCTURE
#define MLPUG_MLSIZE_T       MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLPUG_MLFLOAT        MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLPUG_MLDOUBLE       MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLPUG_MLLONGDOUBLE   MLTK_BIGENDIAN_128BIT_DOUBLE

/********************  POINTER  ****************************/
#define MLPOINTER_NUMERICS_ID    "pointer"
#define MLPOINTER_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLPOINTER_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLPOINTER_CLONG          MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLPOINTER_CINT64         MLTK_64BIT_BIGENDIAN_STRUCTURE
#define MLPOINTER_CSIZE_T        MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLPOINTER_CFLOAT         MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLPOINTER_CDOUBLE        MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLPOINTER_CLONGDOUBLE    MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLPOINTER_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLPOINTER_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLPOINTER_MLLONG         MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLPOINTER_MLINT64        MLTK_64BIT_BIGENDIAN_STRUCTURE
#define MLPOINTER_MLSIZE_T       MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLPOINTER_MLFLOAT        MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLPOINTER_MLDOUBLE       MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLPOINTER_MLLONGDOUBLE   MLTK_BIGENDIAN_IEEE754_DOUBLE

/********************  SAMOYED  ****************************/
#define MLSAMOYED_NUMERICS_ID    "samoyed"
#define MLSAMOYED_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLSAMOYED_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLSAMOYED_CLONG          MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLSAMOYED_CINT64         MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLSAMOYED_CSIZE_T        MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLSAMOYED_CFLOAT         MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLSAMOYED_CDOUBLE        MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLSAMOYED_CLONGDOUBLE    MLTK_BIGENDIAN_128BIT_DOUBLE
#define MLSAMOYED_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLSAMOYED_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLSAMOYED_MLLONG         MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLSAMOYED_MLINT64        MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLSAMOYED_MLSIZE_T       MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLSAMOYED_MLFLOAT        MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLSAMOYED_MLDOUBLE       MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLSAMOYED_MLLONGDOUBLE   MLTK_BIGENDIAN_128BIT_DOUBLE

/********************  SIBERIANHUSKY  ****************************/
#define MLSIBERIANHUSKY_NUMERICS_ID    "siberianhusky"
#define MLSIBERIANHUSKY_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLSIBERIANHUSKY_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLSIBERIANHUSKY_CLONG          MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLSIBERIANHUSKY_CINT64         MLTK_64BIT_BIGENDIAN_STRUCTURE
#define MLSIBERIANHUSKY_CSIZE_T        MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLSIBERIANHUSKY_CFLOAT         MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLSIBERIANHUSKY_CDOUBLE        MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLSIBERIANHUSKY_CLONGDOUBLE    MLTK_BIGENDIAN_128BIT_DOUBLE
#define MLSIBERIANHUSKY_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLSIBERIANHUSKY_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLSIBERIANHUSKY_MLLONG         MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLSIBERIANHUSKY_MLINT64        MLTK_64BIT_BIGENDIAN_STRUCTURE
#define MLSIBERIANHUSKY_MLSIZE_T       MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLSIBERIANHUSKY_MLFLOAT        MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLSIBERIANHUSKY_MLDOUBLE       MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLSIBERIANHUSKY_MLLONGDOUBLE   MLTK_BIGENDIAN_128BIT_DOUBLE

/********************  SHIBAINU  ****************************/
#define MLSHIBAINU_NUMERICS_ID    "shibainu"
#define MLSHIBAINU_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLSHIBAINU_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLSHIBAINU_CLONG          MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLSHIBAINU_CINT64         MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLSHIBAINU_CSIZE_T        MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLSHIBAINU_CFLOAT         MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLSHIBAINU_CDOUBLE        MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLSHIBAINU_CLONGDOUBLE    MLTK_LITTLEENDIAN_128BIT_DOUBLE
#define MLSHIBAINU_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLSHIBAINU_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLSHIBAINU_MLLONG         MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLSHIBAINU_MLINT64        MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLSHIBAINU_MLSIZE_T       MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLSHIBAINU_MLFLOAT        MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLSHIBAINU_MLDOUBLE       MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLSHIBAINU_MLLONGDOUBLE   MLTK_LITTLEENDIAN_128BIT_DOUBLE

/********************  NEWFOUNDLAND  ****************************/
#define MLNEWFOUNDLAND_NUMERICS_ID    "newfoundland"
#define MLNEWFOUNDLAND_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLNEWFOUNDLAND_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLNEWFOUNDLAND_CLONG          MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLNEWFOUNDLAND_CINT64         MLTK_64BIT_BIGENDIAN_STRUCTURE
#define MLNEWFOUNDLAND_CSIZE_T        MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLNEWFOUNDLAND_CFLOAT         MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLNEWFOUNDLAND_CDOUBLE        MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLNEWFOUNDLAND_CLONGDOUBLE    MLTK_BIGENDIAN_128BIT_DOUBLE
#define MLNEWFOUNDLAND_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLNEWFOUNDLAND_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLNEWFOUNDLAND_MLLONG         MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLNEWFOUNDLAND_MLINT64        MLTK_64BIT_BIGENDIAN_STRUCTURE
#define MLNEWFOUNDLAND_MLSIZE_T       MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLNEWFOUNDLAND_MLFLOAT        MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLNEWFOUNDLAND_MLDOUBLE       MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLNEWFOUNDLAND_MLLONGDOUBLE   MLTK_BIGENDIAN_128BIT_DOUBLE

/********************  AFFENPINSCHER  ****************************/
#define MLAFFENPINSCHER_NUMERICS_ID    "affenpinscher"
#define MLAFFENPINSCHER_CSHORT         MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLAFFENPINSCHER_CINT           MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLAFFENPINSCHER_CLONG          MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLAFFENPINSCHER_CINT64         MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLAFFENPINSCHER_CSIZE_T        MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLAFFENPINSCHER_CFLOAT         MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLAFFENPINSCHER_CDOUBLE        MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLAFFENPINSCHER_CLONGDOUBLE    MLTK_BIGENDIAN_128BIT_DOUBLE
#define MLAFFENPINSCHER_MLSHORT        MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLAFFENPINSCHER_MLINT          MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLAFFENPINSCHER_MLLONG         MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLAFFENPINSCHER_MLINT64        MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLAFFENPINSCHER_MLSIZE_T       MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#define MLAFFENPINSCHER_MLFLOAT        MLTK_BIGENDIAN_IEEE754_SINGLE
#define MLAFFENPINSCHER_MLDOUBLE       MLTK_BIGENDIAN_IEEE754_DOUBLE
#define MLAFFENPINSCHER_MLLONGDOUBLE   MLTK_BIGENDIAN_128BIT_DOUBLE

/********************  BEAUCERON  ****************************/
#define MLBEAUCERON_NUMERICS_ID        "beauceron"
#define MLBEAUCERON_CSHORT             MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBEAUCERON_CINT               MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBEAUCERON_CLONG              MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBEAUCERON_CINT64             MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBEAUCERON_CSIZE_T            MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBEAUCERON_CFLOAT             MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLBEAUCERON_CDOUBLE            MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLBEAUCERON_CLONGDOUBLE        MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLBEAUCERON_MLSHORT            MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBEAUCERON_MLINT              MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBEAUCERON_MLLONG             MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBEAUCERON_MLINT64            MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBEAUCERON_MLSIZE_T           MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBEAUCERON_MLFLOAT            MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLBEAUCERON_MLDOUBLE           MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLBEAUCERON_MLLONGDOUBLE       MLTK_LITTLEENDIAN_IEEE754_DOUBLE


/********************  BERGAMASCO  ****************************/
#define MLBERGAMASCO_NUMERICS_ID       "bergamasco"
#define MLBERGAMASCO_CSHORT             MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBERGAMASCO_CINT               MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBERGAMASCO_CLONG              MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBERGAMASCO_CINT64             MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBERGAMASCO_CSIZE_T            MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBERGAMASCO_CFLOAT             MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLBERGAMASCO_CDOUBLE            MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLBERGAMASCO_CLONGDOUBLE        MLTK_96BIT_HIGHPADDED_INTEL_80BIT_EXTENDED
#define MLBERGAMASCO_MLSHORT            MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBERGAMASCO_MLINT              MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBERGAMASCO_MLLONG             MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBERGAMASCO_MLINT64            MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBERGAMASCO_MLSIZE_T           MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBERGAMASCO_MLFLOAT            MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLBERGAMASCO_MLDOUBLE           MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLBERGAMASCO_MLLONGDOUBLE       MLTK_96BIT_HIGHPADDED_INTEL_80BIT_EXTENDED


/********************  BOERBOEL  ****************************/
#define MLBOERBOEL_NUMERICS_ID       "boerboel"
#define MLBOERBOEL_CSHORT             MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBOERBOEL_CINT               MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBOERBOEL_CLONG              MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBOERBOEL_CINT64             MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBOERBOEL_CSIZE_T            MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBOERBOEL_CFLOAT             MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLBOERBOEL_CDOUBLE            MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLBOERBOEL_CLONGDOUBLE        MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLBOERBOEL_MLSHORT            MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBOERBOEL_MLINT              MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBOERBOEL_MLLONG             MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBOERBOEL_MLINT64            MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBOERBOEL_MLSIZE_T           MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLBOERBOEL_MLFLOAT            MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLBOERBOEL_MLDOUBLE           MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLBOERBOEL_MLLONGDOUBLE       MLTK_LITTLEENDIAN_IEEE754_DOUBLE


/********************  CHINOOK  ****************************/
#define MLCHINOOK_NUMERICS_ID       "chinook"
#define MLCHINOOK_CSHORT             MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLCHINOOK_CINT               MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLCHINOOK_CLONG              MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLCHINOOK_CINT64             MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLCHINOOK_CSIZE_T            MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLCHINOOK_CFLOAT             MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLCHINOOK_CDOUBLE            MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLCHINOOK_CLONGDOUBLE        MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLCHINOOK_MLSHORT            MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLCHINOOK_MLINT              MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLCHINOOK_MLLONG             MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLCHINOOK_MLINT64            MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLCHINOOK_MLSIZE_T           MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#define MLCHINOOK_MLFLOAT            MLTK_LITTLEENDIAN_IEEE754_SINGLE
#define MLCHINOOK_MLDOUBLE           MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#define MLCHINOOK_MLLONGDOUBLE       MLTK_LITTLEENDIAN_IEEE754_DOUBLE



/********************  OLD_WIN_ENV  ****************************/
#define MLOLD_WIN_ENV_NUMERICS_ID    "Sep 13 1996, 13:46:34"
#define MLOLD_WIN_ENV_CSHORT         MLTK_CSHORT_P
#define MLOLD_WIN_ENV_CINT           MLTK_CINT_P
#define MLOLD_WIN_ENV_CLONG          MLTK_CLONG_P
#define MLOLD_WIN_ENV_CINT64         MLTK_64BIT_LITTLEENDIAN_STRUCTURE
#define MLOLD_WIN_ENV_CSIZE_T        MLTK_CLONG_P
#define MLOLD_WIN_ENV_CFLOAT         MLTK_CFLOAT_P
#define MLOLD_WIN_ENV_CDOUBLE        MLTK_CDOUBLE_P
#define MLOLD_WIN_ENV_CLONGDOUBLE    MLTK_CLONGDOUBLE_P
#define MLOLD_WIN_ENV_MLSHORT        MLTK_CSHORT_P
#define MLOLD_WIN_ENV_MLINT          MLTK_CINT_P
#define MLOLD_WIN_ENV_MLLONG         MLTK_CLONG_P
#define MLOLD_WIN_ENV_MLINT64        MLTK_64BIT_LITTLEENDIAN_STRUCTURE
#define MLOLD_WIN_ENV_MLSIZE_T       MLTK_CLONG_P
#define MLOLD_WIN_ENV_MLFLOAT        MLTK_CFLOAT_P
#define MLOLD_WIN_ENV_MLDOUBLE       MLTK_CDOUBLE_P
#define MLOLD_WIN_ENV_MLLONGDOUBLE   MLTK_CLONGDOUBLE_P


#define MLTK_CUCHAR  MLTK_8BIT_UNSIGNED_INTEGER
#define MLTK_MLUCHAR MLTK_8BIT_UNSIGNED_INTEGER

#if UNIX_MATHLINK /* hueristic that works for now */
	typedef unsigned int _uint32_nt;
	typedef signed int _sint32_nt;
#else
	typedef unsigned long _uint32_nt;
	typedef signed long _sint32_nt;
#endif


#if WINDOWS_MATHLINK
#define NEW_WIN32_NUMENV 1
#endif


/* #	define MATHLINK_NUMERICS_ENVIRONMENT_ID "Sep 16 1996, 23:14:20" M68KMACINTOSH_MATHLINK */
/* #	define MATHLINK_NUMERICS_ENVIRONMENT_ID_NUMB 33 */

/* #	define MATHLINK_NUMERICS_ENVIRONMENT_ID "newdog" POWERMACINTOSH_MATHLINK */
/* #	define MATHLINK_NUMERICS_ENVIRONMENT_ID_NUMB 24 */



#if WINDOWS_MATHLINK && NEW_WIN32_NUMENV
#if WIN32_MATHLINK
#define MATHLINK_NUMERICS_ENVIRONMENT_ID  MLBOERBOEL_NUMERICS_ID

#define MLTK_CSHORT        MLBOERBOEL_CSHORT
#define MLTK_CINT          MLBOERBOEL_CINT
#define MLTK_CLONG         MLBOERBOEL_CLONG
#define MLTK_CINT64        MLBOERBOEL_CINT64
#define MLTK_CSIZE_T       MLBOERBOEL_CSIZE_T
#define MLTK_CFLOAT        MLBOERBOEL_CFLOAT
#define MLTK_CDOUBLE       MLBOERBOEL_CDOUBLE
#define MLTK_CLONGDOUBLE   MLBOERBOEL_CLONGDOUBLE

#define MLTK_MLSHORT       MLBOERBOEL_MLSHORT
#define MLTK_MLINT         MLBOERBOEL_MLINT
#define MLTK_MLLONG        MLBOERBOEL_MLLONG
#define MLTK_MLINT64       MLBOERBOEL_MLINT64
#define MLTK_MLSIZE_T      MLBOERBOEL_MLSIZE_T
#define MLTK_MLFLOAT       MLBOERBOEL_MLFLOAT
#define MLTK_MLDOUBLE      MLBOERBOEL_MLDOUBLE
#define MLTK_MLLONGDOUBLE  MLBOERBOEL_MLLONGDOUBLE
#elif WIN64_MATHLINK
#define MATHLINK_NUMERICS_ENVIRONMENT_ID  MLBICHON_FRISE_NUMERICS_ID

#define MLTK_CSHORT        MLBICHON_FRISE_CSHORT
#define MLTK_CINT          MLBICHON_FRISE_CINT
#define MLTK_CLONG         MLBICHON_FRISE_CLONG
#define MLTK_CINT64        MLBICHON_FRISE_CINT64
#define MLTK_CSIZE_T       MLBICHON_FRISE_CSIZE_T
#define MLTK_CFLOAT        MLBICHON_FRISE_CFLOAT
#define MLTK_CDOUBLE       MLBICHON_FRISE_CDOUBLE
#define MLTK_CLONGDOUBLE   MLBICHON_FRISE_CLONGDOUBLE

#define MLTK_MLSHORT       MLBICHON_FRISE_MLSHORT
#define MLTK_MLINT         MLBICHON_FRISE_MLINT
#define MLTK_MLLONG        MLBICHON_FRISE_MLLONG
#define MLTK_MLINT64       MLBICHON_FRISE_MLINT64
#define MLTK_MLSIZE_T      MLBICHON_FRISE_MLSIZE_T
#define MLTK_MLFLOAT       MLBICHON_FRISE_MLFLOAT
#define MLTK_MLDOUBLE      MLBICHON_FRISE_MLDOUBLE
#define MLTK_MLLONGDOUBLE  MLBICHON_FRISE_MLLONGDOUBLE
#endif /* WINDOWS_MATHLINK */

#elif DARWIN_MATHLINK

#if X86_DARWIN_MATHLINK
#define MATHLINK_NUMERICS_ENVIRONMENT_ID MLBEAUCERON_NUMERICS_ID

#define MLTK_CSHORT        MLBEAUCERON_CSHORT
#define MLTK_CINT          MLBEAUCERON_CINT
#define MLTK_CLONG         MLBEAUCERON_CLONG
#define MLTK_CINT64        MLBEAUCERON_CINT64
#define MLTK_CSIZE_T       MLBEAUCERON_CSIZE_T
#define MLTK_CFLOAT        MLBEAUCERON_CFLOAT
#define MLTK_CDOUBLE       MLBEAUCERON_CDOUBLE
#define MLTK_CLONGDOUBLE   MLBEAUCERON_CLONGDOUBLE

#define MLTK_MLSHORT       MLBEAUCERON_MLSHORT
#define MLTK_MLINT         MLBEAUCERON_MLINT
#define MLTK_MLLONG        MLBEAUCERON_MLLONG
#define MLTK_MLINT64       MLBEAUCERON_MLINT64
#define MLTK_MLSIZE_T      MLBEAUCERON_MLSIZE_T
#define MLTK_MLFLOAT       MLBEAUCERON_MLFLOAT
#define MLTK_MLDOUBLE      MLBEAUCERON_MLDOUBLE
#define MLTK_MLLONGDOUBLE  MLBEAUCERON_MLLONGDOUBLE

#elif X86_64_DARWIN_MATHLINK
#define MATHLINK_NUMERICS_ENVIRONMENT_ID  MLGREATDANE_NUMERICS_ID

#define MLTK_CSHORT        MLGREATDANE_CSHORT
#define MLTK_CINT          MLGREATDANE_CINT
#define MLTK_CLONG         MLGREATDANE_CLONG
#define MLTK_CINT64        MLGREATDANE_CINT64
#define MLTK_CSIZE_T       MLGREATDANE_CSIZE_T
#define MLTK_CFLOAT        MLGREATDANE_CFLOAT
#define MLTK_CDOUBLE       MLGREATDANE_CDOUBLE
#define MLTK_CLONGDOUBLE   MLGREATDANE_CLONGDOUBLE

#define MLTK_MLSHORT       MLGREATDANE_MLSHORT
#define MLTK_MLINT         MLGREATDANE_MLINT
#define MLTK_MLLONG        MLGREATDANE_MLLONG
#define MLTK_MLINT64       MLGREATDANE_MLINT64
#define MLTK_MLSIZE_T      MLGREATDANE_MLSIZE_T
#define MLTK_MLFLOAT       MLGREATDANE_MLFLOAT
#define MLTK_MLDOUBLE      MLGREATDANE_MLDOUBLE
#define MLTK_MLLONGDOUBLE  MLGREATDANE_MLLONGDOUBLE

#elif ARM_DARWIN_MATHLINK
#define MATHLINK_NUMERICS_ENVIRONMENT_ID MLFRENCH_BULLDOG_NUMERICS_ID

#define MLTK_CSHORT        MLFRENCH_BULLDOG_CSHORT
#define MLTK_CINT          MLFRENCH_BULLDOG_CINT
#define MLTK_CLONG         MLFRENCH_BULLDOG_CLONG
#define MLTK_CINT64        MLFRENCH_BULLDOG_CINT64
#define MLTK_CSIZE_T       MLFRENCH_BULLDOG_CSIZE_T
#define MLTK_CFLOAT        MLFRENCH_BULLDOG_CFLOAT
#define MLTK_CDOUBLE       MLFRENCH_BULLDOG_CDOUBLE
#define MLTK_CLONGDOUBLE   MLFRENCH_BULLDOG_CLONGDOUBLE

#define MLTK_MLSHORT       MLFRENCH_BULLDOG_MLSHORT
#define MLTK_MLINT         MLFRENCH_BULLDOG_MLINT
#define MLTK_MLLONG        MLFRENCH_BULLDOG_MLLONG
#define MLTK_MLINT64       MLFRENCH_BULLDOG_MLINT64
#define MLTK_MLSIZE_T      MLFRENCH_BULLDOG_MLSIZE_T
#define MLTK_MLFLOAT       MLFRENCH_BULLDOG_MLFLOAT
#define MLTK_MLDOUBLE      MLFRENCH_BULLDOG_MLDOUBLE
#define MLTK_MLLONGDOUBLE  MLFRENCH_BULLDOG_MLLONGDOUBLE
#elif ARM64_DARWIN_MATHLINK
#define MATHLINK_NUMERICS_ENVIRONMENT_ID MLAKITAINU_NUMERICS_ID

#define MLTK_CSHORT        MLAKITAINU_CSHORT
#define MLTK_CINT          MLAKITAINU_CINT
#define MLTK_CLONG         MLAKITAINU_CLONG
#define MLTK_CINT64        MLAKITAINU_CINT64
#define MLTK_CSIZE_T       MLAKITAINU_CSIZE_T
#define MLTK_CFLOAT        MLAKITAINU_CFLOAT
#define MLTK_CDOUBLE       MLAKITAINU_CDOUBLE
#define MLTK_CLONGDOUBLE   MLAKITAINU_CLONGDOUBLE

#define MLTK_MLSHORT       MLAKITAINU_MLSHORT
#define MLTK_MLINT         MLAKITAINU_MLINT
#define MLTK_MLLONG        MLAKITAINU_MLLONG
#define MLTK_MLINT64       MLAKITAINU_MLINT64
#define MLTK_MLSIZE_T      MLAKITAINU_MLSIZE_T
#define MLTK_MLFLOAT       MLAKITAINU_MLFLOAT
#define MLTK_MLDOUBLE      MLAKITAINU_MLDOUBLE
#define MLTK_MLLONGDOUBLE  MLAKITAINU_MLLONGDOUBLE
#endif

#elif I86_LINUX_MATHLINK

#define MATHLINK_NUMERICS_ENVIRONMENT_ID  MLBERGAMASCO_NUMERICS_ID

#define MLTK_CSHORT        MLBERGAMASCO_CSHORT
#define MLTK_CINT          MLBERGAMASCO_CINT
#define MLTK_CLONG         MLBERGAMASCO_CLONG
#define MLTK_CINT64        MLBERGAMASCO_CINT64
#define MLTK_CSIZE_T       MLBERGAMASCO_CSIZE_T
#define MLTK_CFLOAT        MLBERGAMASCO_CFLOAT
#define MLTK_CDOUBLE       MLBERGAMASCO_CDOUBLE
#define MLTK_CLONGDOUBLE   MLBERGAMASCO_CLONGDOUBLE

#define MLTK_MLSHORT       MLBERGAMASCO_MLSHORT
#define MLTK_MLINT         MLBERGAMASCO_MLINT
#define MLTK_MLLONG        MLBERGAMASCO_MLLONG
#define MLTK_MLINT64       MLBERGAMASCO_MLINT64
#define MLTK_MLSIZE_T      MLBERGAMASCO_MLSIZE_T
#define MLTK_MLFLOAT       MLBERGAMASCO_MLFLOAT
#define MLTK_MLDOUBLE      MLBERGAMASCO_MLDOUBLE
#define MLTK_MLLONGDOUBLE  MLBERGAMASCO_MLLONGDOUBLE

#elif X86_64_LINUX_MATHLINK || ARM64_LINUX_MATHLINK
#define MATHLINK_NUMERICS_ENVIRONMENT_ID  MLBOXER_NUMERICS_ID

#define MLTK_CSHORT        MLBOXER_CSHORT
#define MLTK_CINT          MLBOXER_CINT
#define MLTK_CLONG         MLBOXER_CLONG
#define MLTK_CINT64        MLBOXER_CINT64
#define MLTK_CSIZE_T       MLBOXER_CSIZE_T
#define MLTK_CFLOAT        MLBOXER_CFLOAT
#define MLTK_CDOUBLE       MLBOXER_CDOUBLE
#define MLTK_CLONGDOUBLE   MLBOXER_CLONGDOUBLE

#define MLTK_MLSHORT       MLBOXER_MLSHORT
#define MLTK_MLINT         MLBOXER_MLINT
#define MLTK_MLLONG        MLBOXER_MLLONG
#define MLTK_MLINT64       MLBOXER_MLINT64
#define MLTK_MLSIZE_T      MLBOXER_MLSIZE_T
#define MLTK_MLFLOAT       MLBOXER_MLFLOAT
#define MLTK_MLDOUBLE      MLBOXER_MLDOUBLE
#define MLTK_MLLONGDOUBLE  MLBOXER_MLLONGDOUBLE

#elif ARM_LINUX_MATHLINK

#define MATHLINK_NUMERICS_ENVIRONMENT_ID  MLCHINOOK_NUMERICS_ID

#define MLTK_CSHORT        MLCHINOOK_CSHORT
#define MLTK_CINT          MLCHINOOK_CINT
#define MLTK_CLONG         MLCHINOOK_CLONG
#define MLTK_CINT64        MLCHINOOK_CINT64
#define MLTK_CSIZE_T       MLCHINOOK_CSIZE_T
#define MLTK_CFLOAT        MLCHINOOK_CFLOAT
#define MLTK_CDOUBLE       MLCHINOOK_CDOUBLE
#define MLTK_CLONGDOUBLE   MLCHINOOK_CLONGDOUBLE

#define MLTK_MLSHORT       MLCHINOOK_MLSHORT
#define MLTK_MLINT         MLCHINOOK_MLINT
#define MLTK_MLLONG        MLCHINOOK_MLLONG
#define MLTK_MLINT64       MLCHINOOK_MLINT64
#define MLTK_MLSIZE_T      MLCHINOOK_MLSIZE_T
#define MLTK_MLFLOAT       MLCHINOOK_MLFLOAT
#define MLTK_MLDOUBLE      MLCHINOOK_MLDOUBLE
#define MLTK_MLLONGDOUBLE  MLCHINOOK_MLLONGDOUBLE

#else
#define MATHLINK_NUMERICS_ENVIRONMENT_ID  MLOLD_WIN_ENV_NUMERICS_ID

#define MLTK_CSHORT        MLOLD_WIN_ENV_CSHORT
#define MLTK_CINT          MLOLD_WIN_ENV_CINT
#define MLTK_CLONG         MLOLD_WIN_ENV_CLONG
#define MLTK_CINT64        MLOLD_WIN_ENV_CINT64
#define MLTK_CSIZE_T       MLOLD_WIN_ENV_CSIZE_T
#define MLTK_CFLOAT        MLOLD_WIN_ENV_CFLOAT
#define MLTK_CDOUBLE       MLOLD_WIN_ENV_CDOUBLE
#define MLTK_CLONGDOUBLE   MLOLD_WIN_ENV_CLONGDOUBLE

#define MLTK_MLSHORT       MLOLD_WIN_ENV_MLSHORT
#define MLTK_MLINT         MLOLD_WIN_ENV_MLINT
#define MLTK_MLLONG        MLOLD_WIN_ENV_MLLONG
#define MLTK_MLINT64       MLOLD_WIN_ENV_MLINT64
#define MLTK_MLSIZE_T      MLOLD_WIN_ENV_MLSIZE_T
#define MLTK_MLFLOAT       MLOLD_WIN_ENV_MLFLOAT
#define MLTK_MLDOUBLE      MLOLD_WIN_ENV_MLDOUBLE
#define MLTK_MLLONGDOUBLE  MLOLD_WIN_ENV_MLLONGDOUBLE

#endif

/* Objects of these numeric types exist in MathLink only in the numerics
 * environment and, unfortunately, in the "stack frames" of the functions that
 * put atomic numbers like MLPutInteger.  These C types are used by client
 * programs solely for type-checked access to the BinaryNumber functions.
 */
typedef unsigned char uchar_nt;
typedef uchar_nt     * ucharp_nt;
typedef ucharp_nt    * ucharpp_nt;

typedef short              short_nt;
typedef short_nt     * shortp_nt;
typedef shortp_nt    * shortpp_nt;

typedef int                int_nt;
typedef int_nt       * intp_nt;
typedef intp_nt      * intpp_nt;

typedef long               long_nt;
typedef long_nt      * longp_nt;
typedef longp_nt     * longpp_nt;

typedef mlint64            int64_nt;
typedef int64_nt         * int64p_nt;
typedef int64p_nt        * int64pp_nt;

typedef float              float_nt;
typedef float_nt     * floatp_nt;
typedef floatp_nt    * floatpp_nt;

typedef double             double_nt;
typedef double_nt    * doublep_nt;
typedef doublep_nt   * doublepp_nt;

#ifndef I87EXTENDED_NT_DEFINED
#define I87EXTENDED_NT_DEFINED
struct _i87extended_nt { unsigned short w[5];};
#endif

#ifndef __extended_nt__
#if WINDOWS_MATHLINK && (MLTK_CLONGDOUBLE != MLTK_MLLONGDOUBLE) /* subtle predicate that works for old and new windows numenvs */
#define __extended_nt__ struct _i87extended_nt
#define EXTENDED_NT_TO_I87_EXTENDED(a,b) a = b
#define I87_EXTENDED_TO_EXTENDED_NT(a,b) a = b
#else
#define __extended_nt__ long double
#define EXTENDED_NT_TO_I87_EXTENDED(a,b) \
				{ \
					int i; \
					unsigned short *c = (unsigned short *)&b; \
					for(i = 0; i < 5; i++) a.w[i] = 0; \
					for(i = 1; i < 5; i++) a.w[i] = *(c + i); \
				}
#define I87_EXTENDED_TO_EXTENDED_NT(a,b) \
				{ \
					int i; \
					unsigned short *c = (unsigned short *)&a; \
					a = 0; \
					for(i = 1; i < 5; i++) *(c + i) = b.w[i]; \
				}
#endif
#endif

	typedef __extended_nt__    mlextended_double;

	typedef __extended_nt__    extended_nt;
	typedef extended_nt  * extendedp_nt;
	typedef extendedp_nt * extendedpp_nt;

#endif /* _MLNTYPES_H */



ML_EXTERN_C

#ifndef _MLSTDDEV_H
#define _MLSTDDEV_H








typedef void * dev_world;
typedef MLINK dev_cookie;

typedef dev_world * dev_worldp;
typedef dev_cookie * dev_cookiep;


typedef  MLAllocatorUPP dev_allocator;
#define call_dev_allocator CallMLAllocatorProc
#define new_dev_allocator NewMLAllocatorProc

typedef  MLDeallocatorUPP dev_deallocator;
#define call_dev_deallocator CallMLDeallocatorProc
#define new_dev_deallocator NewMLDeallocatorProc

typedef dev_main_type world_main_type;

#define MLSTDWORLD_INIT        16
#define MLSTDWORLD_DEINIT      17
#define MLSTDWORLD_MAKE        18

#if UNIX_MATHLINK
#define MLSTDWORLD_GET_SIGNAL_HANDLERS      29
#define MLSTDWORLD_RELEASE_SIGNAL_HANDLERS  30
#endif

#define MLSTDWORLD_PROTOCOL        31
#define MLSTDWORLD_MODES           32
#define MLSTDWORLD_STREAMCAPACITY  33
#define MLSTDWORLD_ID              34

#define MLSTDDEV_CONNECT_READY 19
#define MLSTDDEV_CONNECT       20
#define MLSTDDEV_DESTROY       21

#define MLSTDDEV_SET_YIELDER   22
#define MLSTDDEV_GET_YIELDER   23

#define MLSTDDEV_WRITE_MSG     24
#define MLSTDDEV_HAS_MSG       25
#define MLSTDDEV_READ_MSG      26
#define MLSTDDEV_SET_HANDLER   27
#define MLSTDDEV_GET_HANDLER   28


#define T_WORLD_INIT        MLSTDWORLD_INIT
#define T_WORLD_DEINIT      MLSTDWORLD_DEINIT
#define T_WORLD_MAKE        MLSTDWORLD_MAKE
#define T_DEV_CONNECT_READY MLSTDDEV_CONNECT_READY
#define T_DEV_CONNECT       MLSTDDEV_CONNECT
#define T_DEV_DESTROY       MLSTDDEV_DESTROY

#define T_DEV_SET_YIELDER   MLSTDDEV_SET_YIELDER
#define T_DEV_GET_YIELDER   MLSTDDEV_GET_YIELDER

#define T_DEV_WRITE_MSG     MLSTDDEV_WRITE_MSG
#define T_DEV_HAS_MSG       MLSTDDEV_HAS_MSG
#define T_DEV_READ_MSG      MLSTDDEV_READ_MSG
#define T_DEV_SET_HANDLER   MLSTDDEV_SET_HANDLER
#define T_DEV_GET_HANDLER   MLSTDDEV_GET_HANDLER


typedef unsigned long dev_mode;
/* edit here and in mathlink.r */
#define NOMODE           ((dev_mode)0x0000)
#define LOOPBACKBIT      ((dev_mode)0x0001)
#define LISTENBIT        ((dev_mode)0x0002)
#define CONNECTBIT       ((dev_mode)0x0004)
#define LAUNCHBIT        ((dev_mode)0x0008)
#define PARENTCONNECTBIT ((dev_mode)0x0010)
#define READBIT          ((dev_mode)0x0020)
#define WRITEBIT         ((dev_mode)0x0040)
#define SERVERBIT        ((dev_mode)0x0080)
#define ANYMODE          (~(dev_mode)0)

typedef dev_mode * dev_modep;





typedef unsigned long dev_options;

#define _DefaultOptions              ((dev_options)0x00000000)

#define _NetworkVisibleMask          ((dev_options)0x00000003)
#define _BrowseMask                  ((dev_options)0x00000010)
#define _NonBlockingMask             ((dev_options)0x00000020)
#define _InteractMask                ((dev_options)0x00000100)
#define _YieldMask                   ((dev_options)0x00000200)
#define _UseIPV6Mask                 ((dev_options)0x00010000)
#define _UseIPV4Mask                 ((dev_options)0x00020000)
#define _VersionMask                 ((dev_options)0x0F000000)
#define _UseNewTCPIPConnectionMask   ((dev_options)0x00100000)
#define _UseOldTCPIPConnectionMask   ((dev_options)0x00200000)
#define _UseUUIDTCPIPConnectionMask  ((dev_options)0x00000004)
#define _UseAnyNetworkAddressMask    ((dev_options)0x00000008)

#define _NetworkVisible              ((dev_options)0x00000000)
#define _LocallyVisible              ((dev_options)0x00000001)
#define _InternetVisible             ((dev_options)0x00000002)

#define _Browse                      ((dev_options)0x00000000)
#define _DontBrowse                  ((dev_options)0x00000010)

#define _NonBlocking                 ((dev_options)0x00000000)
#define _Blocking                    ((dev_options)0x00000020)

#define _Interact                    ((dev_options)0x00000000)
#define _DontInteract                ((dev_options)0x00000100)

#define _ForceYield                  ((dev_options)0x00000200)
#define _UseIPV6                     ((dev_options)0x00010000)
#define _UseIPV4                     ((dev_options)0x00020000)
#define _UseNewTCPIPConnection       ((dev_options)0x00100000)
#define _UseOldTCPIPConnection       ((dev_options)0x00200000)
#define _UseUUIDTCPIPConnection      ((dev_options)0x00000004)
#define _UseAnyNetworkAddress        ((dev_options)0x00000008)


/* DEVICE selector and WORLD selector masks */
#define INFO_MASK (1UL << 31)
#define INFO_TYPE_MASK ((1UL << 31) - 1UL)
#define INFO_SWITCH_MASK (1UL << 30)
#define MLDEVICE_MASK INFO_MASK
#define WORLD_MASK (INFO_MASK | (1UL << 30))

/* values returned by selector MLDEVICE_TYPE */
#define UNREGISTERED_TYPE  0
#define UNIXPIPE_TYPE      1
#define UNIXSOCKET_TYPE    2
#define LOOPBACK_TYPE      5
#define WINLOCAL_TYPE      9
#define WINFMAP_TYPE       10
#define WINSHM_TYPE        11
#define SOCKET2_TYPE       12
#define GENERIC_TYPE	   13  /* Internal use only, not valid for MLDeviceInformation */
#define UNIXSHM_TYPE       14
#define INTRAPROCESS_TYPE  15

/* selectors */
#define MLDEVICE_TYPE          MLDEVICE_MASK + 0UL                                        /* long */
#define MLDEVICE_NAME          MLDEVICE_MASK + 1UL                                        /* char */
#define MLDEVICE_NAME_SIZE     MLDEVICE_MASK + 2UL                                   /* long */
#define MLDEVICE_WORLD_ID      MLDEVICE_MASK + 5UL                                    /* char */
#define SHM_FD                 MLDEVICE_MASK + (UNIXSHM_TYPE * 256UL + 0UL)      /* int */
#define PIPE_FD                MLDEVICE_MASK + (UNIXPIPE_TYPE * 256UL + 0UL)     /* int */
#define PIPE_CHILD_PID         MLDEVICE_MASK + (UNIXPIPE_TYPE * 256UL + 1UL)     /* int */
#define SOCKET_FD              MLDEVICE_MASK + (UNIXSOCKET_TYPE * 256UL + 0UL)   /* int */
#define INTRA_FD               MLDEVICE_MASK + (INTRAPROCESS_TYPE * 256UL + 0UL) /* int */
#define SOCKET_PARTNER_ADDR    MLDEVICE_MASK + (UNIXSOCKET_TYPE * 256UL + 1UL)   /* unsigned long */
#define SOCKET_PARTNER_PORT    MLDEVICE_MASK + (UNIXSOCKET_TYPE * 256UL + 2UL)   /* unsigned short */
#define LOOPBACK_FD            MLDEVICE_MASK + (LOOPBACK_TYPE * 256UL + 0UL)     /* int */
#define INTRAPROCESS_FD        MLDEVICE_MASK + (INTRAPROCESS_TYPE * 256 + 0)     /* int */

#define	WINDOWS_SET_NOTIFY_WINDOW     MLDEVICE_MASK + 2330UL /* HWND */
#define	WINDOWS_REMOVE_NOTIFY_WINDOW  MLDEVICE_MASK + 2331UL /* HWND */
#define WINDOWS_READY_CONDITION       MLDEVICE_MASK + 2332UL /* HANDLE */

/* info selectors */
#define WORLD_THISLOCATION (1UL + WORLD_MASK)        /* char */
#define WORLD_MODES (2UL + WORLD_MASK)               /* dev_mode */
#define WORLD_PROTONAME (3UL + WORLD_MASK)           /* char */
#define WORLD_STREAMCAPACITY (4UL + WORLD_MASK)      /* long */ /*this belongs in mlolddev.h*/
#define WORLD_ID (5UL + WORLD_MASK)    /* char */


#ifndef MATHLINK_DEVICE_WORLD_ID
#define MATHLINK_DEVICE_WORLD_ID (__DATE__ ", " __TIME__)
#endif


#define MLDEVICE_MODE      MLDEVICE_MASK + 6UL                                    /* long */
#define MLDEVICE_OPTIONS   MLDEVICE_MASK + 7UL                                    /* long */



#define YIELDVERSION 1

typedef long devyield_result;
typedef long devyield_place;
typedef long devyield_count;
typedef unsigned long devyield_sleep;

#define INTERNAL_YIELDING 0
#define MAKE_YIELDING 1
#define CONNECT_YIELDING 2
#define READ_YIELDING 3
#define WRITE_YIELDING 4
#define DESTROY_YIELDING 5
#define READY_YIELDING 6


typedef struct MLYieldParams * MLYieldParameters;


#define MAX_SLEEP (600)
typedef struct MLYieldData{
	union {long l; double d; void * p;} private_data[8];
} * MLYieldDataPointer;

void MLNewYieldData( MLYieldDataPointer ydp   /* , dev_allocator, dev_deallocator */);
void MLFreeYieldData( MLYieldDataPointer ydp);
MLYieldParameters MLResetYieldData( MLYieldDataPointer ydp, devyield_place func_id);
int   MLSetYieldParameter( MLYieldParameters yp, unsigned long selector, void* data, unsigned long* len);
int   MLYieldParameter( MLYieldParameters yp, unsigned long selector, void* data, unsigned long* len);
devyield_sleep MLSetSleepYP( MLYieldParameters yp, devyield_sleep sleep);
devyield_count MLSetCountYP( MLYieldParameters yp, devyield_count count);


enum { MLSleepParameter = 1, MLCountParameter, MLPlaceParameter};





MLYPROC( int, MLYielderProcPtr, (MLINK mlp, MLYieldParameters yp));
typedef	MLYielderProcPtr MLDeviceYielderProcPtr;

typedef MLYielderProcPtr MLYielderUPP, MLDeviceYielderUPP;
#define NewMLYielderProc(userRoutine) (userRoutine)

#define NewMLDeviceYielderProc NewMLYielderProc

typedef  MLYielderUPP MLYieldFunctionType;

typedef MLYielderUPP MLYieldFunctionObject;

typedef  MLYieldFunctionObject dev_yielder;
typedef dev_yielder * dev_yielderp;


typedef unsigned long dev_message;
typedef dev_message * dev_messagep;


MLMPROC( void, MLHandlerProcPtr, (MLINK mlp, int m, int n));
typedef MLHandlerProcPtr MLDeviceHandlerProcPtr;


typedef MLHandlerProcPtr MLHandlerUPP, MLDeviceHandlerUPP;
#define NewMLHandlerProc(userRoutine) (userRoutine)

#define NewMLDeviceHandlerProc NewMLHandlerProc

typedef  MLHandlerUPP MLMessageHandlerType;

typedef MLHandlerUPP MLMessageHandlerObject;


typedef  MLMessageHandlerObject dev_msghandler;
typedef dev_msghandler * dev_msghandlerp;



#endif /* _MLSTDDEV_H */

#ifndef MLINTERFACE
#error
#endif

/* explicitly not protected by _MLSTDDEV_H in case MLDECL is redefined for multiple inclusion */

/*bugcheck //should the rest of YP stuff be exported? */
MLDECL( devyield_sleep,         MLSleepYP,               ( MLYieldParameters yp));
MLDECL( devyield_count,         MLCountYP,               ( MLYieldParameters yp));

MLDECL( MLYieldFunctionObject,  MLCreateYieldFunction,   ( MLEnvironment ep, MLYieldFunctionType yf, void* reserved)); /* reserved must be 0 */

MLDECL( MLYieldFunctionType,    MLDestroyYieldFunction,  ( MLYieldFunctionObject yfo));

MLDECL( int,        MLCallYieldFunction,     ( MLYieldFunctionObject yfo, MLINK mlp, MLYieldParameters p));

MLDECL( MLMessageHandlerObject, MLCreateMessageHandler,  ( MLEnvironment ep, MLMessageHandlerType mh, void* reserved)); /* reserved must be 0 */


MLDECL( MLMessageHandlerType,   MLDestroyMessageHandler, ( MLMessageHandlerObject mho));

MLDECL( void,                   MLCallMessageHandler,    ( MLMessageHandlerObject mho, MLINK mlp, int m, int n));


/* just some type-safe casts */
MLDECL( __MLProcPtr__, MLYielderCast, ( MLYielderProcPtr yp));
MLDECL( __MLProcPtr__, MLHandlerCast, ( MLHandlerProcPtr mh));

ML_END_EXTERN_C






#ifndef MLSIGNAL_H
#define MLSIGNAL_H

MLYPROC( void, MLSigHandlerProcPtr, (int signal));

typedef MLSigHandlerProcPtr MLSignalHandlerType;
typedef void * MLSignalHandlerObject;

#endif /* MLSIGNAL_H */





#ifndef _MLMAKE_H
#define _MLMAKE_H










/* --binding layer-- */
/*************** Starting MathLink ***************/

#define MLPARAMETERSIZE_R1 256
#define MLPARAMETERSIZE 356

typedef char * MLParametersPointer;
typedef char MLParameters[MLPARAMETERSIZE];

#define MLLoopBackOpen MLLoopbackOpen



ML_EXTERN_C
MLUPROC( void, MLUserProcPtr, (MLINK));
ML_END_EXTERN_C

typedef MLUserProcPtr MLUserUPP;
#define NewMLUserProc(userRoutine) (userRoutine)

typedef MLUserUPP MLUserFunctionType;
typedef MLUserFunctionType * MLUserFunctionTypePointer;

typedef MLUserUPP MLUserFunction;


/* The following defines are
 * currently for internal use only.
 */


/* edit here and in mldevice.h and mathlink.r */
#define MLNetworkVisibleMask         ((unsigned long)0x00000003)          /* 00000000000000000000011 */
#define MLBrowseMask                 ((unsigned long)0x00000010)          /* 00000000000000000010000 */
#define MLNonBlockingMask            ((unsigned long)0x00000020)          /* 00000000000000000110000 */
#define MLInteractMask               ((unsigned long)0x00000100)          /* 00000000000000100000000 */
#define MLYieldMask                  ((unsigned long)0x00000200)          /* 00000000000001000000000 */
#define MLUseIPV6Mask                ((unsigned long)0x00010000)          /* 00000010000000000000000 */
#define MLUseIPV4Mask                ((unsigned long)0x00020000)          /* 00000100000000000000000 */
#define MLVersionMask                ((unsigned long)0x0000F000)          /* 00000001111000000000000 */
#define MLUseNewTCPIPConnectionMask  ((unsigned long)0x00100000)          /* 00100000000000000000000 */
#define MLUseOldTCPIPConnectionMask  ((unsigned long)0x00200000)          /* 01000000000000000000000 */
#define MLUseUUIDTCPIPConnectionMask ((unsigned long)0x00000004)          /* 00000000000000000000110 */
#define MLUseAnyNetworkAddressMask   ((unsigned long)0x00000008)          /* 00000000000000000001000 */

                                                                          /* 01100111111001100111111 */

#define MLDefaultOptions             ((unsigned long)0x00000000)
#define MLNetworkVisible             ((unsigned long)0x00000000)
#define MLLocallyVisible             ((unsigned long)0x00000001)
#define MLInternetVisible            ((unsigned long)0x00000002)

#define MLBrowse                     ((unsigned long)0x00000000)
#define MLDontBrowse                 ((unsigned long)0x00000010)

#define MLNonBlocking                ((unsigned long)0x00000000)
#define MLBlocking                   ((unsigned long)0x00000020)

#define MLInteract                   ((unsigned long)0x00000000)
#define MLDontInteract               ((unsigned long)0x00000100)

#define MLForceYield                 ((unsigned long)0x00000200)
#define MLUseIPV6                    ((unsigned long)0x00010000)
#define MLUseIPV4                    ((unsigned long)0x00020000)

#define MLUseNewTCPIPConnection      ((unsigned long)0x00100000)
#define MLUseOldTCPIPConnection      ((unsigned long)0x00200000)
#define MLUseUUIDTCPIPConnection     ((unsigned long)0x00000004)

#define MLUseAnyNetworkAddress       ((unsigned long)0x00000008)

/* Encoding types for use with MLSetEncodingParameter */
#define MLASCII_ENC		1
#define MLBYTES_ENC		2
#define MLUCS2_ENC		3
#define MLOLD_ENC		4
#define MLUTF8_ENC		5
#define MLUTF16_ENC		6
#define MLUTF32_ENC		8

#define MLTOTAL_TEXT_ENCODINGS 8

#define MLLOGERROR              0
#define MLLOGWARNING            1
#define MLLOGNOTICE             2
#define MLLOGINFO               3
#define MLLOGDEBUG              4
#define MLLOGDEBUG1             5
#define MLLOGDEBUG2             6
#define MLLOGDEBUG3             7
#define MLLOGDEBUG4             8

#endif /* _MLMAKE_H */


/* explicitly not protected by _MLMAKE_H in case MLDECL is redefined for multiple inclusion */


ML_EXTERN_C

MLDECL( MLEnvironmentParameter, MLNewParameters, (unsigned long rev, unsigned long apirev));
MLDECL( void,                   MLReleaseParameters, (MLEnvironmentParameter ep));
MLDECL( void, MLSetAllocParameter, (MLEnvironmentParameter ep, MLAllocator allocator, MLDeallocator deallocator));
MLDECL( long, MLSetThreadSafeLinksParameter, (MLEnvironmentParameter ep));
MLDECL( int,  MLAllocParameter,       (MLEnvironmentParameter ep, MLAllocator* allocator, MLDeallocator* deallocator));
MLDECL( long, MLSetResourceParameter, (MLEnvironmentParameter ep, const char *path));
MLDECL( long, MLSetDeviceParameter,   (MLEnvironmentParameter ep, const char *devspec));
MLDECL( long, MLErrorParameter,       (MLEnvironmentParameter ep));
MLDECL( long, MLSetEncodingParameter, (MLEnvironmentParameter ep, unsigned int etype));
MLDECL( long, MLDoNotHandleSignalParameter,    (MLEnvironmentParameter ep, int signum));




MLDECL( void,          MLStopHandlingSignal, (MLEnvironment env, int signum));
MLDECL( void,          MLHandleSignal,       (MLEnvironment env, int signum));

MLDECL( long,          MLSetEnvironmentData,           ( MLEnvironment env, void *cookie));
MLDECL( void *,        MLEnvironmentData,              ( MLEnvironment env));
MLDECL( int,           MLSetSignalHandler,             ( MLEnvironment env, int signum, void *so));
MLDECL( int,           MLSetSignalHandlerFromFunction, ( MLEnvironment env, int signum, MLSignalHandlerType sigfunc));
MLDECL( int,           MLUnsetSignalHandler,           ( MLEnvironment env, int signum, MLSignalHandlerType sigfunc));

MLDECL( long,          MLSetSymbolReplacement,         ( MLINK mlp, const char *priv, int prlen, const char *pub, int pblen));
MLDECL( int,           MLClearSymbolReplacement,       ( MLINK mlp, long index));
MLDECL( void,          MLClearAllSymbolReplacements,   ( MLINK mlp));


MLDECL(MLEnvironment,  MLInitialize,   ( MLEnvironmentParameter ep));


MLDECL( void,          MLDeinitialize, ( MLEnvironment env));

/*************** MathLink Revsion Number/Interface Number ************/

MLDECL( void,          MLVersionNumbers, ( MLEnvironment env, int *inumb, int *rnumb, int *bnumb));

MLDECL( int,               MLCompilerID, (MLEnvironment env, const char **id));
MLDECL( void,       MLReleaseCompilerID, (MLEnvironment env, const char *id));

MLDECL( int,           MLUCS2CompilerID, (MLEnvironment env, unsigned short **id, int *length));
MLDECL( void,   MLReleaseUCS2CompilerID, (MLEnvironment env, unsigned short *id, int length));

MLDECL( int,           MLUTF8CompilerID, (MLEnvironment env, unsigned char **id, int *length));
MLDECL( void,   MLReleaseUTF8CompilerID, (MLEnvironment env, unsigned char *id, int length));

MLDECL( int,          MLUTF16CompilerID, (MLEnvironment env, unsigned short **id, int *length));
MLDECL( void,  MLReleaseUTF16CompilerID, (MLEnvironment env, unsigned short *id, int length));

MLDECL( int,          MLUTF32CompilerID, (MLEnvironment env, unsigned int **id, int *length));
MLDECL( void,  MLReleaseUTF32CompilerID, (MLEnvironment env, unsigned int *id, int length));

/********************************************************************/

/* or, if you use MLOpenArgcArgv, ...*/

MLDECL( MLEnvironment, MLBegin, (MLEnvironmentParameter ep));


MLDECL( void,          MLEnd,   ( MLEnvironment env));

/*************** Environment Identification Interface ***************/

MLDECL( int, MLSetEnvIDString, ( MLEnvironment ep, const char *environment_id)); /* APPIDSERV */
MLDECL( int, MLGetLinkedEnvIDString, (MLINK mlp, const char **environment_id)); /* APPIDSERV */
MLDECL( void, MLReleaseEnvIDString, (MLINK mlp, const char *environment_id));


/*********************************************************************/


/**************** Network Interface List API *************************/
MLDECL( char **,    MLGetNetworkAddressList, ( MLEnvironment ep, unsigned long *size ));
MLDECL( void,   MLReleaseNetworkAddressList, ( MLEnvironment ep, char **addresses, unsigned long size));


MLDECL( char **,        MLGetDomainNameList, ( MLEnvironment ep, unsigned long *size ));
MLDECL( void,       MLReleaseDomainNameList, ( MLEnvironment ep, char **dnsnames, unsigned long size));


/*********************************************************************/


/************************* Runtime Device Inspection API ***************************/
MLDECL(int, MLGetAvailableLinkProtocolNames, (MLEnvironment ep, char ***protocolNames, int *length));
MLDECL(void,     MLReleaseLinkProtocolNames, (MLEnvironment ep, char **protocolNames, int length));
/*********************************************************************/


/************************* Enumerate Open Links in an Env *************/
MLDECL(int,       MLGetLinksFromEnvironment, (MLEnvironment ep, MLINK **links, int *length));
MLDECL(void,  MLReleaseLinksFromEnvironment, (MLEnvironment ep, MLINK *links, int length));
/*********************************************************************/



#if MLNTESTPOINTS < 1
#undef MLNTESTPOINTS
#define MLNTESTPOINTS 1
#endif


#if MLNTESTPOINTS < 2
#undef MLNTESTPOINTS
#define MLNTESTPOINTS 2
#endif


#if MLNTESTPOINTS < 3
#undef MLNTESTPOINTS
#define MLNTESTPOINTS 3
#endif


#if MLNTESTPOINTS < 4
#undef MLNTESTPOINTS
#define MLNTESTPOINTS 4
#endif


MLDECL( long, MLNumericsQuery, ( MLEnvironment ep, unsigned long selector, void *p1, void *p2, long *np));


/*************** Connection interface ***************/




MLDECL( int,                            MLValid, ( MLINK mlp));





MLDECL( char **,       MLFilterArgv,   ( MLEnvironment ep, char **argv, char **argv_end));


MLDECL( long,          MLFeatureString, ( MLINK mlp, char *buf, long buffsize));
MLDECL( MLINK,         MLOpenArgv,      ( MLEnvironment ep, char **argv, char **argv_end, int *errp));
MLDECL( MLINK,         MLOpenArgcArgv,  ( MLEnvironment ep, int argc, char **argv, int *errp));
MLDECL( MLINK,         MLOpenString,    ( MLEnvironment ep, const char *command_line, int *errp));
MLDECL( MLINK,         MLLoopbackOpen,  ( MLEnvironment ep, int *errp));
MLDECL( int,           MLStringToArgv,  ( const char *commandline, char *buf, char **argv, int len));
MLDECL( long,          MLScanString,    ( char **argv, char ***argv_end, char **commandline, char **buf));
MLDECL( long,          MLPrintArgv,     ( char *buf, char **buf_endp, char ***argvp, char **argv_end));

MLDECL( const char *,     MLErrorMessage,  ( MLINK mlp));
MLDECL( const char *,     MLErrorString,   ( MLEnvironment env, long err));

MLDECL( const unsigned short *,  MLUCS2ErrorMessage,  (MLINK mlp, int *length));
MLDECL( const unsigned char *,   MLUTF8ErrorMessage,  (MLINK mlp, int *length));
MLDECL( const unsigned short *,  MLUTF16ErrorMessage, (MLINK mlp, int *length));
MLDECL( const unsigned int *,    MLUTF32ErrorMessage, (MLINK mlp, int *length));

MLDECL(void,  MLReleaseErrorMessage,      (MLINK mlp, const char *message));
MLDECL(void,  MLReleaseUCS2ErrorMessage,  (MLINK mlp, const unsigned short *message, int length));
MLDECL(void,  MLReleaseUTF8ErrorMessage,  (MLINK mlp, const unsigned char *message, int length));
MLDECL(void,  MLReleaseUTF16ErrorMessage, (MLINK mlp, const unsigned short *message, int length));
MLDECL(void,  MLReleaseUTF32ErrorMessage, (MLINK mlp, const unsigned int *message, int length));

MLDECL( MLINK,         MLOpen,          ( int argc, char **argv));
MLDECL( MLINK,         MLOpenInEnv,     ( MLEnvironment env, int argc, char **argv, int *errp));


MLDECL( MLINK,         MLDuplicateLink,   ( MLINK parentmlp, const char *name, int *errp ));

MLDECL( int,  MLConnect,         ( MLINK mlp));
MLDECL( int,  MLActivate,        ( MLINK mlp));

#ifndef __feature_setp__
#define __feature_setp__
typedef struct feature_set* feature_setp;
#endif
MLDECL( int,  MLEstablish,       ( MLINK mlp, feature_setp features));

MLDECL( int,  MLEstablishString, ( MLINK mlp, const char *features));

MLDECL( void,          MLClose,           ( MLINK mlp));

MLDECL( void,          MLSetUserData,   ( MLINK mlp, void* data, MLUserFunction f));
MLDECL( void*,         MLUserData,      ( MLINK mlp, MLUserFunctionType *fp));
MLDECL( void,          MLSetUserBlock,  ( MLINK mlp, void* userblock));
MLDECL( void*,         MLUserBlock,     ( MLINK mlp));

/* just a type-safe cast */
MLDECL( __MLProcPtr__, MLUserCast, ( MLUserProcPtr f));


MLDECL(int,            MLLogStreamToFile, (MLINK mlp, const char *name));
MLDECL(int,       MLDisableLoggingStream, (MLINK mlp));
MLDECL(int,        MLEnableLoggingStream, (MLINK mlp));
MLDECL(int,    MLStopLoggingStreamToFile, (MLINK mlp, const char *name));
MLDECL(int,          MLStopLoggingStream, (MLINK mlp));

MLDECL(int,         MLLogFileNameForLink, (MLINK mlp, const char **name));
MLDECL(void, MLReleaseLogFileNameForLink, (MLINK mlp, const char *name));

/* MLLinkName returns a pointer to the link's name.
 * Links are generally named when they are created
 * and are based on information that is potentially
 * useful and is available at that time.
 * Do not attempt to deallocate the name's storage
 * through this pointer.  The storage should be
 * considered in read-only memory.
 */

MLDECL( const char *, MLName,    ( MLINK mlp));
MLDECL( const char *, MLLinkName,    ( MLINK mlp));

MLDECL( const unsigned short *, MLUCS2LinkName,  (MLINK mlp, int *length));
MLDECL( const unsigned char *,  MLUTF8LinkName,  (MLINK mlp, int *length));
MLDECL( const unsigned short *, MLUTF16LinkName, (MLINK mlp, int *length));
MLDECL( const unsigned int *,   MLUTF32LinkName, (MLINK mlp, int *length));

MLDECL(void, MLReleaseLinkName,      (MLINK mlp, const char *name));
MLDECL(void, MLReleaseUCS2LinkName,  (MLINK mlp, const unsigned short *name, int length));
MLDECL(void, MLReleaseUTF8LinkName,  (MLINK mlp, const unsigned char *name, int length));
MLDECL(void, MLReleaseUTF16LinkName, (MLINK mlp, const unsigned short *name, int length));
MLDECL(void, MLReleaseUTF32LinkName, (MLINK mlp, const unsigned int *name, int length));

MLDECL( long,      MLNumber,  ( MLINK mlp));
MLDECL( long,  MLToLinkID,  ( MLINK mlp));
MLDECL( MLINK, MLFromLinkID, ( MLEnvironment ep, long n));

MLDECL( char *,  MLSetName, ( MLINK mlp, const char *name));

/* The following functions are
 * currently for internal use only.
 */

MLDECL( void*, MLInit,   ( MLallocator alloc, MLdeallocator dealloc, void* enclosing_environment));
MLDECL( void,  MLDeinit, ( void* env));
MLDECL( void*, MLEnclosingEnvironment, ( void* ep));
MLDECL( void*, MLinkEnvironment, ( MLINK mlp));


MLDECL(void, MLEnableLinkLock,  (MLINK mlp));
MLDECL(void, MLDisableLinkLock, (MLINK mlp));


MLDECL( MLEnvironment, MLLinkEnvironment, (MLINK mlp));


MLDECL( int,  MLIsLinkLoopback, (MLINK mlp));

/* the following two functions are for internal use only */
MLDECL( MLYieldFunctionObject, MLDefaultYieldFunction,    ( MLEnvironment env));

MLDECL( int,          MLSetDefaultYieldFunction, ( MLEnvironment env, MLYieldFunctionObject yf));


ML_END_EXTERN_C





#ifndef MLLINKSERVER_H
#define MLLINKSERVER_H






ML_EXTERN_C

typedef void * MLLinkServer;

typedef void (*MLNewLinkCallbackFunction)(MLLinkServer server, MLINK link);

MLDECL(MLLinkServer, MLNewLinkServer, (MLEnvironment env, void *context, int *error));

MLDECL(MLLinkServer, MLNewLinkServerWithPort, (MLEnvironment env, unsigned short port, void *context,
    int *error));

MLDECL(MLLinkServer, MLNewLinkServerWithPortAndInterface, (MLEnvironment env, unsigned short port, const char *iface,
    void *context, int *error));

MLDECL(void, MLShutdownLinkServer, (MLLinkServer server));

MLDECL(void, MLRegisterCallbackFunctionWithLinkServer, (MLLinkServer server, MLNewLinkCallbackFunction function));

MLDECL(MLINK, MLWaitForNewLinkFromLinkServer, (MLLinkServer server, int *error));

MLDECL(unsigned short, MLPortFromLinkServer, (MLLinkServer server, int *error));

MLDECL(const char *, MLInterfaceFromLinkServer, (MLLinkServer server, int *error));

MLDECL(void *, MLContextFromLinkServer, (MLLinkServer server, int *error));

MLDECL(void, MLReleaseInterfaceFromLinkServer, (MLLinkServer server, const char *iface));

ML_END_EXTERN_C

#endif /* MLLINKSERVER_H */




#ifndef MLSERVICEDISCOVERYAPI_H
#define MLSERVICEDISCOVERYAPI_H







ML_EXTERN_C

#define MLSDADDSERVICE      0x0001
#define MLSDREMOVESERVICE   0x0002
#define MLSDBROWSEERROR     0x0003
#define MLSDRESOLVEERROR    0x0004
#define MLSDREGISTERERROR   0x0005
#define MLSDMORECOMING      0x0010
#define MLSDNAMECONFLICT    0x0007

typedef void * MLServiceRef;

typedef void (*MLBrowseCallbackFunction)(MLEnvironment env, MLServiceRef ref, int flag,
	const char *serviceName, void *context);

MLDECL(int,     MLBrowseForLinkServices, (MLEnvironment env,
    MLBrowseCallbackFunction callbackFunction, const char *serviceProtocol,
    const char *domain, void *context, MLServiceRef *ref));

MLDECL(void, MLStopBrowsingForLinkServices, (MLEnvironment env, MLServiceRef ref));

typedef void (*MLResolveCallbackFunction)(MLEnvironment env, MLServiceRef ref, const char *serviceName,
	const char *linkName, const char *protocol, int options, void *context);

MLDECL(int, MLResolveLinkService, (MLEnvironment env,
    MLResolveCallbackFunction, const char *serviceProtocol,
    const char *serviceName, void *context, MLServiceRef *ref));

MLDECL(void, MLStopResolvingLinkService, (MLEnvironment env, MLServiceRef ref));

typedef void (*MLRegisterCallbackFunction)(MLEnvironment env, MLServiceRef ref, int flag, const char *serviceName,
	void *context);

MLDECL(MLINK, MLRegisterLinkServiceWithPortAndHostname, (MLEnvironment env, const char *serviceProtocol,
    const char *serviceName, unsigned short port, const char *hostname, MLRegisterCallbackFunction function,
    const char *domain, void *context, MLServiceRef *ref, int *error));

MLDECL(MLINK, MLRegisterLinkServiceWithHostname, (MLEnvironment env, const char *serviceProtocol,
    const char *serviceName, const char *hostname, MLRegisterCallbackFunction function,
    const char *domain, void *context, MLServiceRef *ref, int *error));

MLDECL(MLINK, MLRegisterLinkService, (MLEnvironment env, const char *serviceProtocol,
    const char *serviceName, MLRegisterCallbackFunction function,
    const char *domain, void *context, MLServiceRef *, int *error));

MLDECL(MLINK, MLRegisterLinkServiceUsingLinkProtocol, (MLEnvironment env, const char *serviceProtocol,
	const char *serviceName, unsigned short port, const char *hostname, const char *protocol,
	MLRegisterCallbackFunction function, const char *domain, void *context, MLServiceRef *ref, int *error));

MLDECL(void, MLRegisterLinkServiceFromLinkServer, (MLEnvironment env, const char *serviceProtocol,
    const char *serviceName, MLLinkServer server, MLRegisterCallbackFunction function, const char *domain,
    void *context, MLServiceRef *ref, int *error));

MLDECL(void, MLStopRegisteringLinkService, (MLEnvironment env, MLServiceRef ref));

MLDECL(void, MLStopRegisteringLinkServiceForLink, (MLEnvironment env, MLINK link, MLServiceRef ref));

MLDECL(const char *, MLServiceProtocolFromReference, (MLEnvironment env, MLServiceRef ref));

ML_END_EXTERN_C

#endif /* end of include guard: MLSERVICEDISCOVERYAPI_H */



#ifndef _MLERRORS_H
#define _MLERRORS_H



/*************** MathLink errors ***************/
/*
 * When some problem is detected within MathLink, routines
 * will return a simple indication of failure and store
 * an error code internally. (For routines that have nothing
 * else useful to return, success is indicated by returning
 * non-zero and failure by returning 0.)  MLerror() returns
 * the current error code;  MLErrorMessage returns an English
 * language description of the error.
 * The error MLEDEAD is irrecoverable.  For the others, MLClearError()
 * will reset the error code to MLEOK.
 */



#ifndef _MLERRNO_H
#define _MLERRNO_H

/* edit here and in mlerrstr.h */

#define MLEUNKNOWN          -1
#define MLEOK                0
#define MLEDEAD              1
#define MLEGBAD              2
#define MLEGSEQ              3
#define MLEPBTK              4
#define MLEPSEQ              5
#define MLEPBIG              6
#define MLEOVFL              7
#define MLEMEM               8
#define MLEACCEPT            9
#define MLECONNECT          10
#define MLECLOSED           11
#define MLEDEPTH            12  /* internal error */
#define MLENODUPFCN         13  /* stream cannot be duplicated */

#define MLENOACK            15  /* */
#define MLENODATA           16  /* */
#define MLENOTDELIVERED     17  /* */
#define MLENOMSG            18  /* */
#define MLEFAILED           19  /* */

#define MLEGETENDEXPR       20
#define MLEPUTENDPACKET     21 /* unexpected call of MLEndPacket */
                               /* currently atoms aren't
                                * counted on the way out so this error is raised only when
                                * MLEndPacket is called in the midst of an atom
                                */
#define MLENEXTPACKET       22
#define MLEUNKNOWNPACKET    23
#define MLEGETENDPACKET     24
#define MLEABORT            25
#define MLEMORE             26 /* internal error */
#define MLENEWLIB           27
#define MLEOLDLIB           28
#define MLEBADPARAM         29
#define MLENOTIMPLEMENTED   30


#define MLEINIT             32
#define MLEARGV             33
#define MLEPROTOCOL         34
#define MLEMODE             35
#define MLELAUNCH           36
#define MLELAUNCHAGAIN      37
#define MLELAUNCHSPACE      38
#define MLENOPARENT         39
#define MLENAMETAKEN        40
#define MLENOLISTEN         41
#define MLEBADNAME          42
#define MLEBADHOST          43
#define MLERESOURCE         44  /* a required resource was missing */
#define MLELAUNCHFAILED     45
#define MLELAUNCHNAME       46
#define MLEPDATABAD         47
#define MLEPSCONVERT        48
#define MLEGSCONVERT        49
#define MLENOTEXE           50
#define MLESYNCOBJECTMAKE   51
#define MLEBACKOUT          52
#define MLEBADOPTSYM        53
#define MLEBADOPTSTR        54
#define MLENEEDBIGGERBUFFER 55
#define MLEBADNUMERICSID    56
#define MLESERVICENOTAVAILABLE 57
#define MLEBADARGUMENT      58
#define MLEBADDISCOVERYHOSTNAME         59
#define MLEBADDISCOVERYDOMAINNAME       60
#define MLEBADSERVICENAME               61
#define MLEBADDISCOVERYSTATE            62
#define MLEBADDISCOVERYFLAGS            63
#define MLEDISCOVERYNAMECOLLISION       64
#define MLEBADSERVICEDISCOVERY          65
#define MLELAST MLESERVICENOTAVAILABLE

#define MLETRACEON         996  /* */
#define MLETRACEOFF        997  /* */
#define MLEDEBUG           998  /* */
#define MLEASSERT          999  /* an internal assertion failed */
#define MLEUSER           1000  /* start of user defined errors */


#endif /* _MLERRNO_H */



#endif /* _MLERRORS_H */

/* explicitly not protected by _MLERRORS_H in case MLDECL is redefined for multiple inclusion */

ML_EXTERN_C
MLDECL( int,  MLError,        ( MLINK mlp));
MLDECL( int,  MLClearError,   ( MLINK mlp));
MLDECL( int,  MLSetError,     ( MLINK mlp, int err));
ML_END_EXTERN_C



#ifndef _MLYLDMSG_H
#define _MLYLDMSG_H








enum {	MLTerminateMessage = 1, MLInterruptMessage, MLAbortMessage,
	MLEndPacketMessage, MLSynchronizeMessage, MLImDyingMessage,
	MLWaitingAcknowledgment, MLMarkTopLevelMessage, MLLinkClosingMessage,
	MLAuthenticateFailure, MLSuspendActivitiesMessage, MLResumeActivitiesMessage,
	MLFirstUserMessage = 128, MLLastUserMessage = 255 };

typedef unsigned long devinfo_selector;


#endif /* _MLYLDMSG_H */

/* explicitly not protected by _MLYLDMSG_H in case MLDECL is redefined for multiple inclusion */

ML_EXTERN_C

#ifndef MLINTERFACE
#error
#endif

MLDECL( int,   MLPutMessage,   ( MLINK mlp, int  msg));
MLDECL( int,   MLGetMessage,   ( MLINK mlp, int *mp, int *np));
MLDECL( int,   MLMessageReady, ( MLINK mlp));

MLDECL( int,   MLPutMessageWithArg, ( MLINK mlp, int msg, int arg));


MLDECL( MLMessageHandlerObject, MLGetMessageHandler,    ( MLINK mlp));
MLDECL( MLMessageHandlerObject, MLMessageHandler,    ( MLINK mlp));

MLDECL( MLYieldFunctionObject,  MLGetYieldFunction,     ( MLINK mlp));
MLDECL( MLYieldFunctionObject,  MLYieldFunction,     ( MLINK mlp));

MLDECL( int,  MLSetMessageHandler, ( MLINK mlp, MLMessageHandlerObject h));
MLDECL( int,  MLSetYieldFunction,  ( MLINK mlp, MLYieldFunctionObject yf));


MLDECL( int, MLDeviceInformation, ( MLINK mlp, devinfo_selector selector, void* buf, long *buflen));

MLDECL( int,         MLLowLevelDeviceName, (MLINK mlp, const char **name));
MLDECL( void, MLReleaseLowLevelDeviceName, (MLINK mlp, const char *name));


ML_END_EXTERN_C


/*************** Textual interface ***************/


#ifndef _MLGET_H
#define _MLGET_H




#endif /* _MLGET_H */

/* explicitly not protected by _MLGET_H in case MLDECL is redefined for multiple inclusion */

ML_EXTERN_C

MLDECL( int,   MLGetNext,          ( MLINK mlp));
MLDECL( int,   MLGetNextRaw,       ( MLINK mlp));
MLDECL( int,   MLGetType,          ( MLINK mlp));
MLDECL( int,   MLGetRawType,       ( MLINK mlp));
MLDECL( int,   MLGetRawData,       ( MLINK mlp, unsigned char *data, int size, int *gotp));
MLDECL( int,   MLGetData,          ( MLINK mlp, char *data, int size, int *gotp));
MLDECL( int,   MLGetArgCount,      ( MLINK mlp, int *countp));
MLDECL( int,   MLGetRawArgCount,   ( MLINK mlp, int *countp));
MLDECL( int,   MLBytesToGet,       ( MLINK mlp, int *leftp));
MLDECL( int,   MLRawBytesToGet,    ( MLINK mlp, int *leftp));
MLDECL( int,   MLExpressionsToGet, ( MLINK mlp, int *countp));

MLDECL( int,   MLNewPacket,        ( MLINK mlp));

MLDECL( int,   MLTakeLast,         ( MLINK mlp, int eleft));

MLDECL( int,   MLFill,             ( MLINK mlp));
ML_END_EXTERN_C



#ifndef _MLPUT_H
#define _MLPUT_H





#define MLPutExpression is obsolete, use MLPutComposite

#endif /* _MLPUT_H */

/* explicitly not protected by _MLPUT_H in case MLDECL is redefined for multiple inclusion */

ML_EXTERN_C

MLDECL( int,   MLPutNext,      ( MLINK mlp, int tok));
MLDECL( int,   MLPutType,      ( MLINK mlp, int tok));
MLDECL( int,   MLPutRawSize,   ( MLINK mlp, int size));
MLDECL( int,   MLPutRawData,   ( MLINK mlp, const unsigned char *data, int len));
MLDECL( int,   MLPutArgCount,  ( MLINK mlp, int argc));
MLDECL( int,   MLPutComposite, ( MLINK mlp, int argc));
MLDECL( int,   MLBytesToPut,   ( MLINK mlp, int *leftp));
MLDECL( int,   MLEndPacket,    ( MLINK mlp));
MLDECL( int,   MLFlush,        ( MLINK mlp));

ML_END_EXTERN_C



#ifndef _MLTK_H
#define _MLTK_H


#define	MLTKOLDINT     'I'		/* 73 Ox49 01001001 */ /* integer leaf node */
#define	MLTKOLDREAL    'R'		/* 82 Ox52 01010010 */ /* real leaf node */


#define	MLTKFUNC    'F'		/* 70 Ox46 01000110 */ /* non-leaf node */

#define	MLTKERROR   (0)		/* bad token */
#define	MLTKERR     (0)		/* bad token */

/* text token bit patterns: 0010x01x --exactly 2 bits worth chosen to make things somewhat readable */
#define MLTK__IS_TEXT( tok) ( (tok & 0x00F6) == 0x0022)


#define	MLTKSTR     '"'         /* 34 0x22 00100010 */
#define	MLTKSYM     '#'         /* 35 0x23 # 00100011 */

#define MLTKOPTSYM  'O'       /* 79 00101010 */
#define MLTKOPTSTR  'Q'       /* 81 01010001 */

#define	MLTKREAL    '*'         /* 42 0x2A 00101010 */
#define	MLTKINT     '+'         /* 43 0x2B 00101011 */



/* The following defines are for internal use only */
#define	MLTKPCTEND  ']'     /* at end of top level expression */
#define	MLTKAPCTEND '\n'    /* at end of top level expression */
#define	MLTKEND     '\n'
#define	MLTKAEND    '\r'
#define	MLTKSEND    ','

#define	MLTKCONT    '\\'
#define	MLTKELEN    ' '

#define	MLTKNULL    '.'
#define	MLTKOLDSYM  'Y'     /* 89 0x59 01011001 */
#define	MLTKOLDSTR  'S'     /* 83 0x53 01010011 */


typedef unsigned long decoder_mask;
#define	MLTKPACKED	'P'     /* 80 0x50 01010000 */
#define	MLTKARRAY	'A'     /* 65 0x41 01000001 */
#define	MLTKDIM		'D'     /* 68 0x44 01000100 */

#define MLLENGTH_DECODER        ((decoder_mask) 1<<16)
#define MLTKPACKED_DECODER      ((decoder_mask) 1<<17)
#define MLTKARRAY_DECODER	    ((decoder_mask) 1<<18)
#define MLTKMODERNCHARS_DECODER ((decoder_mask) 1<<19)
#if 0
#define MLTKNULLSEQUENCE_DECODER ((decoder_mask) 1<<20)
#else
#define MLTKNULLSEQUENCE_DECODER ((decoder_mask) 0)
#endif
#define MLTKALL_DECODERS (MLLENGTH_DECODER | MLTKPACKED_DECODER | MLTKARRAY_DECODER | MLTKMODERNCHARS_DECODER | MLTKNULLSEQUENCE_DECODER)

#define MLTK_FIRSTUSER '\x30' /* user token */
#define MLTK_LASTUSER  '\x3F'



#endif /* _MLTK_H */


/*************** Native C types interface ***************/


#ifndef MLGETNUMBERS_HPP
#define MLGETNUMBERS_HPP






#endif /* MLGETNUMBERS_HPP */


/* explicitly not protected by MLGETNUMBERS_HPP in case MLDECL is redefined for multiple inclusion */


ML_EXTERN_C

MLDECL( int,   MLGetBinaryNumber,  ( MLINK mlp, void *np, long type));

/*
As of MLINTERFACE 3 the following functions have been deprecated.  Use the suggested functions in their
place:

MLGetShortInteger   - MLGetInteger16
MLGetInteger        - MLGetInteger32
MLGetLongInteger    - MLGetInteger64 for 64-bit integers or MLGetInteger32 for 32-bit integers
*/
MLDECL( int,   MLGetShortInteger,  ( MLINK mlp, short *hp));
MLDECL( int,   MLGetInteger,       ( MLINK mlp, int *ip));
MLDECL( int,   MLGetLongInteger,   ( MLINK mlp, long *lp));


MLDECL( int,   MLGetInteger16,  ( MLINK mlp, short *hp));
MLDECL( int,   MLGetInteger32,  ( MLINK mlp, int *ip));
MLDECL( int,   MLGetInteger64,  ( MLINK mlp, mlint64 *wp));


MLDECL(int, MLGetInteger8, (MLINK mlp, unsigned char *cp));


/*
As of MLINTERFACE 3 the following functions have been deprecated.  Use the suggested functions in their
place:

MLGetFloat      - MLGetReal32
MLGetDouble     - MLGetReal64
MLGetReal       - MLGetReal64
MLGetLongDouble - MLGetReal128
*/
MLDECL( int,   MLGetFloat,         ( MLINK mlp, float *fp));
MLDECL( int,   MLGetDouble,        ( MLINK mlp, double *dp));
MLDECL( int,   MLGetReal,          ( MLINK mlp, double *dp));
MLDECL( int,   MLGetLongDouble,    ( MLINK mlp, mlextended_double *xp));


MLDECL( int,   MLGetReal32,         ( MLINK mlp, float *fp));
MLDECL( int,   MLGetReal64,        ( MLINK mlp, double *dp));
MLDECL( int,   MLGetReal128,          ( MLINK mlp, mlextended_double *dp));

ML_END_EXTERN_C




#ifndef MLGETSTRINGS_HPP
#define MLGETSTRINGS_HPP






#endif /* MLGETSTRINGS_HPP */

/* explicitly not protected by MLGETSTRINGS_HPP in case MLDECL is redefined for multiple inclusion */

ML_EXTERN_C

/*
As of MLINTERFACE 3 MLGet16BitCharacters has been deprecated.  Use the suggested function in its
place:

MLGet16BitCharacters   - MLGetUCS2Characters
*/

MLDECL( int,   MLGet8BitCharacters,   ( MLINK mlp, long *chars_left, unsigned char *buf, long cardof_buf, long *got, long missing));
MLDECL( int,   MLGet7BitCharacters,   ( MLINK mlp, long *chars_left, char *buf, long cardof_buf, long *got));

MLDECL( int,   MLGetUCS2Characters,   ( MLINK mlp, int *chars_left, unsigned short *buf, int cardof_buf, int *got));
MLDECL( int,   MLGetUTF8Characters,   ( MLINK mlp, int *chars_left, unsigned char *buf, int cardof_buf, int *got));
MLDECL( int,   MLGetUTF16Characters,  ( MLINK mlp, int *chars_left, unsigned short *buf, int cardof_buf, int *got));
MLDECL( int,   MLGetUTF32Characters,  ( MLINK mlp, int *chars_left, unsigned int *buf, int cardof_buf, int *got));

/*
As of MLINTERFACE 3 MLGetUnicodeString has been deprecated.  Use the suggested function in its
place:

MLGetUnicodeString - MLGetUCS2String
*/

MLDECL( int,   MLGetByteString,       ( MLINK mlp, const unsigned char **sp, int *lenp, long missing));
MLDECL( int,   MLGetString,           ( MLINK mlp, const char **sp));

MLDECL( int,   MLGetUCS2String,       ( MLINK mlp, const unsigned short **sp, int *lenp));
MLDECL( int,   MLGetUTF8String,       ( MLINK mlp, const unsigned char **sp, int *bytes, int *chars));
MLDECL( int,   MLGetUTF16String,      ( MLINK mlp, const unsigned short **sp, int *ncodes, int *chars));
MLDECL( int,   MLGetUTF32String,      ( MLINK mlp, const unsigned int **sp, int *len));

#ifndef MLINTERFACE
#error
#endif

/*
As of MLINTERFACE 3 MLGetUnicodeString0 has been deprecated.  Use the suggested function in its
place:

MLGetUnicodeString0 - MLGetUCS2String0
*/


MLDECL( int,   MLGetNumberAsByteString,      ( MLINK mlp, const unsigned char **sp, long *lenp, long missing));
MLDECL( int,   MLGetNumberAsString,          ( MLINK mlp, const char **sp));

MLDECL( int,   MLGetNumberAsUCS2String,      ( MLINK mlp, const unsigned short **sp, int *lenp));
MLDECL( int,   MLGetNumberAsUTF8String,      ( MLINK mlp, const unsigned char **sp, int *bytes, int *chars));
MLDECL( int,   MLGetNumberAsUTF16String,     ( MLINK mlp, const unsigned short **sp, int *ncodes, int *chars));
MLDECL( int,   MLGetNumberAsUTF32String,     ( MLINK mlp, const unsigned int **sp, int *lenp));


/*
As of MLINTERFACE 3 the following functions have been deprecated.  Use the suggested functions in their
place:

MLDisownUnicodeString - MLReleaseUCS2String
MLDisownByteString    - MLReleaseByteString
MLDisownString        - MLReleaseString
*/

MLDECL( void,           MLReleaseUCS2String,   ( MLINK mlp, const unsigned short *s,   int len));
MLDECL( void,           MLReleaseUTF8String,   ( MLINK mlp, const unsigned char *s, int len));
MLDECL( void,           MLReleaseUTF16String,  ( MLINK mlp, const unsigned short *s, int len));
MLDECL( void,           MLReleaseUTF32String,  ( MLINK mlp, const unsigned int *s, int len));
MLDECL( void,           MLReleaseByteString,   ( MLINK mlp, const unsigned char * s,   int len));
MLDECL( void,           MLReleaseString,       ( MLINK mlp, const char *s));


MLDECL( int,    MLTestString,      ( MLINK mlp, const char *name));
MLDECL( int,    MLTestUCS2String,  ( MLINK mlp, const unsigned short *name, int length));
MLDECL( int,    MLTestUTF8String,  ( MLINK mlp, const unsigned char *name, int length));
MLDECL( int,    MLTestUTF16String, ( MLINK mlp, const unsigned short *name, int length));
MLDECL( int,    MLTestUTF32String, ( MLINK mlp, const unsigned int *name, int length));

ML_END_EXTERN_C




#ifndef MLGETSYMBOLS_HPP
#define MLGETSYMBOLS_HPP






#endif /* MLGETSYMBOLS_HPP */

/* explicitly not protected by MLGETSYMBOLS_HPP in case MLDECL is redefined for multiple inclusion */

ML_EXTERN_C

/*
As of MLINTERFACE 3 MLGetUnicodeSymbol has been deprecated.  Use the suggested function in its
place:

MLGetUnicodeSymbol - MLGetUCS2Symbol
*/

MLDECL( int,   MLGetByteSymbol,       ( MLINK mlp, const unsigned char ** sp, int *lenp, long missing));
MLDECL( int,   MLGetSymbol,           ( MLINK mlp, const char **          sp));

MLDECL( int,   MLGetUCS2Symbol,       ( MLINK mlp, const unsigned short **sp, int *lenp));
MLDECL( int,   MLGetUTF8Symbol,       ( MLINK mlp, const unsigned char **sp, int *bytes, int *chars));
MLDECL( int,   MLGetUTF16Symbol,      ( MLINK mlp, const unsigned short **sp, int *ncodes, int *chars));
MLDECL( int,   MLGetUTF32Symbol,      ( MLINK mlp, const unsigned int **sp, int *lenp));

/*
As of MLINTERFACE 3 the following functions have been deprecated.  Use the suggested functions in their
place:
MLDisownUnicodeSymbol - MLReleaseUCS2Symbol
MLDisownByteSymbol    - MLReleaseByteSymbol
MLDisownSymbol        - MLReleaseSymbol
*/

MLDECL( void,           MLReleaseUCS2Symbol,   ( MLINK mlp, const unsigned short *s,   int len));
MLDECL( void,           MLReleaseUTF8Symbol,   ( MLINK mlp, const unsigned char *s, int len));
MLDECL( void,           MLReleaseUTF16Symbol,  ( MLINK mlp, const unsigned short *s, int len));
MLDECL( void,           MLReleaseUTF32Symbol,  ( MLINK mlp, const unsigned int *s, int len));
MLDECL( void,           MLReleaseByteSymbol,   ( MLINK mlp, const unsigned char * s,   int len));
MLDECL( void,           MLReleaseSymbol,       ( MLINK mlp, const char *s));

MLDECL( int,            MLTestSymbol,          ( MLINK mlp, const char *name));
MLDECL( int,            MLTestUCS2Symbol,      ( MLINK mlp, const unsigned short *name, int length));
MLDECL( int,            MLTestUTF8Symbol,      ( MLINK mlp, const unsigned char *name, int length));
MLDECL( int,            MLTestUTF16Symbol,     ( MLINK mlp, const unsigned short *name, int length));
MLDECL( int,            MLTestUTF32Symbol,     ( MLINK mlp, const unsigned int *name, int length));

MLDECL( int,            MLGetFunction,         ( MLINK mlp, const char **sp, int *countp));

MLDECL( int,            MLGetUCS2Function,     ( MLINK mlp, const unsigned short **sp, int *length, int *countp));
MLDECL( int,            MLGetUTF8Function,     ( MLINK mlp, const unsigned char **sp, int *length, int *countp));
MLDECL( int,            MLGetUTF16Function,    ( MLINK mlp, const unsigned short **sp, int *length, int *countp));
MLDECL( int,            MLGetUTF32Function,    ( MLINK mlp, const unsigned int **sp, int *length, int *countp));

/*
As of MLINTERFACE 3 the following functions have been deprecated.  Use the suggested functions in their
place:

MLCheckFunction             - MLTestHead
MLCheckFunctionWithArgCount - MLTestHead
*/

MLDECL( int,   MLCheckFunction, ( MLINK mlp, const char *s, long *countp));
MLDECL( int,   MLCheckFunctionWithArgCount, ( MLINK mlp, const char *s, long *countp));

MLDECL( int,   MLTestHead,      ( MLINK mlp, const char *s, int *countp));


/*
For a limited time convenience define the following:
*/

MLDECL( int,      MLTestHeadWithArgCount, (MLINK mlp, const char *s, int *countp));
MLDECL( int,  MLTestUCS2HeadWithArgCount, (MLINK mlp, const unsigned short *s, int length, int *countp));
MLDECL( int, MLTestUTF16HeadWithArgCount, (MLINK mlp, const unsigned short *s, int length, int *countp));
MLDECL( int, MLTestUTF32HeadWithArgCount, (MLINK mlp, const unsigned int *s, int length, int *countp));
MLDECL( int,  MLTestUTF8HeadWithArgCount, (MLINK mlp, const unsigned char *s, int length, int *countp));

MLDECL( int,   MLTestUCS2Head,  ( MLINK mlp, const unsigned short *s, int length, int *countp));
MLDECL( int,   MLTestUTF8Head,  ( MLINK mlp, const unsigned char *s, int length, int *countp));
MLDECL( int,   MLTestUTF16Head, ( MLINK mlp, const unsigned short *s, int length, int *countp));
MLDECL( int,   MLTestUTF32Head, ( MLINK mlp, const unsigned int *s, int length, int *countp));

ML_END_EXTERN_C






#ifndef MLPUTNUMBERS_HPP
#define MLPUTNUMBERS_HPP






#endif /* MLPUTNUMBERS_HPP */


/* explicitly not protected by MLPUTNUMBERS_HPP in case MLDECL is redefined for multiple inclusion */

ML_EXTERN_C

MLDECL( int,   MLPutBinaryNumber,  ( MLINK mlp, void *np, long type));

/*
As of MLINTERFACE 3 the following functions have been deprecated.  Use the following new functions as their replacement:

MLPutShortInteger  - MLPutInteger16
MLPutInteger       - MLPutInteger32
MLPutLongInteger   - MLPutInteger64 for 64-bit integers or MLPutInteger32 for 32-bit integers.
*/
MLDECL( int,   MLPutShortInteger,  ( MLINK mlp, int h));
MLDECL( int,   MLPutInteger,       ( MLINK mlp, int i));
MLDECL( int,   MLPutLongInteger,   ( MLINK mlp, long l));

MLDECL( int,   MLPutInteger16,     ( MLINK mlp, int h));
MLDECL( int,   MLPutInteger32,     ( MLINK mlp, int i));
MLDECL( int,   MLPutInteger64,     ( MLINK mlp, mlint64 w));

MLDECL( int, MLPutInteger8,   (MLINK mlp, unsigned char i));

/*
As of MLINTERFACE 3 the following functions have been deprecated.  Use the following new functions as their replacement:

MLPutFloat      - MLPutReal32
MLPutDouble     - MLPutReal64
MLPutReal       - MLPutReal64
MLPutLongDouble - MLPutReal128
*/
MLDECL( int,   MLPutFloat,         ( MLINK mlp, double f));
MLDECL( int,   MLPutDouble,        ( MLINK mlp, double d));
MLDECL( int,   MLPutReal,          ( MLINK mlp, double d));
MLDECL( int,   MLPutLongDouble,    ( MLINK mlp, mlextended_double x));

MLDECL( int,   MLPutReal32,         ( MLINK mlp, double f));
MLDECL( int,   MLPutReal64,         ( MLINK mlp, double d));
MLDECL( int,   MLPutReal128,        ( MLINK mlp, mlextended_double x));

ML_END_EXTERN_C




#ifndef MLPUTSTRINGS_HPP
#define MLPUTSTRINGS_HPP







#endif /* MLPUTSTRINGS_HPP */

/* explicitly not protected by MLPUTSTRINGS_HPP in case MLDECL is redefined for multiple inclusion */

ML_EXTERN_C

/*
As of MLINTERFACE 3 MLPut16BitCharacters has been deprecated.  Use the suggested function in its
place:

MLPut16BitCharacters   - MLPutUCS2Characters
*/

MLDECL( int,   MLPut8BitCharacters,  ( MLINK mlp, long chars_left, const unsigned char *bytes, long nbytes));
MLDECL( int,   MLPut7BitCount,       ( MLINK mlp, long count, long size));
MLDECL( int,   MLPut7BitCharacters,  ( MLINK mlp, long chars_left, const char *bytes, long nbytes, long nchars_now));

MLDECL( int,   MLPutUCS2Characters,  ( MLINK mlp, int chars_left, const unsigned short *codes, int ncodes));
MLDECL( int,   MLPutUTF8Characters,  ( MLINK mlp, int chars_left, const unsigned char *codes, int ncodes));
MLDECL( int,   MLPutUTF16Characters, ( MLINK mlp, int chars_left, const unsigned short *codes, int ncodes));
MLDECL( int,   MLPutUTF32Characters, ( MLINK mlp, int chars_left, const unsigned int *codes, int ncodes));

/*
As of MLINTERFACE 3 MLPutUnicodeString has been deprecated.  Use the suggested function in its
place:

MLPutUnicodeString - MLPutUCS2String
*/

MLDECL( int,   MLPutByteString,    ( MLINK mlp, const unsigned char *s, long len));
MLDECL( int,   MLPutString,        ( MLINK mlp, const char *s));

MLDECL( int,   MLPutUCS2String,    ( MLINK mlp, const unsigned short *s, int len));
MLDECL( int,   MLPutUTF8String,    ( MLINK mlp, const unsigned char *s, int len));
MLDECL( int,   MLPutUTF16String,   ( MLINK mlp, const unsigned short *s, int len));
MLDECL( int,   MLPutUTF32String,   ( MLINK mlp, const unsigned int *s, int len));

#ifndef MLINTERFACE
#error
#endif

/*
As of MLINTERFACE 3 MLPutRealUnicodeString0 has been deprecated.  Use the suggested function in its
place:

MLPutRealUnicodeString0 - MLPutRealUCS2String0
*/

MLDECL( int,   MLPutRealNumberAsString,        ( MLINK mlp, const char *s));
MLDECL( int,   MLPutRealNumberAsByteString,    ( MLINK mlp, const unsigned char *s));
MLDECL( int,   MLPutRealNumberAsUCS2String,    ( MLINK mlp, const unsigned short *s));
MLDECL( int,   MLPutRealNumberAsUTF8String,    ( MLINK mlp, const unsigned char *s, int nbytes));
MLDECL( int,   MLPutRealNumberAsUTF16String,   ( MLINK mlp, const unsigned short *s, int ncodes));
MLDECL( int,   MLPutRealNumberAsUTF32String,   ( MLINK mlp, const unsigned int *s, int nchars));


MLDECL( int,   MLPutSize,          ( MLINK mlp, int size));
MLDECL( int,   MLPutData,          ( MLINK mlp, const char *buff, int len));

ML_END_EXTERN_C




#ifndef MLPUTSYMBOLS_HPP
#define MLPUTSYMBOLS_HPP







#endif /* MLPUTSYMBOLS_HPP */

/* explicitly not protected by MLPUTSYMBOLS_HPP in case MLDECL is redefined for multiple inclusion */

ML_EXTERN_C

/*
As of MLINTERFACE 3 MLPutUnicodeSymbol has been deprecated.  Use the suggested function in its
place:

MLPutUnicodeSymbol - MLPutUCS2Symbol
*/


MLDECL( int,   MLPutByteSymbol,    ( MLINK mlp, const unsigned char *s, long len));
MLDECL( int,   MLPutSymbol,        ( MLINK mlp, const char *s));

MLDECL( int,   MLPutUCS2Symbol,    ( MLINK mlp, const unsigned short *s, int len));

MLDECL( int,   MLPutUTF8Symbol,    ( MLINK mlp, const unsigned char *s, int len));
MLDECL( int,   MLPutUTF16Symbol,   ( MLINK mlp, const unsigned short *s, int len));
MLDECL( int,   MLPutUTF32Symbol,   ( MLINK mlp, const unsigned int *s, int len));


MLDECL( int,   MLPutFunction,      ( MLINK mlp, const char *s, int argc));

MLDECL( int,   MLPutUCS2Function,  ( MLINK mlp, const unsigned short *s, int length, int argn));
MLDECL( int,   MLPutUTF8Function,  ( MLINK mlp, const unsigned char *s, int length, int argn));
MLDECL( int,   MLPutUTF16Function, ( MLINK mlp, const unsigned short *s, int length, int argn));
MLDECL( int,   MLPutUTF32Function, ( MLINK mlp, const unsigned int *s, int length, int argn));


ML_END_EXTERN_C




#ifndef _MLSTRING_H
#define _MLSTRING_H








#define MAX_BYTES_PER_OLD_CHARACTER 3
#define MAX_BYTES_PER_NEW_CHARACTER 10

#define ML_MAX_BYTES_PER_CHARACTER MAX_BYTES_PER_NEW_CHARACTER

/* for source code compatibility with earlier versions of MathLink */

typedef struct {
	const char *str;
	const char *end;
} MLStringPosition;

typedef MLStringPosition * MLStringPositionPointer;

#define MLStringFirstPos(s,pos) MLStringFirstPosFun( s, &(pos))

#define MLforString( s, pos) \
	for( MLStringFirstPos(s,pos); MLStringCharacter( (pos).str, (pos).end) >= 0; MLNextCharacter(&(pos).str, (pos).end))

#define MLStringChar( pos) MLStringCharacter( (pos).str, (pos).end)

#define MLPutCharToString MLConvertCharacter


/* for internal use only */

typedef struct {
	unsigned char *cc;
	int  mode;
	int  more;
	unsigned char *head;
} MLOldStringPosition;

typedef MLOldStringPosition * MLOldStringPositionPointer;


#define MLOldforString( s, pos) \
  for ( MLOldStringFirstPos( s, pos); (pos).more; MLOldStringNextPos( pos))

#define MLOldStringChar(pos) \
  ( ((pos).mode <= 1) ? (uint_ct)(*(ucharp_ct)((pos).cc)) : MLOldStringCharFun( &pos) )


#define MLOldStringFirstPos(s,pos) MLOldStringFirstPosFun( s, &(pos))

#define MLOldStringNextPos(pos)  ( \
	((pos).mode == 0) \
		? ((*(*(pos).cc ? ++(pos).cc : (pos).cc) ? 0 : ((pos).more = 0)), (pos).cc) \
		: MLOldStringNextPosFun( &pos) )

#endif /* _MLSTRING_H */




/* explicitly not protected by _MLXDATA_H in case MLDECL is redefined for multiple inclusion */

ML_EXTERN_C
/* assumes *startp aligned on char boundary, if n == -1, returns ~(char_count) */

MLDECL( long, MLCharacterOffset,           ( const char **startp, const char *end, long n));
MLDECL( long, MLStringCharacter,           ( const char * start,  const char *end));
MLDECL( long, MLNextCharacter,             ( const char **startp, const char *end));


#ifndef MLINTERFACE
#error
#endif

MLDECL( long, MLNextCharacterFromStringWithLength, (const char *str, long *indexp, long len));

MLDECL( long, MLConvertNewLine,            ( char **sp));
MLDECL( long, MLConvertCharacter,          ( unsigned long ch, char **sp));
MLDECL( long, MLConvertByteString,         ( unsigned char *codes, long len, char **strp, char *str_end));
MLDECL( long, MLConvertByteStringNL,       ( unsigned char *codes, long len, char **strp, char *str_end, unsigned long nl));
MLDECL( long, MLConvertDoubleByteString,   ( unsigned char *codes, long len, char **strp, char *str_end));
MLDECL( long, MLConvertDoubleByteStringNL, ( unsigned char *codes, long len, char **strp, char *str_end, unsigned long nl));

MLDECL( long, MLConvertUCS2String,         ( unsigned short *codes, long len, char **strp, char *str_end));
MLDECL( long, MLConvertUCS2StringNL,       ( unsigned short *codes, long len, char **strp, char *str_end, unsigned long nl));
MLDECL( long, MLConvertUTF8String,         ( unsigned char *codes, long len, char **strp, char *str_end));
MLDECL( long, MLConvertUTF8StringNL,       ( unsigned char *codes, long len, char **strp, char *str_end, unsigned long nl));
MLDECL( long, MLConvertUTF16String,        ( unsigned short *codes, long len, char **strp, char *str_end));
MLDECL( long, MLConvertUTF16StringNL,      ( unsigned short *codes, long len, char **strp, char *str_end, unsigned long nl));
MLDECL( long, MLConvertUTF32String,        ( unsigned int *codes, long len, char **strp, char *str_end));
MLDECL( long, MLConvertUTF32StringNL,      ( unsigned int *codes, long len, char **strp, char *str_end, unsigned long nl));










/* for source code compatibility with earlier versions of MathLink */
MLDECL( const char *,     MLStringFirstPosFun,  ( const char *s, MLStringPositionPointer p));

/* for internal use only */
MLDECL( int,                MLOldPutCharToString,      ( unsigned int ch, char **sp));
MLDECL( unsigned char *,    MLOldStringNextPosFun,     ( MLOldStringPositionPointer p));
MLDECL( unsigned char *,    MLOldStringFirstPosFun,    ( char *s, MLOldStringPositionPointer p));
MLDECL( unsigned int,       MLOldStringCharFun,        ( MLOldStringPositionPointer p));
MLDECL( long,               MLOldConvertByteString,    ( unsigned char *codes, long len, char **strp, char *str_end));

MLDECL( long,               MLOldConvertUCS2String,    ( unsigned short *codes, long len, char **strp, char *str_end));

/* Internal functions */
MLDECL( long, MLCharOffset,           ( const char **startp, const char *end, long n, int more));
MLDECL( long, MLNextChar,             ( const char **startp, const char *end, int more, int useSurrogates, int *wasSurrogatePair));


ML_END_EXTERN_C



#ifndef _MLCAPUT_H
#define _MLCAPUT_H





#ifndef MLINTERFACE
#error
#endif

#ifndef __array_meterp__
#define __array_meterp__
typedef struct array_meter * array_meterp;
typedef array_meterp * array_meterpp;
#endif


#endif /* _MLCAPUT_H */


/* explicitly not protected by _MLCAPUT_H in case MLDECL is redefined for multiple inclusion */

ML_EXTERN_C
MLDECL( int,   MLPutArray,                  ( MLINK mlp, array_meterp meterp));

MLDECL( int,   MLPutBinaryNumberArrayData,  ( MLINK mlp, array_meterp meterp, const void *         datap, long count, long type));
MLDECL( int,   MLPutByteArrayData,          ( MLINK mlp, array_meterp meterp, const unsigned char *datap, long count));
MLDECL( int,   MLPutShortIntegerArrayData,  ( MLINK mlp, array_meterp meterp, const short *        datap, long count));
MLDECL( int,   MLPutIntegerArrayData,       ( MLINK mlp, array_meterp meterp, const int *          datap, long count));
MLDECL( int,   MLPutLongIntegerArrayData,   ( MLINK mlp, array_meterp meterp, const long *         datap, long count));

MLDECL( int,   MLPutInteger8ArrayData,      ( MLINK mlp, array_meterp meterp, const unsigned char * datap, int count));
MLDECL( int,   MLPutInteger16ArrayData,     ( MLINK mlp, array_meterp meterp, const short *        datap, int count));
MLDECL( int,   MLPutInteger32ArrayData,     ( MLINK mlp, array_meterp meterp, const int *          datap, int count));
MLDECL( int,   MLPutInteger64ArrayData,     ( MLINK mlp, array_meterp meterp, const mlint64 *      datap, int count));




MLDECL( int,   MLPutFloatArrayData,         ( MLINK mlp, array_meterp meterp, const float * datap, long count));
MLDECL( int,   MLPutDoubleArrayData,        ( MLINK mlp, array_meterp meterp, const double *datap, long count));
MLDECL( int,   MLPutLongDoubleArrayData,    ( MLINK mlp, array_meterp meterp, const mlextended_double *datap, long count));

MLDECL( int,   MLPutReal32ArrayData,        ( MLINK mlp, array_meterp meterp, const float * datap, int count));
MLDECL( int,   MLPutReal64ArrayData,        ( MLINK mlp, array_meterp meterp, const double *datap, int count));
MLDECL( int,   MLPutReal128ArrayData,       ( MLINK mlp, array_meterp meterp, const mlextended_double *datap, int count));

#ifndef ML_USES_NEW_PUTBYTEARRAY_API
#define ML_USES_NEW_PUTBYTEARRAY_API 1
#endif


MLDECL( int,   MLPutBinaryNumberArray,  ( MLINK mlp, const void *         data, const long *dimp, const char **heads, long depth, long type));
MLDECL( int,   MLPutByteArray,          ( MLINK mlp, const unsigned char *data, const int *dims, const char **heads, int depth));
MLDECL( int,   MLPutShortIntegerArray,  ( MLINK mlp, const short *        data, const long *dims, const char **heads, long depth));
MLDECL( int,   MLPutIntegerArray,       ( MLINK mlp, const int *          data, const long *dims, const char **heads, long depth));
MLDECL( int,   MLPutLongIntegerArray,   ( MLINK mlp, const long *         data, const long *dims, const char **heads, long depth));

MLDECL( int,   MLPutInteger8Array,      ( MLINK mlp, const unsigned char *data, const int *dims, const char **heads, int depth));
MLDECL( int,   MLPutInteger16Array,     ( MLINK mlp, const short *        data, const int *dims, const char **heads, int depth));
MLDECL( int,   MLPutInteger32Array,     ( MLINK mlp, const int *          data, const int *dims, const char **heads, int depth));
MLDECL( int,   MLPutInteger64Array,     ( MLINK mlp, const mlint64 *      data, const int *dims, const char **heads, int depth));


MLDECL( int,   MLPutFloatArray,         ( MLINK mlp, const float * data, const long *dims, const char **heads, long depth));
MLDECL( int,   MLPutDoubleArray,        ( MLINK mlp, const double *data, const long *dims, const char **heads, long depth));
MLDECL( int,   MLPutRealArray,          ( MLINK mlp, const double *data, const long *dims, const char **heads, long depth));
MLDECL( int,   MLPutLongDoubleArray,    ( MLINK mlp, const mlextended_double *data, const long *dims, const char **heads, long depth));

MLDECL( int,   MLPutReal32Array,        ( MLINK mlp, const float * data, const int *dims, const char **heads, int depth));
MLDECL( int,   MLPutReal64Array,        ( MLINK mlp, const double *data, const int *dims, const char **heads, int depth));
MLDECL( int,   MLPutReal128Array,       ( MLINK mlp, const mlextended_double *data, const int *dims, const char **heads, int depth));


MLDECL( int,   MLPutBinaryNumberList, ( MLINK mlp, const void *  data, long count, long type));
MLDECL( int,   MLPutIntegerList,      ( MLINK mlp, const int *   data, long count));
MLDECL( int,   MLPutRealList,         ( MLINK mlp, const double *data, long count));

MLDECL( int,   MLPutInteger8List,     ( MLINK mlp, const unsigned char *data, int count));
MLDECL( int,   MLPutInteger16List,    ( MLINK mlp, const short *   data, int count));
MLDECL( int,   MLPutInteger32List,    ( MLINK mlp, const int *     data, int count));
MLDECL( int,   MLPutInteger64List,    ( MLINK mlp, const mlint64 * data, int count));

MLDECL( int,   MLPutReal32List,       ( MLINK mlp, const float * data, int count));
MLDECL( int,   MLPutReal64List,       ( MLINK mlp, const double *data, int count));
MLDECL( int,   MLPutReal128List,      ( MLINK mlp, const mlextended_double *data, int count));


MLDECL( int, MLPutArrayType,             ( MLINK mlp, MLINK heads, long depth, array_meterpp meterpp));
MLDECL( int, MLReleasePutArrayState,     ( MLINK mlp, MLINK heads, array_meterp meterp));

MLDECL( int, MLPutArrayLeaves,           ( MLINK mlp, MLINK heads, array_meterp meterp, MLINK leaves, long count));
MLDECL( int, MLPutBinaryNumberArrayDataWithHeads, ( MLINK mlp, MLINK heads, array_meterp meterp, const void *datap, long count, long type));


ML_END_EXTERN_C



#ifndef _MLCAGET_H
#define _MLCAGET_H





#ifndef MLINTERFACE
#error
#endif

#ifndef __array_meterp__
#define __array_meterp__
typedef struct array_meter * array_meterp;
typedef array_meterp * array_meterpp;
#endif

#endif /* _MLCAGET_H */



/* explicitly not protected by _MLCAGET_H in case MLDECL is redefined for multiple inclusion */

ML_EXTERN_C

MLDECL( int,   MLGetArrayDimensions,       ( MLINK mlp, array_meterp meterp));
MLDECL( int,   MLGetArrayType,             ( MLINK mlp, array_meterp meterp));

MLDECL( int,  MLGetBinaryNumberList, ( MLINK mlp, void **datap, long *countp, long type));

/*
As of MLINTERFACE 3 the following functions have been deprecated.  Use the suggested functions in their
place:

MLGetIntegerList  - MLGetInteger32List
MLGetRealList     - MLGetReal64List
*/
MLDECL( int,  MLGetIntegerList,      ( MLINK mlp, int **datap, long *countp));
MLDECL( int,  MLGetRealList,         ( MLINK mlp, double **datap, long *countp));

MLDECL( int,  MLGetInteger16List,      ( MLINK mlp, short **   datap, int *countp));
MLDECL( int,  MLGetInteger32List,      ( MLINK mlp, int **     datap, int *countp));
MLDECL( int,  MLGetInteger64List,      ( MLINK mlp, mlint64 ** datap, int *countp));

MLDECL( int,  MLGetReal32List,         ( MLINK mlp, float **                 datap, int *countp));
MLDECL( int,  MLGetReal64List,         ( MLINK mlp, double **                datap, int *countp));
MLDECL( int,  MLGetReal128List,         ( MLINK mlp, mlextended_double **    datap, int *countp));

/*
As of MLINTERFACE 3 the following functions have been deprecated.  Use the suggested functions in their
place:

MLDisownBinaryNumberList  - MLReleaseBinaryNumberList
MLDisownIntegerList       - MLReleaseInteger32List
MLDisownRealList          - MLReleaseReal64List
*/

MLDECL( void, MLReleaseIntegerList,     ( MLINK mlp, int *data, long count));
MLDECL( void, MLReleaseRealList,        ( MLINK mlp, double *data, long count));

MLDECL( void, MLReleaseBinaryNumberList,   ( MLINK mlp, void *data, int count, long type));
MLDECL( void, MLReleaseInteger16List,      ( MLINK mlp, short *data, int count));
MLDECL( void, MLReleaseInteger32List,      ( MLINK mlp, int *data, int count));
MLDECL( void, MLReleaseInteger64List,      ( MLINK mlp, mlint64 *data, int count));

MLDECL( void, MLReleaseReal32List,         ( MLINK mlp, float *data, int count));
MLDECL( void, MLReleaseReal64List,         ( MLINK mlp, double *data, int count));
MLDECL( void, MLReleaseReal128List,        ( MLINK mlp, mlextended_double *data, int count));

MLDECL( int,   MLGetBinaryNumberArrayData,  ( MLINK mlp, array_meterp meterp, void *datap, long count, long type));
MLDECL( int,   MLGetByteArrayData,          ( MLINK mlp, array_meterp meterp, unsigned char * datap, long count));

/*
As of MLINTERFACE 3 the following functions have been deprecated.  Use the suggested functions in their
place:

MLGetShortIntegerArrayData  - MLGetInteger16ArrayData
MLGetIntegerArrayData       - MLGetInteger32ArrayData
MLGetLongIntegerArrayData   - MLGetInteger64ArrayData for 64-bit integers or MLGetInteger32ArrayData for 32-bit integers.
*/

MLDECL( int,   MLGetShortIntegerArrayData,  ( MLINK mlp, array_meterp meterp, short *         datap, long count));
MLDECL( int,   MLGetIntegerArrayData,       ( MLINK mlp, array_meterp meterp, int *           datap, long count));
MLDECL( int,   MLGetLongIntegerArrayData,   ( MLINK mlp, array_meterp meterp, long *          datap, long count));

MLDECL( int,   MLGetInteger16ArrayData,     ( MLINK mlp, array_meterp meterp, short *         datap, int count));
MLDECL( int,   MLGetInteger32ArrayData,     ( MLINK mlp, array_meterp meterp, int *           datap, int count));
MLDECL( int,   MLGetInteger64ArrayData,     ( MLINK mlp, array_meterp meterp, mlint64 *       datap, int count));

/*
As of MLINTERFACE 3 the following functions have been deprecated.  Use the suggested functions in their
place:

MLGetFloatArrayData        - MLGetReal32ArrayData
MLGetDoubleArrayData       - MLGetReal64ArrayData
MLGetLongDoubleArrayData   - MLGetReal128ArrayData
*/

MLDECL( int,   MLGetFloatArrayData,         ( MLINK mlp, array_meterp meterp, float *datap, long count));
MLDECL( int,   MLGetDoubleArrayData,        ( MLINK mlp, array_meterp meterp, double *datap, long count));
MLDECL( int,   MLGetLongDoubleArrayData,   ( MLINK mlp, array_meterp meterp, mlextended_double *datap, long count));

MLDECL( int,   MLGetReal32ArrayData,         ( MLINK mlp, array_meterp meterp, float *datap, int count));
MLDECL( int,   MLGetReal64ArrayData,        ( MLINK mlp, array_meterp meterp, double *datap, int count));
MLDECL( int,   MLGetReal128ArrayData,   ( MLINK mlp, array_meterp meterp, mlextended_double *datap, int count));


MLDECL(int, MLGetInteger8List, (MLINK mlp, unsigned char **datap, int *countp));
MLDECL(int, MLGetInteger8ArrayData, (MLINK mlp, array_meterp meterp, unsigned char *datap, int count));
MLDECL(void, MLReleaseInteger8List, (MLINK mlp, unsigned char *data, int count));

MLDECL( int,   MLGetArrayTypeWithDepthAndLeafType, ( MLINK mlp, MLINK heads, array_meterpp meterpp,
	long *depthp, mlapi__token *leaf_tokp));

MLDECL( int,   MLGetBinaryNumberArrayDataWithHeads, ( MLINK mlp, MLINK heads, array_meterp  meterp, void *datap, long *countp, long type));

MLDECL( void,  MLReleaseGetArrayState,     ( MLINK mlp, MLINK heads, array_meterp  meterp));

MLDECL( int,   MLGetBinaryNumberArrayWithLeafType,   ( MLINK mlp, void **datap, long **dimpp, char ***headsp, long *depthp, long type, mlapi__token *leaf_tokp));



MLDECL( int,   MLGetBinaryNumberArray,    ( MLINK mlp, void **          datap, long **dimpp, char ***headsp, long *depthp, long type));
MLDECL( int,   MLGetByteArray,            ( MLINK mlp, unsigned char ** datap, int **dimsp, char ***headsp, int *depthp));

/*
As of MLINTERFACE 3 the following functions have been deprecated.  Use the suggested functions in their
place:

MLGetShortIntegerArray   - MLGetInteger16Array
MLGetIntegerArray        - MLGetInteger32Array
MLGetLongIntegerArray    - MLGetInteger64Array for 64-bit integers or MLGetInteger32Array from 32-bit integers.
*/

MLDECL( int,   MLGetShortIntegerArray,    ( MLINK mlp, short **         datap, long **dimsp, char ***headsp, long *depthp));
MLDECL( int,   MLGetIntegerArray,         ( MLINK mlp, int **           datap, long **dimsp, char ***headsp, long *depthp));
MLDECL( int,   MLGetLongIntegerArray,     ( MLINK mlp, long **          datap, long **dimsp, char ***headsp, long *depthp));

MLDECL( int,   MLGetInteger16Array,       ( MLINK mlp, short **         datap, int **dimsp, char ***headsp, int *depthp));
MLDECL( int,   MLGetInteger32Array,       ( MLINK mlp, int **           datap, int **dimsp, char ***headsp, int *depthp));
MLDECL( int,   MLGetInteger64Array,       ( MLINK mlp, mlint64 **       datap, int **dimsp, char ***headsp, int *depthp));


MLDECL(int,  MLGetInteger8Array,  (MLINK mlp, unsigned char **datap, int **dimsp, char ***headsp, int *depthp));

/*
As of MLINTERFACE 3 the following functions have been deprecated.  Use the suggested functions in their
place:

MLGetFloatArray      - MLGetReal32Array
MLGetDoubleArray     - MLGetReal64Array
MLGetRealArray       - MLGetReal64Array
MLGetLongDoubleArray - MLGetReal128Array
*/

MLDECL( int,   MLGetFloatArray,           ( MLINK mlp, float ** datap, long **dimsp, char ***headsp, long *depthp));
MLDECL( int,   MLGetDoubleArray,          ( MLINK mlp, double **datap, long **dimsp, char ***headsp, long *depthp));
MLDECL( int,   MLGetRealArray,            ( MLINK mlp, double **datap, long **dimsp, char ***headsp, long *depthp));
MLDECL( int,   MLGetLongDoubleArray,      ( MLINK mlp, mlextended_double **datap, long **dimsp, char ***headsp, long *depthp));

MLDECL( int,   MLGetReal32Array,          ( MLINK mlp, float ** datap, int **dimsp, char ***headsp, int *depthp));
MLDECL( int,   MLGetReal64Array,          ( MLINK mlp, double **datap, int **dimsp, char ***headsp, int *depthp));
MLDECL( int,   MLGetReal128Array,         ( MLINK mlp, mlextended_double **datap, int **dimsp, char ***headsp, int *depthp));


/*
As of MLINTERFACE 3 the following functions have been deprecated.  Use the following new functions as their replacement:

MLDisownBinaryNumberArray - MLReleaseBinaryNumberArray
MLDisownByteArray         - MLReleaseByteArray
MLDisownShortIntegerArray - MLReleaseInteger16Array
MLDisownIntegerArray      - MLReleaseInteger32Array
MLDisownLongIntegerArray  - Use MLReleaseInteger32Array for 32-bit integers or MLReleaseInteger64Array for 64-bit integers.
*/

MLDECL( void,           MLReleaseShortIntegerArray, ( MLINK mlp, short *        data, long *dims, char **heads, long depth));
MLDECL( void,           MLReleaseIntegerArray,      ( MLINK mlp, int *          data, long *dims, char **heads, long depth));
MLDECL( void,           MLReleaseLongIntegerArray,  ( MLINK mlp, long *         data, long *dims, char **heads, long depth));

MLDECL( void,           MLReleaseBinaryNumberArray,  ( MLINK mlp, void *         data, int *dimp, char **heads, int len, long type));
MLDECL( void,           MLReleaseByteArray,          ( MLINK mlp, unsigned char *data, int *dims, char **heads, int depth));
MLDECL( void,           MLReleaseInteger16Array,     ( MLINK mlp, short *        data, int *dims, char **heads, int depth));
MLDECL( void,           MLReleaseInteger32Array,     ( MLINK mlp, int *          data, int *dims, char **heads, int depth));
MLDECL( void,           MLReleaseInteger64Array,     ( MLINK mlp, mlint64 *      data, int *dims, char **heads, int depth));

MLDECL(void,   MLReleaseInteger8Array,   (MLINK mlp, unsigned char *data, int *dimp, char **heads, int depth));

/*
As of MLINTERFACE 3 the following functions have been deprecated.  Use the following new functions as their replacement:

MLDisownFloatArray  - MLReleaseReal32Array
MLDisownDoubleArray - MLReleaseReal64Array
MLDisownRealArray   - MLReleaseReal64Array
*/

MLDECL( void,           MLReleaseFloatArray,         ( MLINK mlp, float * data, long *dims, char **heads, long depth));
MLDECL( void,           MLReleaseDoubleArray,        ( MLINK mlp, double *data, long *dims, char **heads, long depth));
MLDECL( void,           MLReleaseRealArray,          ( MLINK mlp, double *data, long *dims, char **heads, long depth));

MLDECL( void,           MLReleaseReal32Array,         ( MLINK mlp, float * data, int *dims, char **heads, int depth));
MLDECL( void,           MLReleaseReal64Array,          ( MLINK mlp, double *data, int *dims, char **heads, int depth));


/*
As of MLINTERFACE 3 the following functions have been deprecated.  Use the following new functions as their replacement:

MLDisownLongDoubleArray - MLReleaseReal128Array
*/

MLDECL( void,           MLReleaseReal128Array,     ( MLINK mlp, mlextended_double *data, int *dims, char **heads, int depth));

MLDECL( void,           MLReleaseLongDoubleArray,  ( MLINK mlp, mlextended_double *data, long *dims, char **heads, long depth));



ML_END_EXTERN_C


/*************** Unicode Container for mprep templates ***************/


#ifndef MLUNICODECONTAINER_HPP
#define MLUNICODECONTAINER_HPP






ML_EXTERN_C

enum MLUnicodeContainerType
{
	UCS2ContainerType,
	UTF8ContainerType,
	UTF16ContainerType,
	UTF32ContainerType
};

typedef struct _MLUnicodeContainer
{
	union _pointer
	{
		unsigned short *ucs2;
		unsigned char *utf8;
		unsigned short *utf16;
		unsigned int *utf32;
	} pointer;

	int length;
	enum MLUnicodeContainerType type;
} MLUnicodeContainer;

#define MLUCS2String(container)   container->pointer.ucs2
#define MLUTF8String(container)   container->pointer.utf8
#define MLUTF16String(container)  container->pointer.utf16
#define MLUTF32String(container)  container->pointer.utf32
#define MLUnicodeStringLength(container) container->length
#define MLUnicodeStringType(container) container->type

MLDECL(MLUnicodeContainer *, MLNewUnicodeContainer, (const void *string, int length, enum MLUnicodeContainerType type));
MLDECL(void,    MLReleaseUnicodeContainer, (MLUnicodeContainer *string));

ML_END_EXTERN_C


#endif /* MLUNICODECONTAINER_HPP */



/*************** seeking, transfering  and synchronization ***************/

#ifndef _MLMARK_H
#define _MLMARK_H




#endif /* _MLMARK_H */

/* explicitly not protected by _MLMARK_H in case MLDECL is redefined for multiple inclusion */

ML_EXTERN_C
MLDECL( MLINKMark,  MLCreateMark,  ( MLINK mlp));
MLDECL( MLINKMark,  MLSeekToMark,  ( MLINK mlp, MLINKMark mark, int index));
MLDECL( MLINKMark,  MLSeekMark,    ( MLINK mlp, MLINKMark mark, int index));
MLDECL( void,       MLDestroyMark, ( MLINK mlp, MLINKMark mark));
ML_END_EXTERN_C



#ifndef _MLXFER_H
#define _MLXFER_H





#endif /* _MLXFER_H */

/* explicitly not protected by _MLXFER_H in case MLDECL is redefined for multiple inclusion */

ML_EXTERN_C

#ifndef MLINTERFACE
#error
#endif

MLDECL( int, MLTransferExpression, ( MLINK dmlp, MLINK smlp));
MLDECL( int, MLTransferToEndOfLoopbackLink, ( MLINK dmlp, MLINK smlp));


ML_END_EXTERN_C



#ifndef _MLSYNC_H
#define _MLSYNC_H




/* export mls__wait and mls__align(mlsp) */

#endif /* _MLSYNC_H */

/* explicitly not protected by _MLSYNC_H in case MLDECL is redefined for multiple inclusion */

ML_EXTERN_C
/* in response to a reset message */
MLDECL( int, MLForwardReset, ( MLINK mlp, unsigned long marker));


MLDECL( int, MLAlign,        ( MLINK lmlp, MLINK rmlp));
ML_END_EXTERN_C


/*************************************************************/


#ifndef _MLPKT_H
#define _MLPKT_H

/*************** Mathematica packet interface ***************/

			/* MLNextPacket returns one of... */


/* edit here and in mlpktstr.h */

#ifndef _MLPKTNO_H
#define _MLPKTNO_H

#define ILLEGALPKT      0

#define CALLPKT         7
#define EVALUATEPKT    13
#define RETURNPKT       3

#define INPUTNAMEPKT    8
#define ENTERTEXTPKT   14
#define ENTEREXPRPKT   15
#define OUTPUTNAMEPKT   9
#define RETURNTEXTPKT   4
#define RETURNEXPRPKT  16

#define DISPLAYPKT     11
#define DISPLAYENDPKT  12

#define MESSAGEPKT      5
#define TEXTPKT         2

#define INPUTPKT        1
#define INPUTSTRPKT    21
#define MENUPKT         6
#define SYNTAXPKT      10

#define SUSPENDPKT     17
#define RESUMEPKT      18

#define BEGINDLGPKT    19
#define ENDDLGPKT      20

#define FIRSTUSERPKT  128
#define LASTUSERPKT   255


#endif /* _MLPKTNO_H */



#endif /* _MLPKT_H */

/* explicitly not protected by _MLPKT_H in case MLDECL is redefined for multiple inclusion */

ML_EXTERN_C
MLDECL( int,  MLNextPacket, ( MLINK mlp));
ML_END_EXTERN_C



#ifndef _MLALERT_H
#define _MLALERT_H






ML_EXTERN_C
/*************** User interaction--for internal use only ***************/

#if WIN64_MATHLINK
typedef __int64 mldlg_result;
#else
typedef long mldlg_result;
#endif

MLDPROC( mldlg_result, MLAlertProcPtr,             ( MLEnvironment env, const char *message));
MLDPROC( mldlg_result, MLRequestProcPtr,           ( MLEnvironment env, const char *prompt, char *response, long sizeof_response));
MLDPROC( mldlg_result, MLConfirmProcPtr,           ( MLEnvironment env, const char *question, mldlg_result default_answer));
MLDPROC( mldlg_result, MLRequestArgvProcPtr,       ( MLEnvironment env, char **argv, long cardof_argv, char *buf, long sizeof_buf));
MLDPROC( mldlg_result, MLRequestToInteractProcPtr, ( MLEnvironment env, mldlg_result wait_for_permission));
MLDPROC( mldlg_result, MLDialogProcPtr,            ( MLEnvironment env));


typedef MLDialogProcPtr MLDialogUPP;
typedef MLAlertProcPtr MLAlertUPP;
typedef MLRequestProcPtr MLRequestUPP;
typedef MLConfirmProcPtr MLConfirmUPP;
typedef MLRequestArgvProcPtr MLRequestArgvUPP;
typedef MLRequestToInteractProcPtr MLRequestToInteractUPP;
#define NewMLAlertProc(userRoutine) MLAlertCast((userRoutine))
#define NewMLRequestProc(userRoutine) MLRequestCast((userRoutine))
#define NewMLConfirmProc(userRoutine) MLConfirmCast((userRoutine))
#define NewMLRequestArgvProc(userRoutine) MLRequestArgvCast((userRoutine))
#define NewMLRequestToInteractProc(userRoutine) MLRequestToInteractCast((userRoutine))

typedef MLAlertUPP MLAlertFunctionType;
typedef MLRequestUPP MLRequestFunctionType;
typedef MLConfirmUPP MLConfirmFunctionType;
typedef MLRequestArgvUPP MLRequestArgvFunctionType;
typedef MLRequestToInteractUPP MLRequestToInteractFunctionType;
typedef MLDialogUPP MLDialogFunctionType;



/*
	MLDDECL( mldlg_result, alert_user, ( MLEnvironment env, kcharp_ct message));
	MLDDEFN( mldlg_result, alert_user, ( MLEnvironment env, kcharp_ct message))
	{
		fprintf( stderr, "%s\n", message);
	}


	...
	MLDialogFunctionType f = NewMLAlertProc(alert_user);
	MLSetDialogFunction( ep, MLAlertFunction, f);
	...
	or
	...
	MLSetDialogFunction( ep, MLAlertFunction, NewMLAlertProc(alert_user));
	...
*/



enum {	MLAlertFunction = 1, MLRequestFunction, MLConfirmFunction,
	MLRequestArgvFunction, MLRequestToInteractFunction };


#define ML_DEFAULT_DIALOG ( (MLDialogFunctionType) 1)
#define ML_IGNORE_DIALOG ( (MLDialogFunctionType) 0)
#define ML_SUPPRESS_DIALOG ML_IGNORE_DIALOG



#if WINDOWS_MATHLINK

#ifndef _MLWIN_H
#define _MLWIN_H





ML_EXTERN_C
MLDDECL( mldlg_result, MLAlert_win,   ( MLEnvironment ep, const char *alertstr));
MLDDECL( mldlg_result, MLRequest_win, ( MLEnvironment ep, const char *prompt, char *response, long n));
MLDDECL( mldlg_result, MLConfirm_win, ( MLEnvironment ep, const char *okcancelquest, mldlg_result default_answer));
MLDDECL( mldlg_result, MLPermit_win,  ( MLEnvironment ep, mldlg_result wait));
ML_END_EXTERN_C

/* edit here and in mlwin.rc -- in both places because of command-line length limitations in dos */
#define DLG_LINKNAME                101
#define DLG_TEXT                    102
#define RIDOK                       1
#define RIDCANCEL                   104

#endif /* _MLWIN_H */

#define MLALERT         MLAlert_win
#define MLREQUEST       MLRequest_win
#define MLCONFIRM       MLConfirm_win
#define MLPERMIT        MLPermit_win
#define MLREQUESTARGV	default_request_argv
#endif

#if UNIX_MATHLINK
#if DARWIN_MATHLINK  && ! defined(IOS_MATHLINK) && defined (USE_CF_DIALOGS)

#ifndef _MLDARWIN_H
#define _MLDARWIN_H




ML_EXTERN_C

MLDDECL( mldlg_result, MLAlert_darwin,   ( MLEnvironment env, const char *message));
MLDDECL( mldlg_result, MLRequest_darwin, ( MLEnvironment env, const char *prompt, char *response, long sizeof_response));
MLDDECL( mldlg_result, MLConfirm_darwin, ( MLEnvironment env, const char *question, mldlg_result default_answer));
MLDDECL( mldlg_result, MLPermit_darwin,  ( MLEnvironment env, mldlg_result wait_for_permission));
MLDDECL( mldlg_result, MLDontPermit_darwin, ( MLEnvironment ep, mldlg_result wait_for_permission));

ML_END_EXTERN_C

#endif /* _MLDARWIN_H */

#define MLALERT  	MLAlert_darwin
#define MLREQUEST	MLRequest_darwin
#define MLCONFIRM	MLConfirm_darwin
#define MLPERMIT 	MLPermit_darwin
#define MLREQUESTARGV	default_request_argv
#else /* !(DARWIN_MATHLINK && ! defined(IOS_MATHLINK) && defined (USE_CF_DIALOGS)) */

#ifndef _MLUNIX_H
#define _MLUNIX_H




ML_EXTERN_C

MLDDECL( mldlg_result, MLAlert_unix,   ( MLEnvironment env, const char *message));
MLDDECL( mldlg_result, MLRequest_unix, ( MLEnvironment env, const char *prompt, char *response, long sizeof_response));
MLDDECL( mldlg_result, MLConfirm_unix, ( MLEnvironment env, const char *question, mldlg_result default_answer));
MLDDECL( mldlg_result, MLPermit_unix,  ( MLEnvironment env, mldlg_result wait_for_permission));

ML_END_EXTERN_C

#endif /* _MLUNIX_H */

#define MLALERT  	MLAlert_unix
#define MLREQUEST	MLRequest_unix
#define MLCONFIRM	MLConfirm_unix
#define MLPERMIT 	MLPermit_unix
#define MLREQUESTARGV	default_request_argv
#endif /* DARWIN_MATHLINK && ! defined(IOS_MATHLINK) && defined (USE_CF_DIALOGS) */
#endif


MLDDECL( mldlg_result, default_request_argv, ( MLEnvironment ep, char **argv, long len, char *buff, long size));
ML_END_EXTERN_C

#endif /* _MLALERT_H */


/* explicitly not protected by _MLXDATA_H in case MLDECL is redefined for multiple inclusion */
ML_EXTERN_C

MLDECL( mldlg_result,  MLAlert,             ( MLEnvironment env, const char *message));
MLDECL( mldlg_result,  MLRequest,           ( MLEnvironment env, const char *prompt, char *response, long sizeof_response)); /* initialize response with default*/
MLDECL( mldlg_result,  MLConfirm,           ( MLEnvironment env, const char *question, mldlg_result default_answer));
MLDECL( mldlg_result,  MLRequestArgv,       ( MLEnvironment env, char **argv, long cardof_argv, char *buff, long size));

MLDECL( mldlg_result,  MLRequestToInteract, ( MLEnvironment env, mldlg_result wait_for_permission));
MLDECL( int,  MLSetDialogFunction, ( MLEnvironment env, long funcnum, MLDialogFunctionType func));

/* just some type-safe casts */
MLDECL( MLDialogProcPtr, MLAlertCast, ( MLAlertProcPtr f));
MLDECL( MLDialogProcPtr, MLRequestCast, ( MLRequestProcPtr f));
MLDECL( MLDialogProcPtr, MLConfirmCast, ( MLConfirmProcPtr f));
MLDECL( MLDialogProcPtr, MLRequestArgvCast, ( MLRequestArgvProcPtr f));
MLDECL( MLDialogProcPtr, MLRequestToInteractCast, ( MLRequestToInteractProcPtr f));
ML_END_EXTERN_C


/*************************************************************/


#ifndef _MLREADY_H
#define _MLREADY_H

#ifndef _MLTIME_H
#define _MLTIME_H


typedef struct _mltimeval{
	unsigned long tv_sec;
	unsigned long tv_usec;
} mltimeval;


#endif /* _MLTIME_H */



ML_EXTERN_C
MLDECL( int,   MLReady,            ( MLINK mlp));
ML_END_EXTERN_C



#if defined(__cplusplus)
#define MLINFINITEWAIT static_cast<mltimeval>({-1, -1})
#else
#define MLINFINITEWAIT {-1, -1}
#endif /* __cplusplus */

#define MLREADYPARALLELERROR -1
#define MLREADYPARALLELTIMEDOUT -2
#define MLREADYPARALLELINVALIDARGUMENT -3

typedef void *MLREADYPARALLELENV;

ML_EXTERN_C
MLDECL(int,              MLReadyParallel, (MLENV, MLINK *, int, mltimeval));
ML_END_EXTERN_C


ML_EXTERN_C
MLWPROC(int, MLLinkWaitCallBackObject, (MLINK, void *));

MLDECL(int, MLWaitForLinkActivity, (MLINK mlp));
MLDECL(int, MLWaitForLinkActivityWithCallback, (MLINK mlp, MLLinkWaitCallBackObject callback));
ML_END_EXTERN_C


#define MLWAITSUCCESS         1
#define MLWAITERROR           2
#define MLWAITTIMEOUT         3
#define MLWAITCALLBACKABORTED 4

#endif /* _MLREADY_H */


/********************************************************/


#ifndef _MLTM_H
#define _MLTM_H




/*************** Template interface ***************/

/* The following are useful only when using template files as
 * their definitions are produced by mprep.
 */

ML_EXTERN_C
extern MLINK stdlink;
extern MLEnvironment stdenv;

extern MLYieldFunctionObject stdyielder;
extern MLMessageHandlerObject stdhandler;

extern int MLMain(int, char **); /* pass in argc and argv */
extern int MLMainString( char *commandline);

extern int MLMainArgv( char** argv, char** argv_end); /* note not FAR pointers */

extern int MLInstall(MLINK);
extern mlapi_packet MLAnswer(MLINK);
extern int MLDoCallPacket(MLINK);
extern int MLEvaluate( MLINK, char *);
extern int MLEvaluateString( MLINK, char *);
ML_END_EXTERN_C

MLMDECL( void, MLDefaultHandler, ( MLINK, int, int));

MLYDECL( int, MLDefaultYielder, ( MLINK, MLYieldParameters));

ML_EXTERN_C
extern int MLAbort, MLDone;
extern long MLSpecialCharacter;
ML_END_EXTERN_C

#endif /* _MLTM_H */


/*************** Backwards compatibility ***************/

#ifndef P
#define P(s) s
#endif

// Kept for backwards compatibility.
#ifndef FAR
#define FAR
#endif


#endif /* _MATHLINK_H */
