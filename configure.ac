AC_INIT([KAMU], [0.01], [renwick.james.hudspith@gmail.com])
AM_INIT_AUTOMAKE([-Wall subdir-objects])
AC_PROG_CC
AM_PROG_CC_C_O
AC_LANG([C]) 
AC_PROG_RANLIB ## include a check for libtool if not use ranlib?
AM_PROG_AR
AC_C_BIGENDIAN

## SET up the m4
AC_CONFIG_MACRO_DIR([m4])

## prefix!!
AC_DEFINE_UNQUOTED([HAVE_PREFIX],["$prefix"],[Look for where we prefix to])

## tell the code where the PRECOMP library is -> Here & only here!
AC_DEFINE_UNQUOTED([LIB_PATH],["$PWD"],[Where the library and precomps are])
	  
## Tells us where we are prefixed to
AC_MSG_NOTICE([Configuring for location "$prefix"])
## If we prefix away from the default /usr/ then create binary GLU
## will still try and create shared object libraries there
AM_CONDITIONAL([PREF],[ test "$prefix" == "NONE" ])
   if test "$prefix" == "NONE" ; then
       AC_MSG_NOTICE([Not creating binary "URFIT"])
   else
       AC_MSG_NOTICE([Creating binary "URFIT" in $prefix/bin])
       mkdir $prefix  
       mkdir $prefix/bin/
       cp -vp PRECOMP/REG* $prefix/bin
   fi
AC_SUBST(PREF)

## C standard libraries
AC_CHECK_HEADERS([math.h],[LDFLAGS="$LDFLAGS -lm"],
        [AC_MSG_ERROR([We need math.h])])

## Check for sys/time.h
AC_CHECK_HEADERS([sys/time.h],
	[AC_MSG_NOTICE([We are intrinsically timing our operations.])],[])

## Check for the time.h, which we use for the date
AC_CHECK_HEADERS([time.h],
	[AC_MSG_NOTICE([We are providing an accurate date.])],[])

## Check for the time.h, which we use for the date
AC_CHECK_HEADERS([immintrin.h],
	[AC_MSG_NOTICE([We are vectorising the crc32c checksums.])],[])

## I doubt we will not have these as they are all c-standard, check anyway
AC_CHECK_HEADERS([stdbool.h stdint.h strings.h stdio.h stdlib.h],[],[
AC_MSG_ERROR([Sorry, but we need these headers (strings.h , stdlib.h , stdio.h) somewhere in the code])])

## Some compiler checks
## My code is littered with consts
AC_C_CONST([])

## check that we have KQED library
KQED=false 
AC_ARG_WITH([KQED],
	AC_HELP_STRING([--with-KQED=<name>],
		       [Specify KQED location]),
            [kqed=${with_KQED}
	    ## test to see if something is there ...
	    if [ test -a ${kqed}include/KQED.h ]; then
	    	AC_MSG_NOTICE([User set KQED location $kqed])
	  	## check if the gsl static lib is available
	    	if [ test -a ${kqed}lib/libKQED.a ]; then 
		   AC_MSG_NOTICE([Statically linking to the user's KQED lib])
	    	   CFLAGS="$CFLAGS -I${kqed}include"
	    	   LDFLAGS="$LDFLAGS ${kqed}lib/libKQED.a"	
		   AC_DEFINE([HAVE_KQED],,[Do we have the KQED libs])
	    	else
	           AC_MSG_ERROR([KQED.h found, but cannot find the static or shared libraries needed to link against!])
	    	fi
		echo " "
	    else
		AC_MSG_ERROR([Cannot find KQED.h in ${KQED}/include/])
	    fi
		],[ AC_MSG_ERROR([Code requires KQED library.]) ])
AM_CONDITIONAL([DEFKQED],[$KQED])
AC_SUBST(KQED)

## Set up config.h, one of the most important headers of them all
AC_CONFIG_HEADERS([config.h])

## where are my files?
AC_CONFIG_FILES([
Makefile
src/Makefile
])

## tell us what flags we are using
echo CFLAGS used :: $CFLAGS
echo LDFLAGS used :: $LDFLAGS

AC_OUTPUT
