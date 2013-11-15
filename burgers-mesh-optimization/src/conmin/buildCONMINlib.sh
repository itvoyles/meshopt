#!/bin/sh
######################
# Build Conmin shared object/dynamic library Script
# This script will work for linux and Mac
######################
################################################################################
# Script Banner
echo "";
echo "================= ";
echo " BUILDING Conmin Library ";
echo "================= ";
################################################################################
# Base directory
BASEDIR=`pwd`;

RELEASE="RELEASE";
DEBUG="DEBUG";
BUILD=${RELEASE};

# Source directory
SRCDIR=${BASEDIR}/src;

# Kernel type
KERNEL=`uname`;

# Machine dependent directory
DIST="CONMIN-${KERNEL}";
DISTDIR="${BASEDIR}/${DIST}";


# Compilers, compiler flags, include paths, library paths, and libraries
if [ ${KERNEL} = 'Linux' ] || [ ${KERNEL}  = 'Darwin' ]; then	
	# C compiler
	CC="gcc";
	# FORTRAN 77 compiler
	F77="gfortran";
	# Compiler flags
	CFLAGS="-c -m64 -fPIC -O3";
	F77FLAGS="-c -m64 -fPIC -O3";
        LIBNAME="libconmin.dylib";
	if [ ${KERNEL} = 'Linux' ]; then
		F77FLAGS="${F77FLAGS} ";
                LIBNAME="libconmin.so";
	fi
	# Include path
	CPPFLAGS="-DLOWER_UNDERSCORE";
	# Library path
	LDFLAGS="-m64 -fPIC -shared -O3";
	# Libraries
	LIBS=gfortran
# not tested as of 19 Aug 2010
elif [ ${KERNEL} = 'MINGW32_NT-6.1' ]; then
        # C compiler
        CC="gcc";
        # FORTRAN 77 compiler
        F77="gfortran";
        # Compiler flags
        CFLAGS="-c -fPIC";
        F77FLAGS="-c -fPIC ";
        # Include path
        # Library path
        LDFLAGS=-shared
        # Libraries
        LIBS=gfortran
# not tested as of 19 Aug 2010
elif [ ${KERNEL} = 'IRIX64' ]; then 
	# C compiler 
	CC="cc";
	# FORTRAN 77 compiler
	F77="f77";
	# Compiler flags
	CFLAGS="-c -fPIC -64 ";
	F77FLAGS="-c  -64 -fPIC";
	# Include path
	CPPFLAGS="-DLOWER_UNDERSCORE -D PP_SGI";
	# Library path
	LDFLAGS=
	# Libraries
	LIBS=
	# Ranlib
	RANLIB="ar -ts";
# not tested as of 19 Aug 2010
else
	# C compiler
	CC="gcc";
	# FORTRAN 77 compiler
	F77="gfortran";
	# Compiler flags
	CFLAGS="-c -fPIC -64 ";
	F77FLAGS="-c  -64 -fPIC";
	# Include path
	CPPFLAGS="-DLOWER_UNDERSCORE"
	# Library path
	LDFLAGS=
	# Libraries
	LIBS=
	# Ranlib
	RANLIB="ar -ts";
fi

# Exporting variables
export CC F77 CFLAGS F77FLAGS CPPFLAGS LDFLAGS LIBS RANLIB DEBUG RELEASE BUILD LIBNAME DISTDIR ;
# make the distribution directory if it does not exists
mkdir ${DISTDIR}
# Compiling library

./compileCONMIN.sh 
# Copying library to $IGRID_HOME/lib/local
echo "============================================ ";
echo " Copying CONMIN Shared Object/Dynamic Library   to ${IGRID_HOME}/lib/local ";
echo " command = cp ${DISTDIR}/${LIBNAME} ${IGRID_HOME}/lib/local/";
echo "============================================ ";
cp  ${DISTDIR}/${LIBNAME}  ${IGRID_HOME}/lib/local/.

exit 0;
# EOF
