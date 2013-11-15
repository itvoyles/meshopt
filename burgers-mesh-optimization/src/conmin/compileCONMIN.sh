#!/bin/sh
#############################################
# Compile CONMIN
#############################################
# Base directory
BASEDIR=`pwd`;

# Source directory
SRCDIR=${BASEDIR}/src;

# Kernel type
KERNEL=`uname`;
################################################################################
# Script Banner
echo "";
echo "============================================ ";
echo " Compile CONMIN Source Files  ";
echo "${F77} -c ${F77FLAGS} *.f"
echo "============================================ ";
echo "Machine kernel:   ${KERNEL}";
echo "Base directory:   ${BASEDIR}";
echo "Source directory: ${SRCDIR}";
echo "";
################################################################################
# remove existing library in the DISTDIR folder
rm ${DISTDIR}/*
# Source directory list
cd ${SRCDIR}
${F77}  ${F77FLAGS} *.f
echo "";
echo "============================================ ";
echo " Create  CONMIN Shared Object/Dynamic Library  ";
echo " Command = ${CC}  ${CFLAGS} ${LDFLAGS} ${SRCDIR}/*.o -o ${DISTDIR}/${LIBNAME} -lgfortran ";
echo "============================================ ";
${CC}  ${LDFLAGS} ${SRCDIR}/*.o -o ${DISTDIR}/${LIBNAME} -lgfortran
rm ${SRCDIR}/*.o
#
exit 0;
#EOF
