# ccx_prool
Prool's modifications of CalculiX CrunchiX (ccx)

download binaryes for windows
-----------------------------

http://files.calculix.kharkov.org/ccx_prool/4/ccx_5_jul_2015.zip

cygwin
------

For cygwin exe file need this DLLs from /cygwin/bin:

cyggcc_s-1.dll

cyggfortran-3.dll

cyggomp-1.dll

cygquadmath-0.dll

cygwin1.dll

Use native Linux libarpack (recommended for 32 bit Linuses)
--------------------------------------------------------------

for Centos:

yum install arpack-devel

yum install lapack-devel

yum install blas-devel

(or 'apt-get install libarpack2-dev' in Debian)

and use Makefile.arpack

DIRTY HACK only for 32 bit Linuses
----------------------------------

in file ARPACK/LAPACK/dlamch.f add LBETA=2 after string 'LBETA=C+QTR'

Note: NO WARRANTIES! It's FOOLISH hack bikoz i'm mad programmer

prool contacts
--------------

e-mail proolix@gmail.com

sites:

http://calculixforwin.com

http://calculix.kharkov.org

http://prool.kharkov.org
