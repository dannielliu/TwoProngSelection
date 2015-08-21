if test "${CMTROOT}" = ""; then
  CMTROOT=/software/BES/bes6.6.4.p01/CMT/v1r20p20090520; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
tempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then tempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt cleanup -sh -pack=TwoProngAlg -version=TwoProngAlg-00-01 -path=/home/liud/workarea $* >${tempfile}; . ${tempfile}
/bin/rm -f ${tempfile}

