# echo "setup TwoProngAlg TwoProngAlg-00-01 in /home/liud/workarea665"

if test "${CMTROOT}" = ""; then
  CMTROOT=/software/BES/bes6.6.5/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtTwoProngAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtTwoProngAlgtempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt setup -sh -pack=TwoProngAlg -version=TwoProngAlg-00-01 -path=/home/liud/workarea665  -no_cleanup $* >${cmtTwoProngAlgtempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt setup -sh -pack=TwoProngAlg -version=TwoProngAlg-00-01 -path=/home/liud/workarea665  -no_cleanup $* >${cmtTwoProngAlgtempfile}"
  cmtsetupstatus=2
  /bin/rm -f ${cmtTwoProngAlgtempfile}
  unset cmtTwoProngAlgtempfile
  return $cmtsetupstatus
fi
cmtsetupstatus=0
. ${cmtTwoProngAlgtempfile}
if test $? != 0 ; then
  cmtsetupstatus=2
fi
/bin/rm -f ${cmtTwoProngAlgtempfile}
unset cmtTwoProngAlgtempfile
return $cmtsetupstatus

