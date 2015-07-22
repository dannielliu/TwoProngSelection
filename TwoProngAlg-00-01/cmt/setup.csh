# echo "setup TwoProngAlg TwoProngAlg-00-01 in /home/liud/workarea665"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /software/BES/bes6.6.5/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtTwoProngAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtTwoProngAlgtempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt setup -csh -pack=TwoProngAlg -version=TwoProngAlg-00-01 -path=/home/liud/workarea665  -no_cleanup $* >${cmtTwoProngAlgtempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt setup -csh -pack=TwoProngAlg -version=TwoProngAlg-00-01 -path=/home/liud/workarea665  -no_cleanup $* >${cmtTwoProngAlgtempfile}"
  set cmtsetupstatus=2
  /bin/rm -f ${cmtTwoProngAlgtempfile}
  unset cmtTwoProngAlgtempfile
  exit $cmtsetupstatus
endif
set cmtsetupstatus=0
source ${cmtTwoProngAlgtempfile}
if ( $status != 0 ) then
  set cmtsetupstatus=2
endif
/bin/rm -f ${cmtTwoProngAlgtempfile}
unset cmtTwoProngAlgtempfile
exit $cmtsetupstatus

