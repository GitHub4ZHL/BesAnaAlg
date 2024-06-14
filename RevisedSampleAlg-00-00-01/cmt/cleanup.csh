# echo "cleanup SampleAlg SampleAlg-00-00-01 in /workfs/bes/leo591653959/BOSS705/workarea/Analysis"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtSampleAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtSampleAlgtempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt cleanup -csh -pack=SampleAlg -version=SampleAlg-00-00-01 -path=/workfs/bes/leo591653959/BOSS705/workarea/Analysis  $* >${cmtSampleAlgtempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt cleanup -csh -pack=SampleAlg -version=SampleAlg-00-00-01 -path=/workfs/bes/leo591653959/BOSS705/workarea/Analysis  $* >${cmtSampleAlgtempfile}"
  set cmtcleanupstatus=2
  /bin/rm -f ${cmtSampleAlgtempfile}
  unset cmtSampleAlgtempfile
  exit $cmtcleanupstatus
endif
set cmtcleanupstatus=0
source ${cmtSampleAlgtempfile}
if ( $status != 0 ) then
  set cmtcleanupstatus=2
endif
/bin/rm -f ${cmtSampleAlgtempfile}
unset cmtSampleAlgtempfile
exit $cmtcleanupstatus

