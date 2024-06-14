# echo "cleanup SampleAlg SampleAlg-00-00-01 in /workfs/bes/leo591653959/BOSS705/workarea/Analysis"

if test "${CMTROOT}" = ""; then
  CMTROOT=/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtSampleAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtSampleAlgtempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt cleanup -sh -pack=SampleAlg -version=SampleAlg-00-00-01 -path=/workfs/bes/leo591653959/BOSS705/workarea/Analysis  $* >${cmtSampleAlgtempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt cleanup -sh -pack=SampleAlg -version=SampleAlg-00-00-01 -path=/workfs/bes/leo591653959/BOSS705/workarea/Analysis  $* >${cmtSampleAlgtempfile}"
  cmtcleanupstatus=2
  /bin/rm -f ${cmtSampleAlgtempfile}
  unset cmtSampleAlgtempfile
  return $cmtcleanupstatus
fi
cmtcleanupstatus=0
. ${cmtSampleAlgtempfile}
if test $? != 0 ; then
  cmtcleanupstatus=2
fi
/bin/rm -f ${cmtSampleAlgtempfile}
unset cmtSampleAlgtempfile
return $cmtcleanupstatus

