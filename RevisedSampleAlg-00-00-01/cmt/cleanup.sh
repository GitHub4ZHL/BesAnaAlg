# echo "cleanup RevisedSampleAlg RevisedSampleAlg-00-00-03 in /workfs2/bes/zhanghaolin/BOSS708/workarea/Analysis"

if test "${CMTROOT}" = ""; then
  CMTROOT=/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtRevisedSampleAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtRevisedSampleAlgtempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt cleanup -sh -pack=RevisedSampleAlg -version=RevisedSampleAlg-00-00-03 -path=/workfs2/bes/zhanghaolin/BOSS708/workarea/Analysis  $* >${cmtRevisedSampleAlgtempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt cleanup -sh -pack=RevisedSampleAlg -version=RevisedSampleAlg-00-00-03 -path=/workfs2/bes/zhanghaolin/BOSS708/workarea/Analysis  $* >${cmtRevisedSampleAlgtempfile}"
  cmtcleanupstatus=2
  /bin/rm -f ${cmtRevisedSampleAlgtempfile}
  unset cmtRevisedSampleAlgtempfile
  return $cmtcleanupstatus
fi
cmtcleanupstatus=0
. ${cmtRevisedSampleAlgtempfile}
if test $? != 0 ; then
  cmtcleanupstatus=2
fi
/bin/rm -f ${cmtRevisedSampleAlgtempfile}
unset cmtRevisedSampleAlgtempfile
return $cmtcleanupstatus

