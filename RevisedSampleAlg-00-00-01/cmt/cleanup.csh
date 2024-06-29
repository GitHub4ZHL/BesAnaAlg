# echo "cleanup RevisedSampleAlg RevisedSampleAlg-00-00-03 in /workfs2/bes/zhanghaolin/BOSS708/workarea/Analysis"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtRevisedSampleAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtRevisedSampleAlgtempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt cleanup -csh -pack=RevisedSampleAlg -version=RevisedSampleAlg-00-00-03 -path=/workfs2/bes/zhanghaolin/BOSS708/workarea/Analysis  $* >${cmtRevisedSampleAlgtempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt cleanup -csh -pack=RevisedSampleAlg -version=RevisedSampleAlg-00-00-03 -path=/workfs2/bes/zhanghaolin/BOSS708/workarea/Analysis  $* >${cmtRevisedSampleAlgtempfile}"
  set cmtcleanupstatus=2
  /bin/rm -f ${cmtRevisedSampleAlgtempfile}
  unset cmtRevisedSampleAlgtempfile
  exit $cmtcleanupstatus
endif
set cmtcleanupstatus=0
source ${cmtRevisedSampleAlgtempfile}
if ( $status != 0 ) then
  set cmtcleanupstatus=2
endif
/bin/rm -f ${cmtRevisedSampleAlgtempfile}
unset cmtRevisedSampleAlgtempfile
exit $cmtcleanupstatus

