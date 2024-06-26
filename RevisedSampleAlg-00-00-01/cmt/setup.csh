# echo "setup RevisedSampleAlg RevisedSampleAlg-00-00-01 in /workfs2/bes/zhanghaolin/BOSS708/workarea/Analysis"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtRevisedSampleAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtRevisedSampleAlgtempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt setup -csh -pack=RevisedSampleAlg -version=RevisedSampleAlg-00-00-01 -path=/workfs2/bes/zhanghaolin/BOSS708/workarea/Analysis  -no_cleanup $* >${cmtRevisedSampleAlgtempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt setup -csh -pack=RevisedSampleAlg -version=RevisedSampleAlg-00-00-01 -path=/workfs2/bes/zhanghaolin/BOSS708/workarea/Analysis  -no_cleanup $* >${cmtRevisedSampleAlgtempfile}"
  set cmtsetupstatus=2
  /bin/rm -f ${cmtRevisedSampleAlgtempfile}
  unset cmtRevisedSampleAlgtempfile
  exit $cmtsetupstatus
endif
set cmtsetupstatus=0
source ${cmtRevisedSampleAlgtempfile}
if ( $status != 0 ) then
  set cmtsetupstatus=2
endif
/bin/rm -f ${cmtRevisedSampleAlgtempfile}
unset cmtRevisedSampleAlgtempfile
exit $cmtsetupstatus

