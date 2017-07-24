if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /software/BES/bes6.6.4.p01/CMT/v1r20p20090520
endif
source ${CMTROOT}/mgr/setup.csh
set tempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set tempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt cleanup -csh -pack=ThresholdAlg -version=Threshold-00-00-01 -path=/home/lliu/boss/Workarea-6.6.4.p01 $* >${tempfile}; source ${tempfile}
/bin/rm -f ${tempfile}

