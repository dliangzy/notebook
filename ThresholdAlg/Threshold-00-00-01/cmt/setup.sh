# echo "Setting ThresholdAlg Threshold-00-00-01 in /home/lliu/boss/Workarea-6.6.4.p01"

if test "${CMTROOT}" = ""; then
  CMTROOT=/software/BES/bes6.6.4.p01/CMT/v1r20p20090520; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh

tempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then tempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt setup -sh -pack=ThresholdAlg -version=Threshold-00-00-01 -path=/home/lliu/boss/Workarea-6.6.4.p01  -no_cleanup $* >${tempfile}; . ${tempfile}
/bin/rm -f ${tempfile}

