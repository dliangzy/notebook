#!/bin/bash
#This script is used to copy file from localhost to remote host
passwd="********\r"
file=$1
expect<<-END
spawn scp $file user@hostname:/home/lliu/scp/
expect "password: "
send $passwd
expect eof
exit
END
