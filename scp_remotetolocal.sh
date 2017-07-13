#!/bin/bash
passwd="LL@ustc0629\r"
file=$1
expect<<-END
spawn scp lliu@ui06.lcg.ustc.edu.cn:$file .
expect "password: "
send $passwd
expect eof
exit
END
