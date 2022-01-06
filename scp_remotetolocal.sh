#!/bin/bash
passwd="*********\r"
file=$1
expect<<-END
spawn scp user@hostname:$file .
expect "password: "
send $passwd
expect eof
exit
END
