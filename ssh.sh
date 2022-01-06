#!/usr/bin/expect
spawn ssh user@hostname
expect "password: "
send "*******\r"
interact


