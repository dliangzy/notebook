#!/usr/bin/expect
spawn ssh lliu@ui06.lcg.ustc.edu.cn
expect "password: "
send "*******\r"
interact


