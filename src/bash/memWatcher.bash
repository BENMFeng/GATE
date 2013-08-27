#!/bin/bash
##########################################################################
#  Copyright (c) 2011-2012 BENM(Binxiao) Feng                            #
#  All Rights Reserved                                                   #
#  Send all comments to BENM - BinxiaoFeng@gmail.com                     #
#                                                                        #
#  This program is free software: you can redistribute it and/or modify  #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation, either version 3 of the License, or     #
#  (at your option) any later version.                                   #
#  This program is distributed in the hope that it will be useful,       #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
#  GNU General Public License for more details.                          #
#  You should have received a copy of the GNU General Public License     #
#  along with this program.  If not, see <http://www.gnu.org/licenses/>. #
##########################################################################

interval=60
times=7200
minmem=2000
threshold=90
reducedspeed=100
alert_p=0
kill_p=0
linuxmark="Your"
while getopts "i:t:f:m:p:r:L:ak" optname
    do
        case "$optname" in
            "i")
                interval=$OPTARG
                ;;
            "t")
                times=$OPTARG
                ;;
            "f")
                maillist=$OPTARG
                ;;
            "m")
                minmem=$OPTARG
                ;;
            "p")
                threshold=$OPTARG
                ;;
            "r")
                reducedspeed=$OPTARG
                ;;
            "L")
                linuxmark=$OPTARG
                ;;
            "a")
                alert_p=1
                ;;
            "k")
                kill_p=1
                ;;
            "?")
                echo "Unknown option $OPTARG"
                ;;
            ":")
                echo "No argument value for option $OPTARG"
                ;;
            *)
            # Should not occur
                echo "Unknown error while processing options"
                ;;
            esac
        done

LANG=C
date=$(date +"%Y-%m-%d %H:%M:%S")
progname=`basename $0`
if [ "$maillist" == "" ];then
    echo "memWatcher -- This Bash shell script is designed to monitor memory of Linux server"
	echo "Author: BENM <binxiaofeng@gmail.com>"
    echo "Version: 2012-03-02 v1.6 Nightly"
	echo "Usage: $progname [Option]"
    echo "  -f <file>   email_list.conf"
    echo "  -a          alert to users[recommand]"
    echo "  -k          kill the largest ram-mem consumed job"
    echo "  -i <int>    interval_time/sec, default: 60"
    echo "  -t <int>    times, default: 7200"
    echo "  -m <int>    minimum memory, Mb: 2000"
    echo "  -p <int>    %MEM theshold: 90"
    echo "  -r <int>    memory reduced speed, Mb/s: 100"
    echo "  -L <str>    Linux mark: Your"
    echo "Note: run as root!"
    echo "Example: $progname -f email_list.conf -i 60 -t 4320 -L DELL -a -k"
    echo ""
    echo "==> email_list.conf Example <=="
    echo "test  binxiaofeng@gmail.com   2Gb"
	exit
fi

users=($(awk '{print $1}' ${maillist}))
email=($(awk '{print $2}' ${maillist}))
memth=($(perl -ne 'chomp;my @t=split/\s+/,$_;if($t[2]=~/(\d+)G/i){print $t[2]*1000,"\n";}elsif($t[2]=~/(\d+)M/i){print "$t[2]\n"}elsif($t[2]=~/(\d+)k/i){print $t[2]/1000,"\n";}elsif($t[2]=~/^\d+$/){print "$t[2]\n";}' ${maillist}))
root_mail=($(awk '$1=="root"&&$2~/\@/{print $2}' ${maillist}))
if [ "$root_mail" == "" ];then
    echo "No root mail or illegal email address!"
    exit
fi

check_user(){
    timemark=`date | sed 's/ /-/g'`
    B=("FENG" 0)
    while [ ${B[0]} != "BENM" ];do
        B=($(ps auxf |awk '{u[$1]+=$4}END{for (m in u){if (u[m]>0){print m,u[m]}}print "BENM\t0"}'|sort -nr -k 2|head -1))
        cmp=`echo "if(${B[1]} > $threshold) print "1" else print "0"" | bc`
        if [ $cmp == 1 ];then
            for (( i = 0; i < ${#users[@]}; i++ )); do
                if [ ${B[0]} ==  ${users[$i]} ];then
                    if [ $alert_p == 1 ];then
                        process_head=$(ps auxf|head -1)
                        user_process=$(ps auxf|grep ${B[0]})
                        kill_pid=$(ps auxf |grep ${B}|sort -nr -k4|awk '{print "kill PID:"$2}')
                        echo -e "${timemark} alert to ${email[$i]}\n${process_head}\n${user_process}\n"
                        echo -e "Warning! Your job consume too much memory, it had to be killed! Date: ${timemark}\n${kill_pid}\n${process_head}\n${user_process}" |mail -s "Jobs killed!" ${email[$i]}
                    fi
                    if [ $kill_p == 1 ];then
                        echo ${timemark}
                        ps auxf |grep ${B[0]}|sort -nr -k4|head -1|awk '{print "kill PID:"$2}'
                        ps auxf |grep ${B[0]}|sort -nr -k4|head -1|awk '{system("kill -9 "$2)}'
                    fi
                fi
            done
        fi
    done
## users specific physical memory check
    for (( i = 0; i < ${#users[@]}; i++ )); do
        if [ -n "${memth[$i]}" ];then
            check_user_mem=($(ps auxf |grep ${users[$i]} |awk '{u[$1]+=$6/1000}END{for (m in u){if (u[m]>0){print m,int(u[m])}}print "BENM\t0"}'|sort -nr -k 2|head -1))
            cmp=`echo "if(${check_user_mem[1]} > ${memth[$i]}) print "1" else print "0"" | bc`
            if [ $cmp == 1 ];then
                if [ $kill_p == 1 ];then
                    echo ${timemark}
                    ps auxf |grep ${check_user_mem[0]}|sort -nr -k4|head -1|awk '{print "kill PID:"$2}'
                    ps auxf |grep ${check_user_mem[0]}|sort -nr -k4|head -1|awk '{system("kill -9 "$2)}'
                fi
                if [ $alert_p == 1 ];then
                    process_head=$(ps auxf|head -1)
                    user_process=$(ps auxf|grep ${check_user_mem[0]})
                    echo -e "${timemark} alert to ${email[$i]}\n${process_head}\n${user_process}\n"
                    echo -e "Warning! Your job consume too much memory, it had to be killed! Date: ${timemark}\n\n${process_head}\n${user_process}" |mail -s "Jobs killed!" ${email[$i]}
                fi
            fi
        fi
    done
}

mem_monitor(){

    mem_free=$(grep MemFree /proc/meminfo|awk '{print int($2/1000)}')
    mem_resource=$(free -m)
    mem_reduce=$(free -m -s1 -c 5|awk 'BEGIN{p=1;}{if($1~/Mem/){t=$4;if ((p>1)&&(t-p)>=0){print (t-p)}p=$4};}END{print 1}'|sort -nr -k 1|head -1)
    mem_check=$(free -m -s1 -c 5|awk 'BEGIN{p=1;}{if($1~/Mem/){t=$4;if ((p>1)&&(t-p)>=0){print "Mem reduced speed:",(t-p),"Mb/s Temporary Memory free: ",t,"Mb"}p=$4}}'|sort -nr -k 4|head -1)
    C=$(ps auxf |sort -nr -k 4|awk '{u[$1]+=$4}END{for (m in u){if (u[m]>0){print m,u[m]}}print "BENM\t0"}'|sort -nr -k 2|awk '{print $1}'|head -1)
    cmp=`echo "if($mem_free < $minmem) print "1" else print "0"" | bc`
    if [ $cmp == 1 ];then
        for (( i = 0; i < ${#users[@]}; i++ )); do
            if [ ${C} ==  ${users[$i]} ];then
                echo -e "${date} alert to ${users[$i]}: ${linuxmark} Linux server ran out of memory: ${mem_free}Mb free!\n${mem_check}\n\n${mem_resource}\n"
                echo -e "${linuxmark} Linux server ran out of memory: ${mem_free}Mb free! Date: ${date}\n\n${mem_check}\n\n${mem_resource}" | mail -s "Mem Alert! Do not reply!" ${email[$i]}
            fi
        done
        echo -e "${linuxmark} Linux server ran out of memory: ${mem_free}Mb free! Date: ${date}\n\n${mem_check}\n\n${mem_resource}" | mail -s "Mem Alert! Do not reply!" ${root_mail}
        echo -e "${date} alert to root: ${linuxmark} Linux server ran out of memory: ${mem_free}Mb free!\n${mem_check}\n\n${mem_resource}\n"
    fi
    
    cmp=`echo "if($mem_reduce > $reducedspeed) print "1" else print "0"" | bc`
    if [ $cmp == 1 ];then
        for (( i = 0; i < ${#users[@]}; i++)); do
            if [ ${C} ==  ${users[$i]} ];then
                echo -e "${linuxmark} Linux server memory dumped too fast! Date: ${date}\n\n${mem_check} Temporary memory: ${mem_free}Mb\n\n${mem_resource}" |mail -s "Mem Dump! Do not reply!" ${email[$i]}
            fi
        done
        echo -e "${linuxmark} Linux server memory dumped too fast! Date: ${date}\n\n${mem_check}\n\n${mem_resource}" |mail -s "Mem Dump! Do not reply!" ${root_mail}
        echo -e "${date} alert to root: ${linuxmark} Linux server memory dumped too fast!\n${mem_check}\n\n${mem_resource}\n"
    fi
    
    check_user

}

for (( i = 0; i < ${times}; i++)); do 
    mem_monitor
    sleep ${interval}
done

