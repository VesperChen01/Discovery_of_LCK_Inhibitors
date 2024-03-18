#!/bin/bash

basedir="/home/jiexuan/slt/IRAK4/autodock-gpu_nrun100/"
ligdir="${basedir}/lig"

for ((j=0;j<8;j++));do
    cd ${ligdir}/"temp"${j}
    {
    for dlg in $(ls | grep dlg );do
        y=$(basename ${dlg} .dlg)
        echo -n "${y}," >> ../../"score"${j}.txt
        a=$(awk -F ' ' '{print$4}' ${dlg} | tac | sed -n '103p')
        if [ $a ];then
        echo ${a} >> ../../"score"${j}.txt
        fi
    done
    }&
done
wait
cd ../..
cat *.txt > score.txt
exit


# cd result
#     for dlg in $(ls | grep dlg );do
#         y=$(basename ${dlg} .dlg)
#         echo -n "${y}," >> ../score2.txt
#         a=$(awk -F ' ' '{print$4}' ${dlg} | tac | sed -n '6p')
#         if [ $a ];then
#         echo ${a} >> ../score2.txt
#         fi
#     done
# wait
# cat *.txt > score.txt
# exit