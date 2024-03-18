#!/bin/bash
basedir="/home/jiexuan/slt/IRAK4/ledock"
ligdir="${basedir}/lig"

## lepro protein
# for pro in $(ls | grep pdb );do
#     proname=$(basename ${pro} .pdb)
#     mkdir ${proname}_dock
#     lepro ${pro}
#     mv dock.in ./${proname}_dock
#     mv pro.pdb ./${proname}_dock    
# done

##ligand split step1
# for ((j=0;j<5;j++));do    
#     mkdir ${ligdir}/"temp"${j}
# done
# wait
# echo "preparation is over"

# cd ${ligdir}
# ###ligand split step2
# a=1    
# for file in $(ls | grep mol2 );do
#     c=`expr $a % 5`
#     mv $file ${ligdir}/"temp"${c}/$file
#     a=`expr $a + 1`
# done
# b=`expr $a / 5`
# echo "spliting is over"

###ledock
# for ((j=0;j<5;j++));do
#     cp ${basedir}/${proname}_dock/dock.in ${ligdir}/"temp"${j}
#     cp ${basedir}/${proname}_dock/pro.pdb ${ligdir}/"temp"${j}
#     cd ${ligdir}/"temp"${j}
#    # {
#     for lig in $(ls | grep mol2 );do
#         ligname=$(basename ${lig} .mol2)
#         ls *.mol2>ligands
#         {    
#             ledock dock.in
#         }&
#         done
# #    }&
# done
# wait

###ledock analysis
for ((j=0;j<5;j++));do
    cd ${ligdir}/"temp"${j}
    for dok in $(ls | grep dok );do
        y=$(basename ${dok} .dok)
        a=$(awk -F ' ' '{print$8}' ${dok} | sed -n '2p')
        if [ $a ];then
            echo -n "${y} " >> ${basedir}/score.txt
            echo ${a} >> ${basedir}/score.txt
        fi
    done
    cd ../
done
exit 0