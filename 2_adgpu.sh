#!/bin/bash
####===========================================================####
####        File Name: adgpu.sh                                ####
####        Author: szy && gj                                  ####
####        E-mail: shenzheyuan1@163.com                       ####
####        Date: 2020-10-15 16.03                             ####
####        Last Modified: 2020-10-18 01.03                    ####
####        Version: 2.0                                       ####
####===========================================================####
####prepare protein && protein autgrid
####change thing : basedir/ pdbname (vrk3)/ gridcenter /

tooldir="/opt/mgltools/1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24"
range="50,50,50"

basedir="/home/jiexuan/slt/IRAK4/autodock-gpu_nrun100/"
ligdir="${basedir}/lig"
center="1.9,-15.32,-61.42"

# cd lig
# for lig in $(ls | grep mol2 );do
#     ligname=$(basename ${lig} .mol2)
# #    /opt/mgltools/1.5.7/bin/pythonsh $tooldir/prepare_ligand4.py -l ${lig} -o ${ligname}.pdbqt -A checkhydrogens
# done
# wait
# cd ..

cd lig
for lig in $(ls | grep pdbqt );do
    ligname=$(basename ${lig} .pdbqt)
#    /opt/mgltools/1.5.7/bin/pythonsh $tooldir/prepare_ligand4.py -l ${lig} -o ${ligname}.pdbqt -A checkhydrogens
done
wait
cd ..

for pdb in $(ls | grep pdb );do
    pdbname=$(basename ${pdb} .pdb)
    /opt/mgltools/1.5.7/bin/pythonsh $tooldir/prepare_receptor4.py -r $pdb -o $pdbname.pdbqt -A checkhydrogens -U nphs_lps_waters
done

# cd lig2
#     for lig in $(ls | grep pdbqt );do
#         # read -u 10
#         {
#             ligname=$(basename ${lig} .pdbqt)
#             # python /opt/modules/mgltools/1.5.6@096b3c/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l ${lig} -o ${ligname}.pdbqt -A checkhydrogens
#             # echo >&6
#         }
# done

# cd ..

/opt/mgltools/1.5.7/bin/pythonsh ${tooldir}/prepare_gpf4.py -r $pdbname.pdbqt -o $pdbname.gpf -p npts="${range}" -p gridcenter="${center}" -p ligand_types="HD,HS,C,A,N,NA,OA,F,SA,S,P,Cl,Br,I"
/opt/autodock/gpu-1.5.3/bin/autogrid4 -p $pdbname.gpf

for ((j=0;j<8;j++));do    ###ligand split step1
    mkdir ${ligdir}/"temp"${j}
done
wait
echo "preparation is over"

cd ${ligdir}

a=1    ###ligand split step2
for file in $(ls | grep pdbqt );do
    c=`expr $a % 8`
    cp $file ${ligdir}/"temp"${c}/$file
    a=`expr $a + 1`
done
b=`expr $a / 8`
echo "spliting is over (8 * ${b})"

echo "over"
cd ${basedir}

for ((j=0;j<8;j++));do
    cd ${ligdir}/"temp"${j}
   # {
    for lig in $(ls | grep pdbqt );do
        ligname=$(basename ${lig} .pdbqt)
        /opt/autodock/gpu-1.5.3/bin/autodock_gpu_128wi -ffile ${basedir}/$pdbname.maps.fld -lfile ./${ligname}.pdbqt -devnum 1 -resnam ${ligname} -nrun 100
        done
#    }&
done
wait


# cd ${basedir}/result
# for dlg in $(ls | grep dlg );do
# y=$(basename ${dlg} .dlg)
# echo -n "${y}," >> ../score.txt
# a=$(awk -F ' ' '{print$4}' ${dlg} | tac | sed -n '4p')
# if [ $a ];then
# echo ${a} >> ../score.txt
# fi
# done

##for dlg in $(ls | grep dlg );do
##{
##dlgname=$(basename ${dlg} .dlg)
##pythonsh /opt/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/write_conformations_from_dlg.py -d ${dlgname}.dlg
##}
##done
exit


