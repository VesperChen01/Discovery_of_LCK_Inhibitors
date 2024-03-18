#!/bin/bash
## module load mgltools/1.5.7 
## module load vina/1.2.3
## write config_XXXX.txt

# prepare ligands
cd lig
for lig in $(ls | grep mol2 );do
    ligname=$(basename ${lig} .mol2)
    # /opt/mgltools/1.5.7/bin/pythonsh /opt/mgltools/1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -l ${lig} -o ${ligname}.pdbqt -A checkhydrogens
done
wait
cd ..

# prepare protein
cd ./protein
for ipdb in $(ls | grep pdb );do
    ipdbname=$(basename ${ipdb} .pdb)
    pythonsh /opt/mgltools/1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_pdb_split_alt_confs.py -r ${ipdb}
    sleep 1s
    if [ -f "./${ipdbname}_A.pdb" ];then
        rm ${ipdb}
        cat config_${ipdbname}.txt>>config_${ipdbname}_A.txt
        cat config_${ipdbname}.txt>>config_${ipdbname}_B.txt
        rm config_${ipdbname}.txt
    else
        echo "${ipdbname}_A.pdb does not exist"
    fi
done
wait

for vpdb in $(ls | grep pdb );do
    vpdbname=$(basename ${vpdb} .pdb)
    mkdir ../${vpdbname}_out
    pythonsh /opt/mgltools/1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r ${vpdb} -o ${vpdbname}.pdbqt -A checkhydrogens -U nphs_lps_waters
    echo "protein trans over"
    ## vina_docking
    cd ../lig
    for vlig in $(ls | grep pdbqt );do
        vligname=$(basename ${vlig} .pdbqt)
        vina --config ../protein/config_${vpdbname}.txt --receptor ../protein/${vpdbname}.pdbqt --ligand ./${vligname}.pdbqt --out ../${vpdbname}_out/${vligname}.pdbqt 
    done
    wait
    ## vina analysis
    cd ../${vpdbname}_out
    for pdbqt in $(ls | grep pdbqt );do
        filename=$(basename ${pdbqt} .pdbqt)
        echo -n "${filename}," >> ../${vpdbname}_score.txt
        a=$(sed -n 2p ${pdbqt})
        echo ${a} >> ../${vpdbname}_score.txt
        echo "analyzing $pdbname is over" 
        echo "${pdbqt}"
    done
    wait
    cd ../protein
done