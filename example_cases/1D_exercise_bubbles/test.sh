#!/bin/bash
curdir=$(pwd)

nbs=( 1 5 9 )
nxs=( 1000 100000 1000000 )
gpuswitch=( 0 1 )

currswitch=0
for gpu in "${gpuswitch[@]}"; do
    if [[ $gpu -eq 0 ]]; then
        if [[ currswitch -ne 1 ]]; then  
            cd ../..
            make clean
            make -j Makefile.user.gpu
            currswitch=1
            cd $curdir 
        fi
    else
        if [[ currswitch -ne 2 ]]; then  
            make clean
            make -j Makefile.user.cpu
            currswitch=2
            cd $curdir 
        fi
    fi

    for nb in "${nbs[@]}"; do
        for nx in "${nxs[@]}"; do
            echo Nb: $nb and Nx: $nx
            rm -f input_new.py

            if [[ $nb -eq 1 ]]; then
                echo Mono
                sed -e "s/XXNX/$nx/g" input_mono.py > input_new.py
            else
                echo Poly
                sed -e "s/XXNB/$nb/g" -e "s/XXNX/$nx/g" input_poly.py > input_new.py
            fi

            ./input_new.py pre_process
            ./input_new.py simulation
        done
    done
done
