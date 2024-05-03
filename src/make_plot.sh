#!/bin/bash

if [ $# -ne 2 ]; then
    echo "make_plot MODEL NMESH" 1>&2
    exit 1
fi

if [ $1 = "sod" ]; then
    CFLAG_INIT="INIT=\"-D__SOD_INIT__\""
    REF_FILE="sod_ref.dat"
    PREFIX="sod_"
    FIG_DIR="sod_$2"
fi
if [ $1 = "lax" ]; then
    CFLAG_INIT="INIT=\"-D__LAX_INIT__\""
    REF_FILE="lax_ref.dat"
    PREFIX="lax_"
    FIG_DIR="lax_$2"
fi
if [ $1 = "shu" ]; then
    CFLAG_INIT="INIT=\"-D__SHU_INIT__\""
    REF_FILE="shu_ref.dat"
    PREFIX="shu_"
    FIG_DIR="shu_$2"
fi

if [ $1 = "blanc" ]; then
    CFLAG_INIT="INIT=\"-D__LE_BLANC_INIT__\""
    REF_FILE="le_blanc_ref.dat"
    PREFIX="le_blanc_"
    FIG_DIR="le_blanc_$2"
fi
CFLAG_NMESH="NMESH=\"-D__NMESH__=$2\""
mkdir $FIG_DIR

make distclean
make $CFLAG_INIT $CFLAG_NMESH fluid_1d
./fluid_1d > hoge
python3 plot.py hoge $REF_FILE
mv dens.png ${FIG_DIR}/${PREFIX}FVS_dens.png
mv velx.png ${FIG_DIR}/${PREFIX}FVS_velx.png
mv pres.png ${FIG_DIR}/${PREFIX}FVS_pres.png

make distclean
make $CFLAG_INIT $CFLAG_NMESH FLUX="-D__FVS_MUSCL__" fluid_1d
./fluid_1d > hoge
python3 plot.py hoge $REF_FILE
mv dens.png ${FIG_DIR}/${PREFIX}FVS_MUSCL_dens.png
mv velx.png ${FIG_DIR}/${PREFIX}FVS_MUSCL_velx.png
mv pres.png ${FIG_DIR}/${PREFIX}FVS_MUSCL_pres.png

make distclean
make $CFLAG_INIT $CFLAG_NMESH FLUX="-D__FVS_MUSCL2__" fluid_1d
./fluid_1d > hoge
python3 plot.py hoge $REF_FILE
mv dens.png ${FIG_DIR}/${PREFIX}FVS_MUSCL2_dens.png
mv velx.png ${FIG_DIR}/${PREFIX}FVS_MUSCL2_velx.png
mv pres.png ${FIG_DIR}/${PREFIX}FVS_MUSCL2_pres.png

make distclean
make $CFLAG_INIT $CFLAG_NMESH FLUX="-D__FVS_MP5__" fluid_1d
./fluid_1d > hoge
python3 plot.py hoge $REF_FILE
mv dens.png ${FIG_DIR}/${PREFIX}FVS_MP5_wo_TVDRK3_dens.png
mv velx.png ${FIG_DIR}/${PREFIX}FVS_MP5_wo_TVDRK3_velx.png
mv pres.png ${FIG_DIR}/${PREFIX}FVS_MP5_wo_TVDRK3_pres.png

make distclean
make $CFLAG_INIT $CFLAG_NMESH FLUX="-D__FVS_MP5__ -D__TVD_RK3__" fluid_1d
./fluid_1d > hoge
python3 plot.py hoge $REF_FILE
mv dens.png ${FIG_DIR}/${PREFIX}FVS_MP5_dens.png
mv velx.png ${FIG_DIR}/${PREFIX}FVS_MP5_velx.png
mv pres.png ${FIG_DIR}/${PREFIX}FVS_MP5_pres.png

make distclean
make $CFLAG_INIT $CFLAG_NMESH FLUX="-D__FVS_MP5_2__" fluid_1d
./fluid_1d > hoge
python3 plot.py hoge $REF_FILE
mv dens.png ${FIG_DIR}/${PREFIX}FVS_MP5_2_wo_TVDRK3_dens.png
mv velx.png ${FIG_DIR}/${PREFIX}FVS_MP5_2_wo_TVDRK3_velx.png
mv pres.png ${FIG_DIR}/${PREFIX}FVS_MP5_2_wo_TVDRK3_pres.png

make distclean
make $CFLAG_INIT $CFLAG_NMESH FLUX="-D__FVS_MP5_2__ -D__TVD_RK3__" fluid_1d
./fluid_1d > hoge
python3 plot.py hoge $REF_FILE
mv dens.png ${FIG_DIR}/${PREFIX}FVS_MP5_2_dens.png
mv velx.png ${FIG_DIR}/${PREFIX}FVS_MP5_2_velx.png
mv pres.png ${FIG_DIR}/${PREFIX}FVS_MP5_2_pres.png

make distclean 
make $CFLAG_INIT $CFLAG_NMESH FLUX="-D__FVS_MP5SL__" fluid_1d
./fluid_1d > hoge
python3 plot.py hoge $REF_FILE
mv dens.png ${FIG_DIR}/${PREFIX}FVS_MP5SL_dens.png
mv velx.png ${FIG_DIR}/${PREFIX}FVS_MP5SL_velx.png
mv pres.png ${FIG_DIR}/${PREFIX}FVS_MP5SL_pres.png

make distclean 
make $CFLAG_INIT $CFLAG_NMESH FLUX="-D__AUSMP__" fluid_1d
./fluid_1d > hoge
python3 plot.py hoge $REF_FILE
mv dens.png ${FIG_DIR}/${PREFIX}AUSMP_dens.png
mv velx.png ${FIG_DIR}/${PREFIX}AUSMP_velx.png
mv pres.png ${FIG_DIR}/${PREFIX}AUSMP_pres.png

make distclean 
make $CFLAG_INIT $CFLAG_NMESH FLUX="-D__AUSMP_MUSCL__" fluid_1d
./fluid_1d > hoge
python3 plot.py hoge $REF_FILE
mv dens.png ${FIG_DIR}/${PREFIX}AUSMP_MUSCL_dens.png
mv velx.png ${FIG_DIR}/${PREFIX}AUSMP_MUSCL_velx.png
mv pres.png ${FIG_DIR}/${PREFIX}AUSMP_MUSCL_pres.png

make distclean 
make $CFLAG_INIT $CFLAG_NMESH FLUX="-D__AUSMP_MP5__" fluid_1d
./fluid_1d > hoge
python3 plot.py hoge $REF_FILE
mv dens.png ${FIG_DIR}/${PREFIX}AUSMP_MP5_dens.png
mv velx.png ${FIG_DIR}/${PREFIX}AUSMP_MP5_velx.png
mv pres.png ${FIG_DIR}/${PREFIX}AUSMP_MP5_pres.png

make distclean 
make $CFLAG_INIT $CFLAG_NMESH FLUX="-D__HLL__" fluid_1d
./fluid_1d > hoge
python3 plot.py hoge $REF_FILE
mv dens.png ${FIG_DIR}/${PREFIX}HLL_dens.png
mv velx.png ${FIG_DIR}/${PREFIX}HLL_velx.png
mv pres.png ${FIG_DIR}/${PREFIX}HLL_pres.png

make distclean 
make $CFLAG_INIT $CFLAG_NMESH FLUX="-D__HLLC__" fluid_1d
./fluid_1d > hoge
python3 plot.py hoge $REF_FILE
mv dens.png ${FIG_DIR}/${PREFIX}HLLC_dens.png
mv velx.png ${FIG_DIR}/${PREFIX}HLLC_velx.png
mv pres.png ${FIG_DIR}/${PREFIX}HLLC_pres.png

make distclean 
make $CFLAG_INIT $CFLAG_NMESH FLUX="-D__HLLC_MUSCL__" fluid_1d
./fluid_1d > hoge
python3 plot.py hoge $REF_FILE
mv dens.png ${FIG_DIR}/${PREFIX}HLLC_MUSCL_dens.png
mv velx.png ${FIG_DIR}/${PREFIX}HLLC_MUSCL_velx.png
mv pres.png ${FIG_DIR}/${PREFIX}HLLC_MUSCL_pres.png

make distclean 
make $CFLAG_INIT $CFLAG_NMESH FLUX="-D__HLLC_MP5__ -D__TVD_RK3__" fluid_1d
./fluid_1d > hoge
python3 plot.py hoge $REF_FILE
mv dens.png ${FIG_DIR}/${PREFIX}HLLC_MP5_dens.png
mv velx.png ${FIG_DIR}/${PREFIX}HLLC_MP5_velx.png
mv pres.png ${FIG_DIR}/${PREFIX}HLLC_MP5_pres.png

