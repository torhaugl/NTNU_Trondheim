gfortran csv_file.f90 simulate.f90 -O3
rm *.mod
directory=$(date '+%y-%m-%d-')$1
mkdir data/${directory}
mv a.out data/${directory}/
cp data/plot.gnuplot data/${directory}/
cp data/sum.sh data/${directory}/
cd data/$directory
mkdir data
./a.out &&
cp sum.sh data/
cp plot.gnuplot data/
cd data
./sum.sh &&
./plot.gnuplot
mv *.png ../
cd ../../..
