for i in *.csv
do
   # Format: dt,cs,cq,biomass,up,eps_count,eps_amount
   cut -d, -f1 $i | head -n 2 | tail -n 1 | tr '\n' ',' >> sum.csv
   cut -d, -f5 $i | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' | paste -sd+ | bc | tr '\n' ',' >> sum.csv
   cut -d, -f6 $i | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' | paste -sd+ | bc | tr '\n' ',' >> sum.csv
   cut -d, -f7 $i | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' | paste -sd+ | bc | tr '\n' ',' >> sum.csv
   cut -d, -f8 $i | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' | paste -sd+ | bc | tr '\n' ',' >> sum.csv
   cut -d, -f9 $i | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' | paste -sd+ | bc | tr '\n' ',' >> sum.csv
   cut -d, -f10 $i | sed -E 's/([+-]?[0-9.]+)[eE]\+?(-?)([0-9]+)/(\1*10^\2\3)/g' | paste -sd+ | bc  >> sum.csv
done
sort -n sum.csv >> sum_proc.csv
