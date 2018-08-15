for i in *.csv; do sed --quiet '2p' $i; done | sort -n | cut -d, -f1,5- >> proc.csv
