Code repository for NTNU Trondheim's iGEM18 team, Film Fighters.

This code is based upon an article by Fozard et al. (2012) about inhibition of
quorom sensing. The inhibitor molecule has been removed, but otherwise this
code should be similar.

An initial prototype was written in python3, while the newest and updated
version is in the fortran/ folder. Make sure to make a folder called data/
before you run any tests, and that there are no files of the type ---.csv in
the folder before the program runs. If you run a bash-shell, then simply run
fortran/run.sh to test the code. gfortran must be installed for this to work.
After running run.sh, you can run fortran/data/sum.sh to get a processed
data-file out, sum_proc.csv.

Plotting is easily done in gnuplot with the following code snippets,

set datafile separator ","

plot "sum_proc.csv" using 1:n with lines,

where n is an integer determining what you want to plot.

2) Substrate concentration

3) Quorom Sensing molecule concentration

4) Biomass

5) Activated bacteria

6) #EPS particles

7) EPS concentration (ignore this)
