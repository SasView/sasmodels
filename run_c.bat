rm xlate.c
python sascomp cylinder,_cylpy -double!
gcc -c xlate.c 2> err.txt