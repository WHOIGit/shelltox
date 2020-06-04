
cd alphagamma

time python alphagamma.py alphagamma.ini --src ../data --alpha 0 0.05 0.001 --gammas 0 0.2 0.005 -o output/alphagamma_fullscan.csv -m 1 2 3

cd ..
