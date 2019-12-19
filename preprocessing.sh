# Written on GMT v6 (http://gmt.soest.hawaii.edu/)

# Resamples rivers and calculates azimuths using gmt mapproject and find_near.py
# Assumes rivers are imported with longitude and latitude in decimal degrees

gmt mapproject -Gk river.ll > distances

python find_near.py # Calls external python script to find next point on river

# Computes resulting file of form long, lat, distance, azimuth

awk '{print $2, $3, $4}' < distances.out | gmt mapproject -Af | awk '{if (NR>1) print $0}' > riv_equid_az.lldv
rm distances.out
# Example output of this is given in repository github.com/AlexLipp/directional-wavelets
# as file `colorado.dat`.
