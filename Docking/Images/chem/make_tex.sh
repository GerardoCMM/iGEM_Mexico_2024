for f in *.smi
do

	mol2chemfig -zwof  "$f"  >  "${f/%.smi/.tex}"

	obabel -ismi "$f" -opng -O "${f/%.smi/.png}" -xp 1000

done