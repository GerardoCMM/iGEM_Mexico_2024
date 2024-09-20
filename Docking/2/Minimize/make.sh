for f in *.sdf
do
	python make.py "$f"
done

for f in *_minimized.mol
do
	obabel -imol "$f" -opdb -O "${f/%.mol/.pdb}"
done

obminimize -ff UFF -o sdf porphyrin.sdf > porphyrin_minimized.sdf
