for f in *.png
do
	mogrify -crop 900x900+510+95 "$f"
done