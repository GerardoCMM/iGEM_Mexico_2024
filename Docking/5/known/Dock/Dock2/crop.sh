for f in *.png
do
	mogrify -crop 800x800+560+120 "$f"
done