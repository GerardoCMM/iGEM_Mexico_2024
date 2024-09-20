import os

os.system("mkdir images")

for i in range(25):
	os.system(f"cp ./{i+1}/circos.png ./images/{i+1}.png")

