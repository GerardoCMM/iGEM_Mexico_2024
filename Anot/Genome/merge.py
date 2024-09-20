import sys
import os

original_stdout = sys.stdout

with open('merged.txt', 'w') as f:
    sys.stdout = f # Change the standard output to the file we created.
    for i in range(9):
    	with open(f"result_{i+1}.txt") as file:
    		for line in file:
    			print(line.rstrip())
    sys.stdout = original_stdout # Reset the standard output to its original value

"""with open('merged_ko.txt', 'w') as f:
    sys.stdout = f # Change the standard output to the file we created.
    for i in range(9):
        with open(f"result_ko_{i+1}.txt") as file:
            for line in file:
                print(line.rstrip())
    sys.stdout = original_stdout # Reset the standard output to its original value"""

