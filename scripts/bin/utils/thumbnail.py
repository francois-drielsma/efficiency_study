import sys
import subprocess

for image_name in sys.argv[1:]:
    image_first = image_name[:-4]
    subprocess.check_output(["convert", "-thumbnail", "x100", image_name, image_first+"_thumb.png"])
