import os
import glob

cwd = os.getcwd()
source = os.path.join(cwd, "output/2017-02-7-Systematics-v1/*")
target = os.path.join(cwd, "output/2017-02-Systematics-5")

for source_filename in glob.glob(source):
    filename = os.path.split(source_filename)[1]
    target_filename = os.path.join(target, filename)
    try:
        #print source_filename, target_filename
        os.symlink(source_filename, target_filename)
    except OSError:
        print "Failed", filename
