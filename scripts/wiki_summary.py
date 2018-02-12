import glob

def do_summary(parent_folder):
   fout = open(parent_folder+"/wiki_summary.txt", "w")
   target = parent_folder+"/*"
   file_list = sorted(glob.glob(target), key = lambda name: name[-5:])
   print "Found", len(file_list), "searching target", target
   print file_list[0:3], "...", file_list[-3:]
   for item in file_list:
      fname = item+"/data_plots/wiki_summary.txt"
      try:
          fin = open(fname)
          fout.write(fin.readline())
      except IOError:
          print "Failed for dir", item

def main():
    do_summary("output/2017-03_summary")

if __name__ == "__main__":
    main()
