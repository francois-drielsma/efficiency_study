import shutil
import os

source = "/home/astg/scarf148/work/amplitude/output/2017-02-drielsma"
target = "/home/astg/scarf148/work/systematics_mc/beams"
analysis_path = "data_recorder/"
file_name = "tku_5.json"
analysis_to_run = {
    "plots_2017-2.7_10-140_lH2_empty":10052,
    "plots_2017-2.7_10-140_lH2_full":9970,
    "plots_2017-2.7_10-140_LiH":10486,
    "plots_2017-2.7_10-140_None":10447,
    "plots_2017-2.7_3-140_lH2_empty":10069,
    "plots_2017-2.7_3-140_lH2_full":9971,
    "plots_2017-2.7_3-140_LiH":10483,
    "plots_2017-2.7_3-140_None":10444,
    "plots_2017-2.7_4-140_lH2_empty":10064,
    "plots_2017-2.7_4-140_lH2_full":9962,
    "plots_2017-2.7_4-140_LiH":10484,
    "plots_2017-2.7_4-140_None":10445,
    "plots_2017-2.7_6-140_lH2_empty":10051,
    "plots_2017-2.7_6-140_lH2_full":9966,
    "plots_2017-2.7_6-140_LiH":10485,
    "plots_2017-2.7_6-140_None":10446,
}

def main():
    for analysis_name, run_number in analysis_to_run.iteritems():
        src_path = os.path.join(source, analysis_name)
        src_path = os.path.join(src_path, analysis_path)
        src_path = os.path.join(src_path, file_name)
        if not os.path.exists(src_path):
            print "Couldn't find", src_path
            continue
        tgt_path = os.path.join(target, str(run_number))
        try:
            os.makedirs(tgt_path)
        except OSError:
            pass
        tgt_path = os.path.join(tgt_path, file_name)
        print "Copying", tgt_path
        print "     to", src_path
        try:
            shutil.rmtree(tgt_path)
        except OSError:
            if os.path.exists(tgt_path):
               print "Failed to clear", tgt_path
               continue
        shutil.copy(src_path, tgt_path)
        
if __name__ == "__main__":
    main()