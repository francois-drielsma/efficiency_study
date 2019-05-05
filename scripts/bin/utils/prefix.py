import os
import glob

input_file_list = glob.glob("/work/ast/cr67/10051_systematics_v3/*/?")+\
                  glob.glob("/work/ast/cr67/10051_systematics_v3/*/??")+\
                  glob.glob("/work/ast/cr67/10051_systematics_v3/*/???")
input_file_list = sorted(input_file_list)

for src_file in input_file_list:
    target_file = list(os.path.split(src_file))
    target_file[1] = target_file[1].rjust(4, '0')
    target_file = os.path.join(target_file[0], target_file[1])
    os.rename(src_file, target_file)
    print src_file, target_file
