import shutil
import glob
import os
import subprocess

def get_image_name(path_to_image):
    return os.path.split(path_to_image)[1][:-4]

def target_list(image_files, directories):
    target_list = []
    for image in image_files:
        for directory in directories:
            target = os.path.join(directory, image)+".png"
            target_list += glob.glob(target)
    image_files = [os.path.split(target)[1] for target in target_list]
    image_files = list(set(image_files))
    target_list = {}
    print image_files
    for image in image_files:
        file_names = [os.path.join(directory, image) for directory in directories]
        file_names = [fname for fname in file_names if os.path.exists(fname)]
        target_list[image] = file_names
    return target_list


def merge(image_files, directories, plot_dir):
    try:
        shutil.rmtree(plot_dir)
    except OSError:
        pass
    os.makedirs(plot_dir)
    if not os.path.isdir(plot_dir):
        raise OSError(str(plot_dir)+" does not exist")
    my_target_list = target_list(image_files, directories)
    for image, image_list in my_target_list.iteritems():
        image_name = os.path.splitext(image)[0]
        merge_target = plot_dir+"/"+image_name+"_merge.png"
        print image_list, merge_target
        subprocess.Popen(["convert"]+image_list+["-append", merge_target])

    proc = subprocess.Popen(["gthumb"]+[plot_dir]) #gthumb_list)
    proc.wait()

def main():
    image_files = [
        "*",
    ]

    recon_directories = [
        "output/2016-04_1.2_reco/plots_3-140/",
        "output/2016-04_1.2_reco/plots_6-140/",
        "output/2016-04_1.2_reco/plots_10-140/",
    ]

    mc_directories = [
        "output/2016-04_1.2_mc/plots_3-140_MC/",
        "output/2016-04_1.2_mc/plots_6-140_MC/",
        "output/2016-04_1.2_mc/plots_10-140_MC/",
    ]

    e3_directories = [
        "output/2016-04_1.2_reco/plots_3-140+M3-Test2_MAUS-v2.8.2/",
        "output/2016-04_1.2_mc/plots_3-140_MC_Scale_D1=1.??_D2=0.93*/",
    ]

    e6_directories = [
        "output/2016-04_1.2_reco/plots_6-140+M3-Test2_MAUS-v2.8.5_Full_hi-stats/",
        "output/2016-04_1.2_mc/plots_3-140_MC_Scale_D1=1.02_D2=1.02_DS=1.00_Br2.0_W1.4/",
        "output/2016-04_1.2_mc/plots_6-140_MC_Scale_D1=1.02_D2=1.02_DS=1.00_Hacked-Geom"
   ]

    e10_directories = [
        "output/2016-04_1.2_reco/plots_10-140+M3-Test3_MAUS-v2.8.5_Full_hi-stats/",
        "output/2016-04_1.2_mc/plots_10-140_MC_Scale_D1=1.02_D2=1.02_DS=1.00_Br12.87_W8.4/",
    ]

    merge(image_files, mc_directories, "plots/mc")

if __name__ == "__main__":
    main()

