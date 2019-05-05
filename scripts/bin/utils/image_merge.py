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


def merge(image_files, directories, plot_dir, will_require_complete, do_gthumb = True):
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
        if will_require_complete and len(image_list) != len(directories):
            continue
        print image_list, merge_target
        subprocess.Popen(["convert"]+image_list+["-append", merge_target])
    if do_gthumb:
        proc = subprocess.Popen(["gthumb"]+[plot_dir]) #gthumb_list)
        proc.wait()

def main():
    image_files = [
        "*",
    ]

    recon_directories = [
        "iteration_6_high-stats/2016-04_1.2_reco/plots_3-140_v2.8.5/",
        "iteration_6_high-stats/2016-04_1.2_reco/plots_6-140/",
        "iteration_6_high-stats/2016-04_1.2_reco/plots_10-140/",
    ]

    mc_directories = [
        "iteration_6_high-stats/2016-04_1.2_mc/plots_3-140_MC/",
        "iteration_6_high-stats/2016-04_1.2_mc/plots_6-140_MC/",
        "iteration_6_high-stats/2016-04_1.2_mc/plots_10-140_MC/",
    ]

    e3_directories = [
        "iteration_6_high-stats/2016-04_1.2_reco/plots_3-140/",
        "iteration_6_high-stats/2016-04_1.2_mc/plots_3-140_MC/",
    ]

    e6_directories = [
        "iteration_6_high-stats/2016-04_1.2_reco/plots_6-140/",
        "iteration_6_high-stats/2016-04_1.2_mc/plots_6-140_MC/",
    ]

    e10_directories = [
        "iteration_6_high-stats/2016-04_1.2_reco/plots_10-140/",
        "iteration_6_high-stats/2016-04_1.2_mc/plots_10-140_MC/",
    ]

    merge(image_files, mc_directories, "plots/mc", True, False)
    merge(image_files, recon_directories, "plots/recon", True, False)

if __name__ == "__main__":
    main()

