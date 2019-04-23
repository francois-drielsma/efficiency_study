import os
import shutil
import subprocess
import sys

def geometry_file(run):
    run_str = str(run).rjust(5, "0")
    geometry_dir = "/home/astg/scarf148/data/geometry/geometry_"+run_str
    os.symlink(geometry_dir, "geometry_"+run_str)
    geometry_name = geometry_dir+"/ParentGeometryFile.dat"
    return geometry_name

def scripts():
    return "/home/astg/scarf148/work/amplitude/scripts/"

LOG = None
def log_file():
    global LOG
    if LOG == None:
        LOG = open("check_all_geometry.log", "w")
    return LOG
    
def dir_hack(run):
    plot_path = "geometry_plots/"+str(run).rjust(5, "0")
    if os.path.exists(plot_path):
        shutil.rmtree(plot_path)
    os.makedirs(plot_path+"/plots")
    return plot_path

def check(run):
    do_fields = False
    do_geometry = True
    log = log_file()
    geometry_name = geometry_file(run)
    print "Checking", run
    sys.stdout.flush()

    if do_fields:
        pid = subprocess.Popen(["python", scripts()+"bin/make_field_map.py",
                                "--simulation_geometry_filename", geometry_name],
                                stdout=log, stderr=subprocess.STDOUT)
        pid.wait()
        if pid.returncode != 0:
            raise RuntimeError("Failed on fields with return code "+str(pid.returncode))
        print "    ...done fields"
    if do_geometry:
        pid = subprocess.Popen(["python", scripts()+"bin/check_geometry.py",
                                "--simulation_geometry_filename", geometry_name],
                                stdout=log, stderr=subprocess.STDOUT)
        pid.wait()
        if pid.returncode != 0:
            raise RuntimeError("Failed on geometry with return code "+str(pid.returncode))
        print "    ...done geometry"

def main():
    log_file()
    run_list = [9962, 9966, 9970, 9971]+[10051, 10052, 10064, 10069]+\
               range(10444, 10448)+range(10483, 10487)
    run_list = [7469]
    for run in run_list:
        here = os.getcwd()
        os.chdir(dir_hack(run))
        try:
            check(run)
        except Exception:
            sys.excepthook(*sys.exc_info())
        finally:
            os.chdir(here)

if __name__ == "__main__":
    main()
 