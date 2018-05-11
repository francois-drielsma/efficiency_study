import os
import shutil
import subprocess
import sys

def geometry_file(run):
    run_str = str(run).rjust(5, "0")
    geometry_dir = "/home/cr67/work/reco/geometry/geometry_"+run_str
    os.symlink(geometry_dir, "geometry_"+run_str)
    geometry_name = geometry_dir+"/ParentGeometryFile.dat"
    return geometry_name

def scripts():
    return "/home/cr67/work/2016-11-18_emittance-analysis/amplitude/scripts/bin/"

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
    os.makedirs(plot_path)
    return plot_path

def check(run):
    log = log_file()
    geometry_name = geometry_file(run)
    print "Checking", run
    sys.stdout.flush()

    pid = subprocess.Popen(["python", scripts()+"make_field_map.py",
                            "--simulation_geometry_filename", geometry_name],
                            stdout=log, stderr=subprocess.STDOUT)
    pid.wait()
    if pid.returncode != 0:
        raise RuntimeError("Failed on fields with return code "+str(pid.returncode))
    print "    ...done fields"
    pid = subprocess.Popen(["python", scripts()+"check_geometry.py",
                            "--simulation_geometry_filename", geometry_name],
                            stdout=log, stderr=subprocess.STDOUT)
    pid.wait()
    if pid.returncode != 0:
        raise RuntimeError("Failed on geometry with return code "+str(pid.returncode))
    print "    ...done geometry"

def main():
    log_file()
    for run in [9962, 10064, 10445, 10484]: #range(10444, 10448)+range(10483, 10487):
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
 