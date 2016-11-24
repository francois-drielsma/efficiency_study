import cdb

def get_beamline(run_number):
    bl = cdb.Beamline()
    beamline_all = bl.get_beamline_for_run(run_number)[run_number]
    return beamline_all

def get_tof_triggers_sum(tof_key, run_list):
    trig_list = [get_beamline(run_number)["scalars"][tof_key] for run_number in run_list]
    return sum(trig_list)

def get_time_sum(run_list):
    start_time = [get_beamline(run_number)["start_time"] for run_number in run_list]
    end_time = [get_beamline(run_number)["end_time"] for run_number in run_list]
    delta_time = [end_time[i] - start for i, start in enumerate(start_time)]
    total_seconds = sum([delta.total_seconds() for delta in delta_time])
    minutes = int(total_seconds)/60
    hours = minutes/60
    minutes = minutes % 60
    return hours, minutes

def parse_one_setting(run_list):
    tof1_triggers = get_tof_triggers_sum("ToF1 Triggers", run_list)
    tof2_triggers = get_tof_triggers_sum("ToF2 Triggers", run_list)
    lmc1234 = get_tof_triggers_sum("LMC-1234 Count", run_list)
    hours, minutes = get_time_sum(run_list)
    return {"runs":run_list, "tof1_triggers":tof1_triggers, "tof2_triggers":tof2_triggers, "lmc1234":lmc1234, "time":[hours, minutes]}

# +100pc 155 13374/25592 0.523
# +ao    118 11575/20535 0.564

if __name__ == "__main__":
    print parse_one_setting([8526] )
