import glob
import numpy
import cdb

optics_keys = ["beta4d", "emit4d", "alpha4d", "beta_x", "emit_x", "alpha_x", "beta_y", "emit_y", "alpha_y"]

def get_run_comment(run):
    bl = cdb.Beamline()
    bl_list = bl.get_beamline_for_run(run)[long(run)]
    #print bl_list.keys()
    return bl_list['end_notes']

def jacobian_row(plus, minus, delta_keys, magnet_list):
    dvar_dict = {}
    for var in plus.keys():
        try:
            dvar_dict[var] = float(plus[var])-float(minus[var])
        except (ValueError, TypeError):
            pass
    matrix_row = [0. for j in range(len(delta_keys))]
    for i, magnet in enumerate(magnet_list):
        top = None
        if abs(dvar_dict[magnet]) > 0.1:
           top = dvar_dict[magnet]
           print magnet, plus[magnet], minus[magnet], dvar_dict[magnet]
           break
    if top == None:
        print [(key, plus[key]) for key in magnet_list]
        print [(key, minus[key]) for key in magnet_list]
        print [(key, dvar_dict[key]) for key in magnet_list]
    for j, optics in enumerate(delta_keys):
        matrix_row[j] = dvar_dict[optics]/top # dJ/d<twiss>
    return i, matrix_row

def make_ellipse_run_pairs(nominal, list_of_pairs):
    matrix = []
    delta_keys = ["beta4d", "alpha4d", "beta_x", "alpha_x", "beta_y", "alpha_y"]
    magnet_list = ["Q4", "Q5", "Q6", "Q7", "Q8", "Q9"]

    for pair in list_of_pairs:
        file_list = get_folder_list(pair)
        plus = get_run_summary(file_list[0])
        minus = get_run_summary(file_list[1])
        matrix.append(jacobian_row(plus, minus, delta_keys, magnet_list))
    matrix = sorted(matrix)
    matrix = [row[1] for row in matrix]
    print "d<optics>.../dJ"

    print "".rjust(10),
    for key in delta_keys:
        print key.rjust(10),
    print
    for i, row in enumerate(matrix):
        print magnet_list[i].rjust(10),
        for element in row:
            print str(round(element, 3)).rjust(10),
        print
    folder = get_folder_list([nominal])[0]
    nominal = get_run_summary(folder)
    delta_currents = [10., -20., -30., 40., 20., -20.]
    new_optics = [float(nominal[key]) for key in delta_keys]
    for i, magnet_row in enumerate(matrix):
        for j, optics_element in enumerate(magnet_row):
            new_optics[j] += delta_currents[i]*optics_element
    for i, current in enumerate(delta_currents):
        print magnet_list[i], nominal[magnet_list[i]]+current
    print new_optics

def get_run_summary(a_file):
    global optics_keys
    keys = optics_keys
    run_number = a_file.split("-")
    run_number = [item for item in run_number if "96" in item or "97" in item][0]
    run_number = run_number.split("/")[0]
    ell_file = a_file+"/ellipse_summary.txt"
    fin = open(ell_file)
    data = {'run':long(run_number)}
    for line in fin.readlines():
        for key in optics_keys:
            if key in line:
                data[key] = line.split()[1]
        if "tku cut: all" in line:
            break
    bl = cdb.Beamline()
    bl_dict = bl.get_beamline_for_run(data['run'])[data['run']]
    for key in bl_dict:
        if type(bl_dict[key]) != type({}):
            data[key] = bl_dict[key]
    for key in bl_dict["magnets"]:
        data[key] = bl_dict["magnets"][key]["set_current"]
    return data

def make_ellipse_summary_table(file_list):
    global optics_keys
    magnet_list = ["Q4", "Q5", "Q6", "Q7", "Q8", "Q9"]
    keys = ["run"]+optics_keys+magnet_list+["end_notes"]
    print "|",
    for key in keys:
        print key.rjust(10), "|",
    print
    summary = []
    for a_file in file_list:
        data = get_run_summary(a_file)
        for key in keys:
            if key not in data:
                data[key] = "None"
            print  "|", str(data[key]).rjust(10),
        print "|"
        summary.append(data)
    print 
    return summary

  

def make_wiki_summary_table(file_list):
    print "| name                   |  runs      | lmc1234  | tof1     | tof2     | time         | nevents  | us cut   | ds cut   |" 
    for a_file in file_list:
        wiki_file = a_file+"/wiki_summary.txt"
        try:
            fin = open(wiki_file)
            print fin.readline().rstrip("\n")
        except:
            pass

def get_folder_list(run_list):
    file_list = []
    for item in run_list:
        file_list += glob.glob("output/2017-02_1_reco/plots*"+str(item)+"*")#plot_3-140-"+item+"/ellipse*.txt")
    file_list = sorted(file_list)
    return file_list

def main():
    run_list = ["200"]
    file_list = get_folder_list(run_list)
    make_wiki_summary_table(file_list)
    make_ellipse_summary_table(file_list)
    return
    make_ellipse_run_pairs(9701, [[9702, 9703], [9704, 9705], [9706, 9707], [9708, 9710], [9711, 9712], [9713, 9714]])
    make_ellipse_summary_table(get_folder_list([9701]))
    print
    #make_ellipse_run_pairs(9634, [[9639, 9640], [9641, 9642], [9643, 9644], [9645, 9646], [9647, 9648], [9649, 9650]])
    #print
    #make_ellipse_summary_table(get_folder_list([9634, 9655, 9656, 9657]))

if __name__ == "__main__":
    main()
