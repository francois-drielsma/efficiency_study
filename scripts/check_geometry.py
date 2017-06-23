import Configuration
import maus_cpp.globals
import maus_cpp.material

def initialise_maus():
    configuration = Configuration.Configuration().\
                                          getConfigJSON(command_line_args=True)
    maus_cpp.globals.birth(configuration)

def print_materials():
    for x in [120, 101, 99, 0]:
        material = None
        for z in range(1300000, 1400000, 1):
            z /= 100.
            maus_cpp.material.set_position(x, 0., z)
            material_data = maus_cpp.material.get_material_data()
            new_material = material_data['name']
            if new_material != material:
                material = new_material
                print str(x).ljust(10), str(z).ljust(10), material.ljust(20), material_data["volume_name"]
        print

def main():
    initialise_maus()
    print_materials()

if __name__ == "__main__":
    main()
