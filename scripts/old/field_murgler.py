import maus_cpp.globals

class Murgler(object): # pylint: disable=R0902
    """Plotting class handles generation of particles and tracking"""
    def __init__(self, verbose):
        """
        Initialise the tracking class
        - configuration should be a valid datacard set in string representation
        """
        self.modules = maus_cpp.globals.get_monte_carlo_mice_modules()
        self.verbose = verbose

    def update(self):
        maus_cpp.globals.set_monte_carlo_mice_modules(self.modules)

    def set_scale_factor(self, mod_name, scale_factor):
        if not self._recursive_set_scale_factor(mod_name, self.modules, scale_factor):
            raise KeyError("Failed to find module "+str(mod_name))

    def misalign_mice_module(self, dx, dxp, dy, dyp, dz, mod_name):
        if not self._recursive_misalign_mice_module(dx, dxp, dy, dyp, dz, mod_name, self.modules):
            raise KeyError("Failed to find module "+str(mod_name))

    def get_scale_factor(self, mod_name):
        return self._recursive_get_scale_factor(mod_name, self.modules)

    def _recursive_misalign_mice_module(self, dx, dxp, dy, dyp, dz, mod_name, mice_mod):
        if mice_mod.get_name() == mod_name:
            position = mice_mod.get_property("Position", "hep3vector")
            position["x"] = dx
            position["y"] = dy
            position["z"] = dz
            mice_mod.set_property("Position", "hep3vector", position)
            rotation = mice_mod.get_property("Rotation", "hep3vector")
            rotation["x"] = dxp
            rotation["y"] = dyp
            mice_mod.set_property("Rotation", "hep3vector", rotation)
            if self.verbose:
                print "Setting", mod_name, "Position to", position, "and Rotation to", rotation
            return True
        children = mice_mod.get_children()
        for mod in children:
            if self._recursive_misalign_mice_module(dx, dxp, dy, dyp, dz, mod_name, mod):
                mice_mod.set_children(children)
                return True
        return False

    def _recursive_get_scale_factor(self, mod_name, mice_mod):
        if mice_mod.get_name() == mod_name:
            try:
                return mice_mod.get_property("ScaleFactor", "double")
            except KeyError:
                return 1.
        else:
            for mod in mice_mod.get_children():
                scale = self._recursive_get_scale_factor(mod_name, mod)
                if scale != None:
                    return scale
            return None                  

    def _recursive_set_scale_factor(self, mod_name, mice_mod, scale_factor):
        if mice_mod.get_name() == mod_name:
            pol = mice_mod.get_property("ScaleFactor", "double")
            mice_mod.set_property("ScaleFactor", "double", scale_factor)
            if self.verbose:
                print "Setting", mod_name, "ScaleFactor to", scale_factor
            return True
        children = mice_mod.get_children()
        for mod in children:
            if self._recursive_set_scale_factor(mod_name, mod, scale_factor):
                mice_mod.set_children(children)
                return True
        return False


