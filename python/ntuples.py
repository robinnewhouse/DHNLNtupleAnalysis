from array import array
import ROOT

# key: https://docs.python.org/2/library/array.html
# val: https://root.cern.ch/doc/master/classTTree.html
data_type = {'c': 'C', 'b': 'B', 'B': 'b', 'u': 'python unicode', 'h': 'I',
             'H': 'i', 'i': 'I', 'I': 'i', 'l': 'L', 'L': 'l', 'f': 'F', 'd': 'D',
             }

class Ntuples:
    def __init__(self, tree_name='ntuple', *args, **kwargs):
        """
        This is a simple class to handle the creation and filling of micro-ntuples
        to give more flexibility in plotting later on.
        Trees will be added automicticall when using the [] operator.
        e.g. micro_ntuples['secVtx_VSI_ntrk'] = n_tracks
        will create the tree and fill the associated array with the n_tracks value.
        micro_ntuples.fill() will write all current values to each root tree.
        """
        self.arrays = {}
        self.ttree = ROOT.TTree(tree_name, 'Micro ntuples')

    def fill(self):
        """Fills root ntuple"""
        self.ttree.Fill()

    def write(self):
        """Writes root ntuple to currently open file"""
        self.ttree.Write()

    def reset(self):
        """Resets ntuple values to non-physical numbers (i.e. not zero)"""
        # TODO test if this is necessary
        # TODO implement this function
        pass

    def add(self, name, dtype='d'):
        """
        Creates new ntuple branch with the name `name` and the data type `dtype`.
        If dtype is not specified it defaults to a signed double.
        Python and root data types are associated using the data_type dictionary.
        :param name: Name used to store ntuple in output tree.
        :param dtype: Python array data type. e.g. d (signed double), i (signed integer).
        """
        self.arrays[name] = array(dtype, [0])
        self.ttree.Branch(name, self.arrays[name], "{}/{}".format(name, data_type[dtype]))

    def __setitem__(self, key, value):
        """
        Primary setter function accessed using [] operator. Will create tree if one does not exist.
        :param key: Ntuple name
        :param value: Associated value.
        :return:
        """
        if key not in self.arrays:  # if exists, update array value
            self.add(key, 'd')  # if doesn't exist
        self.arrays[key][0] = value


if __name__ == '__main__':
    """Run this file individually for testing"""
    file = ROOT.TFile.Open('testing.root', 'RECREATE')
    ntuples = Ntuples('testing')
    ntuples.add(name='ints', dtype='i')  # Specify signed int
    ntuples.add(name='nothings', dtype='i')  # this will be filled with zeros
    for i in range(9):
        ntuples['some_doubles'] = i  # don't care about the data type. it will be double
        ntuples['double_doubles'] = i * 2
        ntuples['some_ints'] = -i * 3
        ntuples.fill()
    ntuples.write()
    file.Close()

    # open and read to see what's filled
    import uproot

    uproot.open('testing.root')['testing'].arrays()
