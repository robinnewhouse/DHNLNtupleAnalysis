from array import array
import ROOT

# key: https://docs.python.org/2/library/array.html
# val: https://root.cern.ch/doc/master/classTTree.html
data_type = {
    'c': 'C',
    'b': 'B',
    'B': 'b',
    'u': 'python unicode',
    'h': 'I',
    'H': 'i',
    'i': 'I',
    'I': 'i',
    'l': 'L',
    'L': 'l',
    'f': 'F',
    'd': 'D',
}

'''
This is a simple class to handle the creation and filling of micro-ntuples
to give more flexibility in plotting later on
'''
class Ntuples():
    def __init__(self, tree_name='ntuple', *args, **kwargs):
        self.arrays = {}
        self.ttree = ROOT.TTree(tree_name, 'Micro ntuples')

    def fill(self):
        self.ttree.Fill()

    def write(self):
        self.ttree.Write()

    def add(self, name, dtype='d'):
        self.arrays[name] = array(dtype, [0])
        self.ttree.Branch(name, self.arrays[name], "{}/{}".format(name, data_type[dtype]))

    def __setitem__(self, key, value):
        if key not in self.arrays:  # if exists, update array value
            self.add(key, 'd')  # if doesn't exist
        self.arrays[key][0] = value


'''
Run this file individually for testing
'''
if __name__ == '__main__':
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
