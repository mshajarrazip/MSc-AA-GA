"""Contain classes and methods to manage problem information.

Class(es):
    info - Record of problem instance information.
    problem - Subclass of info with an additional methods to load the problem.
Method(s):
    parseSCP(source) - Parse ORLib SCP data file and return a dictionary.
"""
from os.path import split


from general import *

class info():
    """Record for information on problem instance.

        Attribute(s):
            name - Name of the problem.
            density - Density of the matrix.
            m - Size of row.
            n - Size of columns
            I - The set of all rows
            J - The set of all columns
            C - Cost of each column j=0,..,n-1
            Ji - The set of columns that cover row i
            Ij - The set of rows covered by column j

        Method(s):
            isSet() - Return False if no field are set.
            set(prob) - Set info fields with prob fields.

        Subclass(es):
            problem
                + A, density
                + parseSCP()

    """
    def __init__(self):
        self.name = -1
        self.m = 0 # Size of row
        self.n = 0 # Size of col
        self.I = []
        self.J = []
        self.C = []
        self.Ji = []
        self.Ij = []

    def isSet(self):
        return (self.name != -1)

    def set(self, prob):
        self.name = prob.name
        self.m = prob.m
        self.n = prob.n
        self.I = prob.I
        self.J = prob.J
        self.Ij = prob.Ij
        self.Ji = prob.Ji
        self.C = map(float,prob.C)


class problem(info):
    """Store problem information.

    Super: info()
    Attribute(s)
        A:list([])
            - The 0-1 matrix Aij

    Method(s):
        load(source) - Parse the data from source file and initialize attributes of problem.

    """

    def __init__(self, source = ''):
        """Load problem from source.

        Parameter(s):
            source:str (Default = '')
                - Path of data file.
        """
        if source != '':
            self.load(source)


    def load(self, source):
        """Loads problem data from source.

        Initialize all variables related to the problem.
        parseSCP(source) stores problem data info in a dictionary and return it to data.
        This function calculates the other values for attributes to this class based on data.

        """
        disp_status("Loading problem data file: {0}".format(source))
        data = parseSCP(source)
        row = data['row']
        col = data['col']
        C = map(float,data['C'])
        Ji = data['Ji']
        I = [i for i in range(data['row'])]  # the set of all rows
        J = [j for j in range(data['col'])]  # the set of all columns
        Ij = [[] for j in J]  # the set or rows covered by column j
        A = [[0]*col for i in I] # the Aij matrix

        # Transposing Ji
        for j in J:
            cov = 0
            for i in Ji:
                if j in i:
                    Ij[j].append(cov)
                cov += 1

        # generating matrices
        for i in range(len(A)):
            for j in Ji[i]:
                A[i][j] = 1

        self.name = data['name']
        self.density = data['density']
        self.m = row
        self.n = col
        self.C = C
        self.Ji = Ji
        self.I = I
        self.J = J
        self.Ij = Ij
        self.A = A

        disp_status("done")

def parseSCP(source):
    """Return a dictionary containing problem information.

    Parameter(s):
        source
            - The directory of SCP data file.
            The format of data file should be based on ORLib's SCP data file for problem sets 4-6:
            Size of rows followed by size of column, and for each column, its costs in increasing order, and for each
            row, the number of columns that cover that row, followed by indices of the column.
            The row and column index is within the range [1,n]. This method converts the range to [0, n).

    Attribute(s):
    //Attributes are accessed directly
        prob_info
            - A dictionary on the problem information:
                'row': size of row
                'col': size of column
                'C': [ Cost of column j ] for j = 1,..,n
                'Ji': [ [ index of column j ] ] for i = 1,..,m and j = [ column j that covers row i ]
                'density': density of ones in the Matrix
    """

    disp('Reading', source)
    prob_info = dict()
    name = split(source)[1]
    prob_info['name'] = name[:name.index(".txt")]

    with open(source) as f:
        data = f.read().split()

    prob_info['row'] = int(data.pop(0))
    m = prob_info['row']
    prob_info['col'] = int(data.pop(0))
    n = prob_info['col']

    cost = []
    for j in range(n):
        cost.append(int(data.pop(0)))
    prob_info['C'] = cost

    A = []
    for i in range(m):
        Ji = int(data.pop(0))
        J = [int(data.pop(0))-1 for j in range(Ji)]
        A.append(J)
    prob_info['Ji'] = A

    q = sum(len(i) for i in A)
    mn = m*n
    prob_info['density'] = float(q)/mn

    disp_dict(prob_info)

    return prob_info






