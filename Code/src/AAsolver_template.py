''' This module containts functions written to solve the SCP problem.

Hence, the tendency to decribe in the context of SCP.

>> solv = solver(prob, '''

from __future__ import print_function


import random
import math

import cplex


class solver():

    def __init__(self, prob, data, mode = 'LP', factor = 0, fullprint = True, save=False, savelp = False, popsave = ""):
        '''The constructor organizes the values to be used in the
        problem.

        data: dictionary containing the fields related to the problem:
            name: name of problem (string)
            row:  row size, m (int)
            col:  column sizen, n (int)
            cost: cost of each column / coefficient of obj function (int)
            J:    dictionary of i:Ji ( i = 0,.., m, Ji is a list of int)
                  i refers to index of rows, and Ji the set of columns that
                  cover row i

        mode: String referring to the method used to derive the solutions.
            'LP':  LP Relaxation
            'MIP': MIP solution (Branch and Bound, exact)
            'AAD': LP Relaxation + Deterministic Rounding AA
            'AAR': LP Relaxation + Randomized Rounding AA

        factor: Real value for range of randomization at obj function value.
                Initially 0.

        Initializes the field in the solution_info, which is a dictionary containing the solution information

        If fullprint is True, solution written to solution file will be in full length.

        If save is true, solution file will be written.

        If savelp is True, lp info will be written.

        popsave is the name of population file. If save is True, we will check if the solution has already existed in
        popsave. If solution has already been found, then solution will not be saved. By default, if popsave is "",
        then solution file will automatically be saved. In other words, this parameter handles duplicates. By default, all duplicates are printed.
            '''


        self.savelp = savelp
        self.mode = mode
        self.name = data["name"]

        if (factor>0):
            factor = random.random() * factor

        self.solution_info = solution(data["name"], mode, factor)
        self.solution_info.update({
            "density": data["density"]
        })

        m = data["row"]
        n = data["col"]
        self.my_real_obj = [c for c in data["C"]]
        self.my_obj = [c + factor for c in data["C"]]
        self.my_ub = [1.0 for j in range(n)]
        self.my_lb = [0.0 for j in range(n)]
        self.my_ctype = ""
        self.my_colnames = ['x' + str(j) for j in range(n)]
        self.my_rhs = [1.0 for i in range(m)]
        self.my_rownames = ['E' + str(i) for i in range(m)]
        self.my_sense = 'G' * m

        # Populate the matrix

        matrix = []
        for i in range(m):
            var_row = [0 for j in range(n)]
            for j in data['J'][i]:
                var_row[j - 1] = 1
            matrix.append(var_row)
        self.my_matrix = matrix

        # prob = cplex.Cplex()
        prob.set_problem_name("{0} Density: {1} Factor: {2}".format(data["name"], data["density"], factor))
        if (mode == 'LP'):
            self.solution = self.solve(prob, Integral=False)
            self.solution_info.update({
                "OBJ_VAL": prob.solution.get_objective_value(),
            })
        elif (mode == 'MIP'):
            self.solveMIP(prob)
        elif (mode == 'AAD'):
            self.solveAAD(prob)
        elif (mode == 'AAR'):
            self.solveAAR(prob)
        else:
            print("Invalid mode. \n Available solvers: LP/ MIP/ AAD/ AAR")

        if (save == True):
            exist = False
            if (popsave != ""):
                fname = popsave[:popsave.index(".txt")] + "_OBJVAL.txt"
                try:
                    with open(fname, 'r') as f:
                        list = (f.read()).split()
                        print("Previous solutions: ", list)
                        sol = str(self.solution_info.info["OBJ_VAL"])
                        print("current solution: ", sol)
                    if (sol in list):
                        exist = True
                except IOError:
                    pass
            print ('Is Exist? ', exist)
            if (exist == False):
                self.solution_info.write(full=fullprint)

        self.solution_info.display()

    def solveMIP(self, prob):
        n = len(self.my_colnames)
        self.my_ctype = 'I' * n  # set type of variables to integer

        self.solve(prob, Integral = True)

        self.solution_info.update({
            "OBJ_VAL": self.solution_info.info["obj_val"],
            "X": self.solution_info.info["x"],
            "feasible": self.checkfeasible(self.solution_info.info["x"])
        })


    def solveAAD(self, prob):
        '''Applies the Deterministic Rounding Algorithm (AAD) after LP relaxation solution.
                '''

        self.solve(prob)

        x = self.solution_info.info["x"]
        X = []
        f = max(v.count(1) for v in self.my_matrix)  # f is the max number of sets in which any elements appear
        for i in range(len(x)):
            if x[i] > (1 / f):
                X.append(1)
            else:
                X.append(0)
        solution = dict()
        solution["X"] = X

        sum = 0
        for c, xi in zip(self.my_real_obj, solution["X"]):
            sum = sum + c * xi
        solution["OBJ_VAL"] = sum
        solution["feasible"] = self.checkfeasible(solution["X"])

        self.solution_info.update(solution)

    def solveAAR(self, prob):
        '''Applies the Randomized Rounding Algorithm (AAD) after LP relaxation solution.
                '''

        self.solve(prob)


        x = self.solution_info.info["x"]
        X = []
        for i in range(len(x)):
            if (self.random_decision(x[i], len(self.my_rownames), 2) == True):
                X.append(1)
            else:
                X.append(0)

        solution = dict()
        solution["X"] = X

        sum = 0
        for c, xi in zip(self.my_real_obj, solution["X"]):
            sum = sum + c * xi
        solution["OBJ_VAL"] = sum
        solution["feasible"] = self.checkfeasible(solution["X"])

        self.solution_info.update(solution)

    def solve(self, prob, Integral = False):
        ''' Set up the problem information, solve and manage the solutions.
        '''
        if (prob.solution.get_status() == -1):
            self.populate_problem(prob)
            prob.solve()

        self.solution_info.update({
            "status": prob.solution.status[prob.solution.get_status()],
            "obj_val": prob.solution.get_objective_value(),
            "x": prob.solution.get_values()
        })


        if (Integral == False):
            self.solution_info.update({
            "slack": prob.solution.get_linear_slacks(),
            "Pi": prob.solution.get_dual_values(),
            "Dj": prob.solution.get_reduced_costs()
            })


    def populate_problem(self, prob):
        '''Set up the problem information.
        '''

        prob.objective.set_sense(prob.objective.sense.minimize)
        prob.variables.add(obj=self.my_obj,
                           lb=self.my_lb,
                           ub=self.my_ub,
                           types=self.my_ctype,
                           names=self.my_colnames)

        rows = [[[name for name in self.my_colnames], v] for v in self.my_matrix]

        prob.linear_constraints.add(lin_expr=rows,
                                    senses=self.my_sense,
                                    rhs=self.my_rhs,
                                    names=self.my_rownames)

        if (self.savelp == True):
            prob.write("../lp/" + self.name + "_f" + str(self.solution_info.info["factor"]) + ".lp")

    def checkfeasible(self, X):

        matrix = self.my_matrix
        covered = [0 for i in range(len(matrix))]

        for j in range(len(X)):
            if X[j]==1:
                for i in range(len(matrix)):
                    if matrix[i][j] == 1:
                        covered[i] = 1

        covered = covered.count(1)
        if covered == len(matrix):
            return True
        else:
            return False

    def random_decision(self, Pr, n, c): # Pr is xj, n is len(rows) and c is some constant
        '''returns true if 1 should be included'''
        decision = False
        event = (int) (c* math.log(n))
        heads = 0
        if (Pr>0):
            heads = 0
            for t in range(event):
                if (random.random() < Pr):
                    heads += 1
        if (heads >=1):
            decision = True
        return decision
















