from random import random
from math import log


import cplex


from general import *

class AA_generator():
    """Solve LP relaxation and perform AA to generate chromosomes.

    """

    def __init__(self, AAparam = (), Prob = (0, 0)):
        """

        Parameter(s):
            AAparam:(Mode:str, factor:float, Randomizer:int, (Obj:list(int), A:list([{0-1}])
                - Parameters to be passed to the solver:
                    Mode: Mode of solution generation (AAR/ AAD) (Default = "AAD")
                    factor: Range of randomization. (Default = 0)
                    Randomizer: Mode of randomization. (1:Pseudorandom, 2:Sobol etc) (Default = 1)
            Prob:(Obj, A)
                - Coefficient of LP problem. Obj - Obj function coefficients, A - The matrix.
        """

        if AAparam == ():
            AAparam = ("AAD", 0, 1)

        Mode = self.mode = AAparam[0]
        factor = AAparam[1]
        Randomizer = self.randomizer = AAparam[2]

        C = Prob[0]
        A = Prob[1]

        if factor > 0:
            factor = random()*factor

        self.factor = factor

        m = len(A)
        n = len(C)
        self.my_real_obj = [c for c in C]
        self.my_obj = [c + factor for c in C]
        self.my_ub = [1.0 for j in range(n)]
        self.my_lb = [0.0 for j in range(n)]
        self.my_ctype = ""
        self.my_colnames = ['x' + str(j) for j in range(n)]
        self.my_rhs = [1.0 for i in range(m)]
        self.my_rownames = ['E' + str(i) for i in range(m)]
        self.my_sense = 'G' * m
        self.my_matrix = A
        prob = cplex.Cplex()
        self.populate_problem(prob)

        disp_head("Solving {0} {1}".format(Mode, factor))
        prob.solve()
        LP_sol = prob.solution.get_values()

        if Mode=="AAD":
            self.solve_AAD(LP_sol)
        elif Mode == "AAR":
            self.solve_AAR(LP_sol)


    def solve_AAD(self, LP_sol):

        x = LP_sol
        X = []
        f = max(v.count(1) for v in self.my_matrix)  # f is the max number of sets in which any elements appear
        for i in range(len(x)):
            if x[i] > (1 / f):
                X.append(1)
            else:
                X.append(0)
        self.int_sol = X

    def solve_AAR(self, LP_sol):
        x = LP_sol
        X = []
        for i in range(len(x)):
            if (self.random_decision(x[i], len(self.my_rownames), 2) == True):
                X.append(1)
            else:
                X.append(0)
        self.int_sol = X

    def get_int_sol(self):
        return self.int_sol

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

    def random_decision(self, Pr, n, c): # Pr is xj, n is len(rows) and c is some constant
        '''returns true if 1 should be included'''
        decision = False
        event = (int) (c* log(n))
        heads = 0
        if (Pr>0):
            heads = 0
            for t in range(event):
                if (random.random() < Pr):
                    heads += 1
        if (heads >=1):
            decision = True
        return decision