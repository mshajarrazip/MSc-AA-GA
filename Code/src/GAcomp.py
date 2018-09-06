"""Genetic Algorithm components.



"""


from random import choice, randint


from general import *
from problem import info
#from AAcomp import AA_generator

def best_chr(pool):
    """Return the cheapest chromosome in pool.

    Parameter(s):
        pool:list(chromosome())
            - A pool of chromosomes.
    """
    min = pool[0]
    for i in range(len(pool)):
        if pool[i].fitness < min.fitness:
            min = pool[i]
    return min

def worst_chr(pool):
    """Return the most expensive chromosome in pool.

        Parameter(s):
            pool:list(chromosome())
                - A pool of chromosomes.
        """
    max = pool[0]
    for i in range(len(pool)):
        if pool[i].fitness > max.fitness:
            max = pool[i]
    return max

class GApool(list):

    def ismember(self, chr):
        """Checks if chr is already in pool.

        """
        stat = False
        pool = self
        for n in range(len(pool)):
            if (chr.genotype == self[n].genotype):
                stat = True
                break
        return stat

class chromosome(object):
    """Chromosome class.

    When an instance of chromosome is created, its info record is filled first.
    To set its genotype, assign a 0-1 vector to genotype.setter or assign a list of solution
    columns to S.setter. fitness and feasibility status will be set automatically.

    Attribute(s):
        info:problem.info
            - Problem instance info
        __genotype:list({0, 1})
            - Individual genotype.
        fitness:float
            - Individual fitness.
        __S:list(int)
            - List of solution columns
        w:list(int)
            - For each element, the number of columns that covers it.
        status:boolean
            - Feasibility status. True if feasible, False otherwise.
        HD:dict
            - A table of hamming distances to other strings. Keyword is the other string's id.
        FD:dict
            - A table of fitness distances to other strings. Keyword is the other string's id.

    """

    def __init__(self, prob):
        self.info = info()
        self.info.set(prob)
        self.HD = dict()
        self.FD = dict()
        self.fitness = 0

    @property
    def genotype(self):
        return self.__genotype

    @genotype.setter
    def genotype(self, g2):
        g = map(int, g2)
        self.__genotype = g
        self.__S = []
        for x in range(len(g)):
            if g[x] == 1:
                self.__S.append(x)
        self.setw()
        self.setfitness()



    @property
    def S(self):
        return self.__S

    @S.setter
    def S(self, s):
        self.__S = s
        self.__genotype = [0]*self.info.n
        for j in self.__S:
            self.genotype[j] = 1
        self.setw()
        self.setfitness()


    def printinfo(self):
        pass
        '''
        disp("*" * 20)
        disp("Chromosome updated: ")
        disp("ID:", self.ID)
        disp("Solution: ", self.S)
        disp("genotype: ", self.genotype)
        disp("Value: ", self.fitness)
        disp("Feasible: ", self.status)
        disp("*" * 20)
        '''
    def setw(self):
        S = self.S
        Ij = self.info.Ij
        w = [0 for i in self.info.I]
        for j in S:
            for i in Ij[j]:
                w[i] += 1
        self.w = w
        self.status = all(w[i]>0 for i in self.info.I)
        if (self.status == False):
            self.makefeasible()

    def setfitness(self):
        S = self.S
        cost = self.info.C
        val = 0
        for j in S:
            val += cost[j]
        self.fitness = val

    def makefeasible(self):
        """Feasibility operator.

        Converts infeasible solutions to feasible ones.

        """
        Ji = self.info.Ji
        Ij = self.info.Ij
        C = self.info.C
        w = list(self.w)
        S = list(self.S)
        U = [i for i in range(len(w)) if w[i] == 0] # set of uncovered rows

        # adding columns to cover uncovered rows
        while len(U) > 0:
            i = U.pop(0)
            # Finding the first column j in Ji that minimizes Cj / |U X Ij| (Cj/number of uncovered rows which j covers)
            UxIj = [] # |U X Ij|
            for j in Ji[i]:
                count_i = 0
                for i2 in Ij[j]:
                    if w[i2] == 0:
                        count_i += 1
                if (count_i == 0):
                    count_i = 0.5
                UxIj.append(count_i)
            Q = [float(C[j])/d for j, d in zip(Ji[i], UxIj)]
            Qj = Q.index(min(Q))
            j = Ji[i][Qj]
            S.append(j)
            #update w
            for i in Ij[j]:
                w[i] += 1

        S = sorted(S)
        # removing redundant columns
        T = list(S)
        while (len(T) > 0):
            j = choice(T)
            T.remove(j)
            stat = True
            for i in (Ij[j]):
                stat = stat and (w[i] > 1)
            if (stat == True):
                S.remove(j)
                for i in (Ij[j]):
                    w[i] = w[i] - 1

        self.S = S

    def setID(self, gen, chrno = 999):
        """ID for each chromosome is set by:
        <Generation at which chromosome was generated>_<fitness value>_<time stamp>

        """
        self.ID = "g{0}ge_f{1}fe_{2}".format(str(gen),str(self.fitness), str(chrno))

    def setHD(self, other):
        """Sets the hamming distance between self and other.

        Parameter(s):
            other:chromosome
                - Another chromosome to compare to
        """

        kw = other.ID
        ostring = other.genotype
        HD = sum([1 for st, so in zip(self.genotype, other.genotype) if (st != so)])
        self.HD[kw] = HD

    def setFD(self, other):
        """Sets the fitness distance between self and other.

                Parameter(s):
                    other:chromosome
                        - Another chromosome to compare to
                """
        kw = other.ID
        FD = abs(self.fitness - other.fitness)
        self.FD[kw] = FD

class population(object):
    """Population class.

    Attribute(s):
        pool:list(:chromosome)
            - The population pool
        BFS:chromosome
            - The Best Feasible Solution
    """

    def __init__(self):
        self.pool = GApool()
        self.chridno = 0

    def initialize(self, prob, N = 100, init_method = "Random", AA_param = (), PopSource = "", AutoN = False, N_seed = 0):
        """Fill pool with N chromosomes generated by the method specified in init_method.

        Parameters:
             prob:problem
                - The problem instance.
             N:int (Default = 100)
                - Size of population.
             init_method:str (Default = "Random")
                - Population generation method.
             AA_param:(Mode:str, factor:float, Randomizer:int, (Obj:list(int), A:list([{0-1}])
                - Parameters to be passed to the solver:
                    Mode: Mode of solution generation (AAR/ AAD) (Default = "AAD")
                    factor: Range of randomization. (Default = 0)
                    Randomizer: Mode of randomization. (1:Pseudorandom, 2:Sobol etc) (Default = 1)
            PopSource:str
                - A readily available population. Data file format is as follows:
                    population size, and for each solution, the solution value and the number of columns included in
                    the solution.
            Max_N: bool (Default = False)
                - If default is True, then the population pool, if not == N, will be populated by random strings
                    until it is reaches size N. Works on LoadPopulation only.
        """

        disp_status("Initializing population")

        pool = self.pool

        BFS = self.BFS = chromosome(prob)
        if init_method == "Random":
            while len(pool) < N:
                chr = gen_random(prob)
                #disp("chr {3} generated | S {0} | genotype {1} | fitness {2} | {4}"
                #     .format(len(chr.S), chr.genotype.count(1), chr.fitness, len(pool) + 1, chr.status))
                if pool.ismember(chr) == False:
                    self.chridno += 1
                    chr.setID(0, self.chridno)
                    chr.printinfo()
                    pool.append(chr)
            self.update_BFS()
            self.update_avgfitness()

        elif init_method == "LoadPopulation":
            poppool = []
            with open(PopSource, "r") as f:
                dat = f.read().split()
            popsize = int(dat.pop(0))

            # Filter the individuals in sampling population
            for n in range(popsize):
                val = int(dat.pop(0))
                I_size = int(dat.pop(0))
                S = [int(dat.pop(0)) for j in range(I_size)]
                if S not in poppool:
                    poppool.append(S)
                    chr = chromosome(prob)
                    chr.S = S
                    chr.setID(0)
                    chr.printinfo()
                    pool.append(chr)
            disp_status("Loaded population size: {0}".format(len(pool)))

            if (N_seed > 0):
                seedsize = int((N_seed) * N)
                if (len(pool) < seedsize):
                    print("SEED ERROR: Not enough strings for seeding.")
                    return 2
                while (len(pool) > seedsize):
                    pool.pop(randint(0, len(pool) - 1))
                print "\tPopulation seeded with {0} individuals from LoadPopulation.".format(len(pool))
                while len(
                        pool) < N:  # if loaded population has <N members, random members are generated to fill the spot
                    chr = gen_random(prob)
                    if pool.ismember(chr) == False:
                        chr.setID(0)
                        chr.printinfo()
                        pool.append(chr)
                print "\tFilled to N with random: Population size is now {0}.".format(len(pool))
            else:
                while (len(
                        pool) > N):  # if loaded population has >N members, we remove size(sampling population) - N randomly to make N-szie population
                    # reduce population size by randomly choosing the strings among AA-pop
                    pool.pop(randint(0, len(pool) - 1))
                    pass

            if (len(pool) < N):
                print "ERROR: size of loaded population < Minimum population size ({0})".format(N)
                return 1

            self.update_BFS()
            self.update_avgfitness()

        '''
        elif init_method == "AA":
            if AA_param == ():
                AA_param = ("AAD", 0, 1)
            prob_coefficient = (prob.C, prob.A)
            while len(pool) < N:
                AA = AA_generator(AA_param, prob_coefficient)
                chr = chromosome(prob)
                chr.genotype = AA.get_int_sol()
                #disp("pool size {5} | chr {3} generated | S {0} | genotype {1} | fitness {2} | {4}"
                #     .format(len(chr.S), chr.genotype.count(1), chr.fitness, len(pool) + 1, chr.status, len(pool)))
                if pool.ismember(chr) == False:
                    chr.setID(0)
                    chr.printinfo()
                    pool.append(chr)
            self.update_BFS()
            self.update_avgfitness()
        '''
        disp_status("done")
        return 0

    def insert(self, chr):
        """Insert chr into the population pool by steady state replacement.

        Replaces a randomly selected individual with chr if the individual has a fitness value larger than the average
         and if chr is better than the individual. Does nothing if chr is already in the pool. Does nothing if chr is
         worst than all individuals in the pool.

        Parameter(s):
            chr:chromosome
                - An instance of chromosome.
        """
        pool = self.pool
        '''
        disp_head("POPULATION REPLACEMENT")
        disp("Current pool:")
        for i in range(len(pool)):
            disp(i, pool[i].fitness, pool[i].status)

        if chr.fitness > worst_chr(pool).fitness:
            disp("chr is worst.")
            return

        if self.pool.ismember(chr)==True:
            disp("{0} is already in pool".format(chr.fitness))
            return
        '''

        avg_fitness = self.avg_fitness
        cand = []
        chosen = 0
        for n in range(len(pool)):
            if (pool[n].fitness > avg_fitness):
                cand.append(n)
        if len(cand) == 0:
            chosen = randint(0, len(pool)-1)
        else:
            chosen = choice(cand)
        pool.pop(chosen)
        pool.append(chr)

        '''
        disp("Steady-state replacement: {0} replaced {1}".format(chr.fitness, pool[chosen].fitness))
        disp("New pool:")
        for i in range(len(pool)):
            disp(i, pool[i].fitness, pool[i].status)
        '''
        self.update_BFS()
        self.update_avgfitness()

    def update_BFS(self):
        self.BFS = best_chr(self.pool)

    def update_avgfitness(self):
        total = 0
        pool = self.pool
        for n in pool:
            total += n.fitness
        avg = total / len(pool)
        self.avg_fitness = avg

def gen_random(prob):
    """ Return a randomly generated chromosome.

    Parameter(s):
        prob:problem
            - Problem information, an instance of problem.

    """
    chr = chromosome(prob)
    info = chr.info

    I = info.I
    J = info.J
    Ji = info.Ji
    Ij = info.Ij
    m = len(I)
    n = len(J)

    k = 5

    # randomizing initial solution
    S = []
    w = [0 for i in I]

    for i in I:
        A = Ji[i][:k]
        j = choice(A)
        if (j in S):
            continue
        S.append(j)
        for ij in Ij[j]:
            w[ij] = w[ij] + 1
    S = sorted(S)

    # column reduction
    T = list(S)
    while (len(T) > 0):
        j = choice(T)
        T.remove(j)
        stat = True
        for i in (Ij[j]):
            stat = stat and (w[i] > 1)
        if (stat == True):
            S.remove(j)
            for i in (Ij[j]):
                w[i] = w[i] - 1

    chr.S = S

    return chr

def gen_AAD(prob, AAparam = ()):
    """Return a chromosome generated by AAD.

    Parameter(s):
        AAparam:(factor:int)
    """
    chr = chromosome(prob)
    info = chr.info

    A = prob.A


