import time
import os
import json



class solution():

    special_data = ["X"]

    def __init__(self, name, mode, factor):
        '''Initializes the fields in info, the table that keeps solution information.

        info:
            (Some information about the problem: )
                Name: Name of problem.
                mode: Mode of solution.
                factor: If randomization added to obj function, factor is > 0.
            (Integral solution info:)
                OBJ_VAL: Objective value of integral solution.
                X_len: Size of feasible solution
                X: Integral solution.
                feasible: If solution is feasible, True, else False.
            (CPLEX solution info:)
                status: status of solution (Cplex().solution.status[Cplex().solution.get_status()])
                obj_val: objective value (Cplex().solution.get_objective_value())
                x = prob.solution.get_values()
                Pi = prob.solution.get_dual_values()
                Dj = prob.solution.get_reduced_costs()
                slack = prob.solution.get_linear_slacks()

        '''

        self.info = {
                "Name": name,
                "mode": mode,
                "factor": factor,
                # Integral solution info
                "OBJ_VAL": "N/A", # obj value of integral solution
                "X_len" : "N/A", # size of X
                "X": "N/A", # integral solution
                "feasible": "N/A",
                # CPLEX solution info
                "status": "N/A",
                "obj_val": "N/A",
                "x": "N/A",
                "Pi": "N/A",
                "Dj": "N/A",
                "slack": "N/A",
                # Other problem information
                "density": "N/A",
                "time": "N/A"
            }

    def update(self, newdata):
        ''' Update the information in info.
            newdata is a dictionary containing the infos to be updated.

            if the keyword in newdata matches an entry in self.special_data, special procedures will be taken:
                "X": This field expects the list of columns that forms a set cover. Hence, uses getcollist on
                     self.info["x"] to get the column list and size of the list that corresponds to it.

        '''

        data = newdata
        key = ""
        specialdata = []
        try:
            for kw in newdata:
                key = kw
                if key in self.special_data:
                    specialdata.append(key)
                    continue
                self.info[kw] = data[kw]
        except KeyError:
            print("Field not updated: ",key, " is not a valid key.")

        for kw in specialdata:
            if kw is "X":
                self.info["X"] = self.getcollist(data[kw])
                self.info["X_len"] = len(self.info["X"])

    def update_integral(self):
        pass

    def display(self):
        '''Prints the current solution info solution.
        '''
        sol = self.info

        print("{0:10}: {1}".format("Name", sol["Name"]))
        print("{0:10}: {1}".format("Mode", sol["mode"]))
        print("{0:10}: {1}".format("Factor", sol["factor"]))
        print("SOLUTION INFO (Integral)")
        print(">> {0:10} ({2}): {1}".format("Obj Value", sol["OBJ_VAL"], sol["mode"]))
        print(">>> {0:10}: {1}".format("Solution size", sol["X_len"]))
        print(">>> {0:10}: {1}".format("Solution", sol["X"]))
        print(">>> {0:10}: {1}".format("Feasible", sol["feasible"]))
        print("CPLEX SOLUTION")
        print(">> {0:10}: {1}".format("Status", sol["status"]))
        print(">> {0:10}: {1}".format("Obj Value", sol["obj_val"]))
        print(">>> {0:10}: {1}".format("Solution", sol["x"]))
        print(">>> {0:10}: {1}".format("Slack", sol["slack"]))
        print(">>> {0:10}: {1}".format("Dual Values", sol["Pi"]))
        print(">>> {0:10}: {1}".format("Reduced Costs", sol["Dj"]))
        print("OTHER PROBLEM INFORMATION")
        print(">> {0:10}: {1}".format("Density", sol["density"]))
        print(">> {0:10}: {1}".format("Time", sol["time"]))


    def getcollist(self, x):
        '''Get list of columns represented by integral solution

        '''
        X = []
        for j in range(len(x)):
            if (x[j] > 0):
                X.append(j)
        return X

    def write(self, outdir = "..", filename = "", full = True, dumpjson = True):
        ''' Writes solution value to a file.

        def_rep is the target directory to store the file in. Default is ../results/

        filename is set to <problem name>_<Mode of solution>_f<factor>_<objective value of solution>_<timestamp>.txt
        '''

        sol = self.info
        data = self.info

        #def_rep ="..\\results\\"
        def_rep ="{1}/results/{0}/".format(data["Name"], outdir)

        if (filename==""):
            timestamp = self.timestamp()
            filename = "{0}{1}_{2}_f{3}_v{4}_{5}.txt"\
                .format(def_rep,
                        data["Name"],
                        data["mode"],
                        str(data["factor"]),
                        str(data["OBJ_VAL"]),
                        str(timestamp))

        # creates def_rep if directory is not yet created
        if not os.path.exists(def_rep):
            os.makedirs(def_rep)

        print("Saving solution info ... {0}".format(filename))
        f = open(filename, 'w')

        f.write("Name: {0} \t Mode: {1} \t  Factor: {2} \t Density: {3} \t \n".format(sol["Name"], sol["mode"], sol["factor"], sol["density"]))
        f.write("Time(s): {0} \n".format(sol["time"]))

        # INTEGRAL SOLUTION INFORMATION: Consist of objective value, size of solution, and column indexes that make up the solution
        f.write("<< INTEGRAL >>\n Value: {0} (Feasible? {2})\t Size: {1} \n Set Cover: \n ".format(sol["OBJ_VAL"], sol["X_len"],\
                sol["feasible"]))
        X = sol["X"]
        if X != "N/A":
            linebreak = 0
            for j in range(len(X)):
                f.write(str(X[j]) + "\t")
                linebreak += 1
                if (linebreak == 10):
                    f.write("\n")
                    linebreak = 0
        else:
            f.write(X)
        f.write("\n")

        # CPLEX SOLUTION INFORMATION: Consists of LP/MIP optimal solution value, the values of the solutions,
        # as well as dual values, slack and reduced costs. If full is False, latter three will not be printed.

        f.write("<< CPLEX >>\n Value: {0} ( {1} ) \n".format(sol["obj_val"], sol["status"]))
        if (full == True) & (sol["status"] == "optimal"):
            f.write("Solution: \n")
            f.write("{0} \t {1:30} {2:10} \n".format("Var", "Values", "Reduced Costs"))
            x = sol["x"]
            D = sol["Dj"]
            for j in range(len(x)):
                f.write("x{0} \t {1:30}  {2:10} \n".format(j, str(x[j]), str(D[j])))
            f.write("\n")
            S = sol["slack"]
            P = sol["Pi"]
            f.write("{0} \t {1:30} {2:10} \n".format("Constraint", "Dual", "Slack"))
            for i in range(len(P)):
                f.write("C{0} \t {1:30} {2:30} \n".format(i, str(P[i]), str(S[i])))
        else:
            f.write("Solution: \n")
            x = sol["x"]
            f.write("{0} \t {1:30} \n".format("Var", "Values", "Reduced Costs"))
            for j in range(len(x)):
                if (x[j]>0):
                    f.write("x{0} \t {1:30} \n".format(j, str(x[j])))

        f.close()

        if (dumpjson == True):
            self.writeasjson(outpdir = outdir, filename = filename)

        print("Solution info saved in {0}".format(filename))

    def writeasjson(self, outpdir = "..", filename = ""):
        '''prints the whole of self.info as json to a file. '''
        data = self.info
        def_rep = "{1}/results/{0}/".format(data["Name"], outpdir)
        if not os.path.exists(def_rep):
            os.makedirs(def_rep)
        if (filename==""):
            timestamp = self.timestamp()
            filename = "{0}_{1}_f{3}_v{4}_{5}.txt"\
                .format(def_rep,
                        data["Name"],
                        data["mode"],
                        str(data["factor"]),
                        str(data["OBJ_VAL"]),
                        str(timestamp))

        filename = "{1}/{0}_JSON.txt".format(filename[:filename.index(".txt")], def_rep)

        print("Printing BFS solution info in JSON to: \n \t {0}".format(filename))
        with open(filename, 'w') as outfile:
            json.dump(data, outfile)

    def addtofile(self, outdir = "..", filename=""):

        ''' Appends some solution information to filename.

        ../population stores various solutions in one file.

        Using f, the number of solutions, and for each solution, the value and the column indexes that make up the
        solution is written to a file named <problem name>_<mode of solution>compact_<timestamp>.txt by default.
        filename should be set by user before running multi-tests.

        Using f2, the values of all solutions is gathered in 1 file, named as in f, with an added OBJVAL.
        For ease of recording purposes.
        '''



        if (filename == ""):
            def_rep = "{0}/population/".format(outdir)
            # creates def_rep if directory is not yet created
            if not os.path.exists(def_rep):
                os.makedirs(def_rep)
            timestamp = self.timestamp()
            filename = def_rep + self.info["Name"] + "_" + self.info["mode"] + "compact" + "_" + str(timestamp) + ".txt"

        # append obj value, solution size, and solution to f
        value = self.info["OBJ_VAL"]
        X_len = self.info["X_len"]
        X = self.info["X"]
        Xrow = str(value) + "\t" + str(X_len) + "\n"
        f = open(filename, "a+")
        for i in range(X_len):
            Xrow += str(X[i]) + "\t"
        Xrow += "\n"
        f.write(Xrow)
        f.close()

        # append solution value to f2
        filename = filename[:filename.index(".txt")] + "_OBJVAL.txt"
        f2 = open(filename, "a+")
        f2.write(str(self.info["OBJ_VAL"]) + "\t")
        f2.close()



    def timestamp(self):
        ''' Convert time.time() to proper timestamp in the format <year><month><day><hour><min><sec>'''
        sec = time.time()
        T = time.localtime(sec)
        tstamp = ""
        for i in range(6):
            t = str(T[i])
            if (len(t)==1):
                t = "0" + t
            tstamp +=  t
        return (tstamp)

