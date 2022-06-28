# -*- coding: utf-8 -*-
import sys, inspect
from datetime import datetime
import numpy          as np
import sympy          as syp
import graph_tool.all as gt

########################################################################

class Vertex(object):                                               ##### Vertex object (connectionpoint)
    def __init__(self, id = None, K = None, G = None, C = None, 
                 R = None):
        self.id = id                                                    # Vertex id (starting from 0)
        self.K = K                                                      # Parent object - Knotdiagram
        self.G = G                                                      # Parent object - Graph
        self.C = C                                                      # Parent object - Crossing
        self.R = R                                                      # Parent object - Rope
        self.V_b = None                                                 # Coming from <- Vertex object ('b'efore)
        self.E_b = None                                                 # Edge   from <- ('b'efore)
        self.E_a = None                                                 # Edge   to   -> ('a'fter)
        self.V_a = None                                                 # Going  to   -> Vertex object ('a'fter)
        self.Fl = None                                                  # Face object on the left side
        self.Fr = None                                                  # Face object on the right side
        self.LX = {}                                                    # Dictionary for Latex properties

########################################################################
        
class Edge(object):                                                 ##### Edge object (connection between vertices)
    def __init__(self, id = None, K = None, G = None, C = None, 
                 R = None, V0 = None, V1 = None):
        self.id = id                                                    # Edge id (starting from 0)
        self.K = K                                                      # Parent object - Knotdiagram
        self.G = G                                                      # Parent object - Graph
        self.C = C                                                      # Parent object - Crossing
        self.R = R                                                      # Parent object - Rope
        self.V0 = V0                                                    # Starting vertex object
        self.V1 = V1                                                    # Ending vertex object
        self.Fl = None                                                  # Face object on the left side
        self.Fr = None                                                  # Face object on the right side
        self.LX = {'curve' : None, 'id' : None}                         # Dictionary for Latex properties
        
########################################################################
        
class Face(object):                                                 ##### Face object (region separated by edges)
    def __init__(self, id = None, K = None, G = None, C = None, 
                 R = None):
        self.id = id                                                    # Face id (starting from 0)
        self.K = K                                                      # Parent object - Knotdiagram
        self.G = G                                                      # Parent object - Graph
        self.V = []                                                     # List of vertices     (in order)
        self.E = []                                                     # List of edges        (in order)
        self.Ed = []                                                    # List of edge directions (+1 or -1)
        self.C = []                                                     # List of crossings    (in order)
        self.c = None                                                   # Color of the face - Used for Seifert
        self.cy = False                                                 # True if face is a directed cycle
        self.e = None                                                   # Number of edges
        self.df0 = None                                                 # Distance from f0 - outside face
        self.LX = {}                                                    # Dictionary for Latex properties
        
########################################################################
        
class Crossing(object):                                             ##### Crossing object
    def __init__(self, id = None, K = None, G = None):
        self.id = id                                                    # Crossing id (starting from 0)
        self.K = K                                                      # Parent object - Knotdiagram
        self.G = G                                                      # Parent object - Graph
        self.ch = None                                                  # Chirality of the crossing
        self.a = Vertex(id = 0, K = self.K, G = self.G, C = self)       # Coming from vertex object on 'a'bove edge
        self.A = Vertex(id = 2, K = self.K, G = self.G, C = self)       # Going to vertex object on 'A'bove edge
        self.E_a = None                                                 # Edge of 'a'bove part
        self.R_a = None                                                 # Brother object - Rope of 'a'bove part
        self.b = Vertex(id = 1, K = self.K, G = self.G, C = self)       # Coming from vertex object on 'b'elow edge
        self.B = Vertex(id = 3, K = self.K, G = self.G, C = self)       # Going to vertex object on 'B'elow edge
        self.E_b = None                                                 # Edge of 'b'elow part
        self.R_b = None                                                 # Brother object - Rope of 'b'elow part
        self.F1 = Face(K = self.K, G = self.G, C = self)                # Face object of region F1
        self.F2 = Face(K = self.K, G = self.G, C = self)                # Face object of region F2
        self.F3 = Face(K = self.K, G = self.G, C = self)                # Face object of region F3
        self.F4 = Face(K = self.K, G = self.G, C = self)                # Face object of region F4
        
########################################################################

class Rope(object):                                                 ##### Rope object
    def __init__(self, id = None, K = None, G = None):
        self.id = id                                                    # Rope id (starting from 0)
        self.K = K                                                      # Parent object - Knotdiagram
        self.G = G                                                      # Parent object - Graph
        self.C = []                                                     # List of crossings objects
        self.V = []                                                     # List of vertex objects
        self.E = []                                                     # List of edge objects

########################################################################
        
class Graph(object):                                                ##### Graph object
    def __init__(self, K = None, c = None, r = None):
        self.K = K                                                      # Parent object - Knotdiagram
        self.c = c                                                      # Number of crossings
        self.r = r                                                      # Number of ropes
        self.C = []                                                     # List of crossing objects (c)
        self.R = []                                                     # List of rope objects (r)
        self.V = []                                                     # List of vertex objects
        self.E = []                                                     # List of edge objects
        self.F = []                                                     # List of face objects
        self.f = None                                                   # Number of face objects
        self.f0= None                                                   # The outside face
        self.n = 0                                                      # Largest of number of face edges
        self.A = None                                                   # Adjacency matrix
        self.sA = None                                                  # Adjacency matrix with (sympy object) A B
        self.dA = None                                                  # Adjacency matrix (directed only)
        self.fA = None                                                  # Adjacency matrix of dual (face <-> vertex)
        #self.I = None                                                   # Incidence matrix
        #self.L = None                                                   # List of layers
        #self.D = None                                                   # Dual graph object
        #self.S = None                                                   # Seifert surface object
        self.LX = {}                                                    # Dictionary for Latex properties
        
########################################################################
        
class KnotDiagram(object):                                          ##### Main object
    def __init__(self, name = None, c = [0], r = [0]):
        self.name = name                                                # The name of the object
        self.c = c                                                      # Number of crossings
        self.r = r                                                      # Number of ropes
        self.G = Graph(K = self, c = self.c, r = self.r)                # Graph object
        self.Gr = None                                                  # Graph-tool object
        self.G_v_e = None                                               # Graph-tool object where V = V + E
        self.g = Graph(K = self, c = self.c, r = self.r)                # Crossings graph object
        self.gr = None                                                  # Crossings graph-tool object
        self.P = None                                                   # Polynomial (sympy object)
        self.W = ''                                                     # List of warnings
        self.w = False                                                  # True if warnings has to be printed
        self.verbose = True                                             # Boolean for printing settings
        self.input = None                                               # Will be set to the filename if given
        self.L = None                                                   # Latex code string for drawing
        self.Lpos = None                                                # Positions of vertices (graph_tool)
    
    ####################################################################
    
    def PrintWarnings(self):
        if self.verbose:                                                # Is True by default
            if self.W != '':                                            # Print only if conflicts were found
                self.W += ('Conflicts in ' +
                           '"{}"'.format(self.name) + 
                           '\n')
                print(self.W)
    
    ####################################################################
    
    def GenerateGraphs(self, InPut = None, SP = None, CM = None, 
                       to_do = None):
        if InPut:                                                       # If InPut is given, parse it
            self.input = InPut
            if self.name == None:
                self.name = InPut                                       # If name was not given, name by default
            SP, CM = self.ParseInPut(InPut)
        else:
            self.W += ('Warning!!! - InPut was not given to ' + 
                       '"{}()" '.format(inspect.stack()[0][3]) + 
                       '- {}.\n'.format(datetime.now()))
            self.PrintWarnings()
        
        if not all([type(SP) == type(None), type(CM) == type(None)]):
            self.GenerateG(SP, CM)                                      # If SP and CM is not None, calculate Graph
        else:
            self.W += ('Warning!!! - ' + 
                       'SP (StartingPoints) and CM (ConnectionsMatrix)'+ 
                       ' could not be found. Were not given and could' + 
                       ' not be parsed from InPut.\n')
            self.PrintWarnings()
        
        if to_do != None:
            if to_do == 'all':
                self.getAM()
                self.getsAM()
                self.getPolynomial()
                self.findFaces()
                self.getGVE()
                self.genLatex()
            elif to_do == 'graph':
                self.getAM(symp = False)
                self.getsAM(symp = False)
                self.findFaces()
                self.getGVE()
                self.genLatex()
    
    ####################################################################
            
    def ParseInPut(self, InPut):
        try:                                                            # Try to read from the given InPut
            with open(InPut, 'r') as f:
                Text = f.readlines()
        except Exception as e:
            self.W += ('{}\n'.format(e) + 
                       'Could not open InPut = ' +
                       '"{}"'.format(InPut) + 
                       ' in ' + 
                       '"{}()" '.format(inspect.stack()[0][3]) + 
                       '- {}.\n'.format(datetime.now()))
            self.PrintWarnings()
            return
        
        SText = []                                                      # Strip lines and check if empty
        for line in Text:                                               # because crossing number will depend
            if line.strip() != '':                                      # on the number of lines
                SText.append(line.strip())
                
        if len(SText) < 2:                                              # There are no crossings in this case
            self.W += ('Invalid file format. - ' +
                       'Number of non empty lines is less than 2.\n' + 
                       'Please read the documentation.\n')
            self.PrintWarnings()
        
        Header = SText[0]                                               # Line  that contains startingpoints
        Body = SText[1:]                                                # Lines that contain the crossing information
        
        S_P = []                                                        # List of starting points <- Header
        C_M = []                                                        # List of connections     <- Body
        
        #Testing for integer values in the first line of inputfile      ############################################
        for val in Header.split():
            try:
                n = int(val.strip())
                S_P.append(n)
            except Exception as e:
                self.W += ('{}\n'.format(e) + 
                           'Could not turn ' + 
                           '"{}"'.format(val) + 
                           ' value from first line of input ' + 
                           'into an integer in ' + 
                           '"{}()" '.format(inspect.stack()[0][3]) + 
                           '- {}.\n'.format(datetime.now()))
                self.PrintWarnings()
                return
        
        S_P = np.asarray(S_P)
        
        #Testing for even number of values in Header                    # Header is the first line of InPut ########
        if len(Header.split()) % 2 != 0:
            self.W += ('Odd number of values in first line of input ' + 
                       '"{}()" '.format(inspect.stack()[0][3]) + 
                       '- {}.\n'.format(datetime.now()))
            self.PrintWarnings()
            sys.exit(0)
        
        S_P = S_P.reshape((int(len(S_P)/2), 2))
        
        #Testing for integer values and length of 8 in Body             # Body is the connections part of InPut ####
        for l, L in enumerate(Body):
            if len(L.split()) != 8:                                     # Length check
                self.W += ('"{}"'.format(len(L.split())) + 
                           ' is the number of values in line ' + 
                           '"{}"'.format(l + 2) + 
                           ' of input instead of 8, ' + 
                           '"{}()" '.format(inspect.stack()[0][3]) + 
                           '- {}.\n'.format(datetime.now()))
                self.PrintWarnings()
                return
            
            ns = []
            for val in L.split():                                       # Integer check
                try:
                    n = int(val.strip())
                    ns.append(n)
                except Exception as e:
                    self.W += ('{}\n'.format(e) + 
                               'Could not turn ' + 
                               '"{}"'.format(val) + 
                               ' value from line ' + 
                               '"{}"'.format(l + 2) + 
                               ' of input into an integer in ' + 
                               '"{}()" '.format(inspect.stack()[0][3]) + 
                               '- {}.\n'.format(datetime.now()))
                    self.PrintWarnings()
                    return
            
            C_M.append(ns)
        
        C_M = np.asarray(C_M)
        C_M = C_M.reshape((int(len(C_M)), 4, 2))
        
        return(S_P, C_M)
    
    ####################################################################
    
    def CheckSP(self, SP):
        for p in SP:
            if len(p) != 2:                                             # Length check
                self.W += ('Length of point ' + 
                           '"{}"'.format(p) + 
                           ' in SP is not 2 ' + 
                           '"{}()" '.format(inspect.stack()[0][3]) + 
                           '- {}.\n'.format(datetime.now()))
                self.PrintWarnings()
                return
            if not all(np.equal(np.mod(p, 1), 0)):                      # Integer check
                self.W += ('Some values in "SP" are not integers ' + 
                           '"{}()" '.format(inspect.stack()[0][3]) + 
                           '- {}.\n'.format(datetime.now()))
                self.PrintWarnings()
                return
    
    ####################################################################
    
    def CheckCM(self, CM):
        c_nr = len(CM)
        Check = np.zeros(shape = (c_nr, 4, 2))
        Duplicates = []
        
        for cr in CM:                                                   # cr = cossing, CM = connectivity matrix
            if len(cr) != 4:                                            # Length check
                self.W += ('Length of crossing array' +  
                           '"{}"'.format(cr) + 
                           ' in CM is not 4 ' + 
                           '"{}()" '.format(inspect.stack()[0][3]) + 
                           '- {}.\n'.format(datetime.now()))
                self.PrintWarnings()
                return
            for ve in cr:                                               # ve = vertex, cr = crossing
                if len(ve) != 2:                                        # Length check
                    self.W += ('Length of vertex array ' + 
                               '"{}"'.format(ve) + 
                               ' in CM is not 2 ' + 
                               '"{}()" '.format(inspect.stack()[0][3]) + 
                               '- {}.\n'.format(datetime.now()))
                    self.PrintWarnings()
                    return
                if not all(np.equal(np.mod(ve, 1), 0)):                 # Integer check
                    self.W += ('Some values in "CM" are not integers ' + 
                               '"{}()" '.format(inspect.stack()[0][3]) + 
                               '- {}.\n'.format(datetime.now()))
                    self.PrintWarnings()
                    return
                
                if Check[ve[0] - 1][ve[1] - 1][0] != 0:                 # Duplicate check #1
                    Duplicates.append(ve.tolist())
                Check[ve[0] - 1][ve[1] - 1][0] = ve[0]
                Check[ve[0] - 1][ve[1] - 1][1] = ve[1]
                
        if np.count_nonzero(Check) != 8 * c_nr:                         # Duplicate check #2
            self.W += ('Duplicate vertex (vertices) ' + 
                       '"{}"'.format(Duplicates) + 
                       ' found in "CM" ' + 
                       '"{}()" '.format(inspect.stack()[0][3]) + 
                       '- {}.\n'.format(datetime.now()))
            self.w = True
    
    ####################################################################
                    
    def GenerateG(self, SP, CM):
                                                                        # Checking values
        self.CheckSP(SP)
        self.CheckCM(CM)
                                                                        # Setting number of ropes and of crossings
        self.r[0] = len(SP)
        self.c[0] = len(CM)
                                                                        # Creating rope objects
        for ro, Ro in enumerate(SP):
            RO = Rope(id = ro, K = self, G = self.G)
            self.G.R.append(RO)
                                                                        # Creating crossing objects
        for cr, Cr in enumerate(CM):
            CR = Crossing(id = cr, K = self, G = self.G)
            self.G.C.append(CR)
            self.G.V.append(self.G.C[-1].a)
            self.G.V.append(self.G.C[-1].b)
            self.G.V.append(self.G.C[-1].A)
            self.G.V.append(self.G.C[-1].B)
        
        C_d = {}                                                        # Connections dictionary
                                                                        # Loop over crossings
        for cr, Cr in enumerate(CM):                                    # to add connections to a dictionary
            C_d[(cr + 1, 1)] = tuple(Cr[0].tolist())
            C_d[(cr + 1, 2)] = tuple(Cr[1].tolist())
            C_d[(cr + 1, 3)] = tuple(Cr[2].tolist())
            C_d[(cr + 1, 4)] = tuple(Cr[3].tolist())
        
        Vc_d = {1:3, 2:4, 3:1, 4:2}                                     # Vertex connections dictionary
        Rope_paths = []
                                                                        # Loop over the crossings dictionary
        for Sp in SP:                                                   # to create rope paths
            This_path = [Sp.tolist(), [Sp[0], Vc_d[Sp[1]]]]
            self.setChirality(list(Sp))
            
            while This_path[0] != list(C_d[tuple(This_path[-1])]):
                This_path.append(list(C_d[tuple(This_path[-1])]))
                This_path.append([list(C_d[tuple(This_path[-2])])[0],
                                  Vc_d[list(C_d[tuple(This_path[-2])])[1]]])
                                                                        # Setting chirality of the crossing
                self.setChirality(This_path[-2])
                
                if len(This_path) > self.c[0] * 4:
                    self.W += ('Infinite loop found ' + 
                               'when following rope path ' +
                               '"{}" '.format(This_path) + 
                               '"{}()" '.format(inspect.stack()[0][3]) + 
                               '- {}.\n'.format(datetime.now()))
                    self.w = True
                    break                                               # Break if stuck in an infinite loop
            
            Rope_paths.append(This_path)
        
        for r, path in enumerate(Rope_paths):
            for i in range(len(path) - 1):
                self.setVE(path[i], path[i + 1], r)
            
            self.setVE(path[-1], path[0], r)
                                                                        # Printing errors, warnings if any
        
        if self.w:
            self.PrintWarnings()
        
    ####################################################################
    
    def setVE(self, v0, v1, r):
        w = False
        
        if self.G.E == []:                                              # Setting edge id
            Id = 0
        else:
            Id = self.G.E[-1].id + 1
            
        if v0[1] == 1:
            if v0[0] == v1[0] and v1[1] == 3:
                self.G.C[v0[0] - 1].a.V_a = self.G.C[v1[0] - 1].A       # a.V_a --> A
                self.G.C[v1[0] - 1].A.V_b = self.G.C[v0[0] - 1].a       # a --> A.V_b
                ed = Edge(id = Id, K = self, G = self.G, 
                          C  = self.G.C[v0[0] - 1],
                          R  = self.G.R[r],
                          V0 = self.G.C[v0[0] - 1].a,
                          V1 = self.G.C[v1[0] - 1].A)
                self.G.C[v0[0] - 1].a.R = self.G.R[r]                   # C[v0].a.R = R
                self.G.C[v1[0] - 1].A.R = self.G.R[r]                   # C[v1].A.R = R
                self.G.E.append(ed)                                     # G.E[add E]
                self.G.R[r].E.append(self.G.E[-1])                      # R.E[add E]
                self.G.C[v0[0] - 1].a.E_a = self.G.E[-1]                # a.E_a
                self.G.C[v1[0] - 1].A.E_b = self.G.E[-1]                # A.E_b
                self.G.R[r].C.append(self.G.C[v0[0] - 1].a.C)           # R.C[add C]
                self.G.R[r].V.append(self.G.C[v0[0] - 1].a)             # R.V[add V0]
                self.G.R[r].V.append(self.G.C[v1[0] - 1].A)             # R.V[add V1]
                self.G.C[v0[0] - 1].E_a = self.G.E[-1]                  # C.E_a
                self.G.C[v0[0] - 1].R_a = self.G.R[r]                   # C.R_a
            else:
                w = True
        elif v0[1] == 2:
            if v1[1] == 1:
                self.G.C[v0[0] - 1].b.V_a = self.G.C[v1[0] - 1].a       # b.V_a --> a
                self.G.C[v1[0] - 1].a.V_b = self.G.C[v0[0] - 1].b       # b --> a.V_b
                ed = Edge(id = Id, K = self, G = self.G, 
                          C  = self.G.C[v0[0] - 1],
                          R  = self.G.R[r],
                          V0 = self.G.C[v0[0] - 1].b,
                          V1 = self.G.C[v1[0] - 1].a)
                self.G.C[v0[0] - 1].b.R = self.G.R[r]                   # C[v0].b.R = R
                self.G.C[v1[0] - 1].a.R = self.G.R[r]                   # C[v1].a.R = R
                self.G.E.append(ed)                                     # G.E[add E]
                self.G.R[r].E.append(self.G.E[-1])                      # R.E[add E]
                self.G.C[v0[0] - 1].b.E_a = self.G.E[-1]                # b.E_a
                self.G.C[v1[0] - 1].a.E_b = self.G.E[-1]                # a.E_b
            elif v1[1] == 2:
                self.G.C[v0[0] - 1].b.V_a = self.G.C[v1[0] - 1].b       # b.V_a --> b
                self.G.C[v1[0] - 1].b.V_b = self.G.C[v0[0] - 1].b       # b --> b.V_b
                ed = Edge(id = Id, K = self, G = self.G,  
                          C  = self.G.C[v0[0] - 1],
                          R  = self.G.R[r],
                          V0 = self.G.C[v0[0] - 1].b,
                          V1 = self.G.C[v1[0] - 1].b)
                self.G.C[v0[0] - 1].b.R = self.G.R[r]                   # C[v0].b.R = R
                self.G.C[v1[0] - 1].b.R = self.G.R[r]                   # C[v1].b.R = R
                self.G.E.append(ed)                                     # G.E[add E]
                self.G.R[r].E.append(self.G.E[-1])                      # R.E[add E]
                self.G.C[v0[0] - 1].b.E_a = self.G.E[-1]                # b.E_a
                self.G.C[v1[0] - 1].b.E_b = self.G.E[-1]                # b.E_b
            elif v1[1] == 3:
                w = True
            elif v1[1] == 4:
                self.G.C[v0[0] - 1].b.V_a = self.G.C[v1[0] - 1].B       # b.V_a --> B
                self.G.C[v1[0] - 1].B.V_b = self.G.C[v0[0] - 1].b       # b --> B.V_b
                ed = Edge(id = Id, K = self, G = self.G, 
                          C  = self.G.C[v0[0] - 1],
                          R  = self.G.R[r],
                          V0 = self.G.C[v0[0] - 1].b,
                          V1 = self.G.C[v1[0] - 1].B)
                self.G.C[v0[0] - 1].b.R = self.G.R[r]                   # C[v0].b.R = R
                self.G.C[v1[0] - 1].B.R = self.G.R[r]                   # C[v1].B.R = R
                self.G.E.append(ed)                                     # G.E[add E]
                self.G.R[r].E.append(self.G.E[-1])                      # R.E[add E]
                self.G.C[v0[0] - 1].b.E_a = self.G.E[-1]                # b.E_a
                self.G.C[v1[0] - 1].B.E_b = self.G.E[-1]                # B.E_b
                if v0[0] == v1[0]:
                    self.G.R[r].C.append(self.G.C[v0[0] - 1].b.C)       # R.C[add C]
                    self.G.R[r].V.append(self.G.C[v0[0] - 1].b)         # R.V[add V0]
                    self.G.R[r].V.append(self.G.C[v1[0] - 1].B)         # R.V[add V1]
                    self.G.C[v0[0] - 1].E_b = self.G.E[-1]              # C.E_b
                    self.G.C[v0[0] - 1].R_b = self.G.R[r]               # C.R_b
        elif v0[1] == 3:
            if v1[1] == 1:
                self.G.C[v0[0] - 1].A.V_a = self.G.C[v1[0] - 1].a       # A.V_a --> a
                self.G.C[v1[0] - 1].a.V_b = self.G.C[v0[0] - 1].A       # A --> a.V_b
                ed = Edge(id = Id, K = self, G = self.G, 
                          C  = self.G.C[v0[0] - 1],
                          R  = self.G.R[r],
                          V0 = self.G.C[v0[0] - 1].A,
                          V1 = self.G.C[v1[0] - 1].a)
                self.G.C[v0[0] - 1].A.R = self.G.R[r]                   # C[v0].A.R = R
                self.G.C[v1[0] - 1].a.R = self.G.R[r]                   # C[v1].a.R = R
                self.G.E.append(ed)                                     # G.E[add E]
                self.G.R[r].E.append(self.G.E[-1])                      # R.E[add E]
                self.G.C[v0[0] - 1].A.E_a = self.G.E[-1]                # A.E_a
                self.G.C[v1[0] - 1].a.E_b = self.G.E[-1]                # a.E_b
            elif v1[1] == 2:
                self.G.C[v0[0] - 1].A.V_a = self.G.C[v1[0] - 1].b       # A.V_a --> b
                self.G.C[v1[0] - 1].b.V_b = self.G.C[v0[0] - 1].A       # A --> b.V_b
                ed = Edge(id = Id, K = self, G = self.G, 
                          C  = self.G.C[v0[0] - 1],
                          R  = self.G.R[r],
                          V0 = self.G.C[v0[0] - 1].A,
                          V1 = self.G.C[v1[0] - 1].b)
                self.G.C[v0[0] - 1].A.R = self.G.R[r]                   # C[v0].A.R = R
                self.G.C[v1[0] - 1].b.R = self.G.R[r]                   # C[v1].b.R = R
                self.G.E.append(ed)                                     # G.E[add E]
                self.G.R[r].E.append(self.G.E[-1])                      # R.E[add E]
                self.G.C[v0[0] - 1].A.E_a = self.G.E[-1]                # A.E_a
                self.G.C[v1[0] - 1].b.E_b = self.G.E[-1]                # b.E_b
            elif v1[1] == 3:
                w = True
            elif v1[1] == 4:
                self.G.C[v0[0] - 1].A.V_a = self.G.C[v1[0] - 1].B       # A.V_a --> B
                self.G.C[v1[0] - 1].B.V_b = self.G.C[v0[0] - 1].A       # A --> B.V_b
                ed = Edge(id = Id, K = self, G = self.G, 
                          C  = self.G.C[v0[0] - 1],
                          R  = self.G.R[r],
                          V0 = self.G.C[v0[0] - 1].A,
                          V1 = self.G.C[v1[0] - 1].B)
                self.G.C[v0[0] - 1].A.R = self.G.R[r]                   # C[v0].A.R = R
                self.G.C[v1[0] - 1].B.R = self.G.R[r]                   # C[v1].B.R = R
                self.G.E.append(ed)                                     # G.E[add E]
                self.G.R[r].E.append(self.G.E[-1])                      # R.E[add E]
                self.G.C[v0[0] - 1].A.E_a = self.G.E[-1]                # A.E_a
                self.G.C[v1[0] - 1].B.E_b = self.G.E[-1]                # B.E_b
        elif v0[1] == 4:
            if v1[1] == 1:
                self.G.C[v0[0] - 1].B.V_a = self.G.C[v1[0] - 1].a       # B.V_a --> a
                self.G.C[v1[0] - 1].a.V_b = self.G.C[v0[0] - 1].B       # B --> a.V_b
                ed = Edge(id = Id, K = self, G = self.G, 
                          C  = self.G.C[v0[0] - 1],
                          R  = self.G.R[r],
                          V0 = self.G.C[v0[0] - 1].B,
                          V1 = self.G.C[v1[0] - 1].a)
                self.G.C[v0[0] - 1].B.R = self.G.R[r]                   # C[v0].B.R = R
                self.G.C[v1[0] - 1].a.R = self.G.R[r]                   # C[v1].a.R = R
                self.G.E.append(ed)                                     # G.E[add E]
                self.G.R[r].E.append(self.G.E[-1])                      # R.E[add E]
                self.G.C[v0[0] - 1].B.E_a = self.G.E[-1]                # B.E_a
                self.G.C[v1[0] - 1].a.E_b = self.G.E[-1]                # a.E_b
            elif v1[1] == 2:
                self.G.C[v0[0] - 1].B.V_a = self.G.C[v1[0] - 1].b       # B.V_a --> b
                self.G.C[v1[0] - 1].b.V_b = self.G.C[v0[0] - 1].B       # B --> b.V_b
                ed = Edge(id = Id, K = self, G = self.G, 
                          C  = self.G.C[v0[0] - 1],
                          R  = self.G.R[r],
                          V0 = self.G.C[v0[0] - 1].B,
                          V1 = self.G.C[v1[0] - 1].b)
                self.G.C[v0[0] - 1].B.R = self.G.R[r]                   # C[v0].B.R = R
                self.G.C[v1[0] - 1].b.R = self.G.R[r]                   # C[v1].b.R = R
                self.G.E.append(ed)                                     # G.E[add E]
                self.G.R[r].E.append(self.G.E[-1])                      # R.E[add E]
                self.G.C[v0[0] - 1].B.E_a = self.G.E[-1]                # B.E_a
                self.G.C[v1[0] - 1].b.E_b = self.G.E[-1]                # b.E_b
                if v0[0] == v1[0]:
                    self.G.R[r].C.append(self.G.C[v0[0] - 1].B.C)       # R.C[add C]
                    self.G.R[r].V.append(self.G.C[v0[0] - 1].B)         # R.V[add V0]
                    self.G.R[r].V.append(self.G.C[v1[0] - 1].b)         # R.V[add V1]
                    self.G.C[v0[0] - 1].E_b = self.G.E[-1]              # C.E_b
                    self.G.C[v0[0] - 1].R_b = self.G.R[r]               # C.R_b
            elif v1[1] == 3:
                w = True
            elif v1[1] == 4:
                self.G.C[v0[0] - 1].B.V_a = self.G.C[v1[0] - 1].B       # B.V_a --> B
                self.G.C[v1[0] - 1].B.V_b = self.G.C[v0[0] - 1].B       # B --> B.V_b
                ed = Edge(id = Id, K = self, G = self.G, 
                          C  = self.G.C[v0[0] - 1],
                          R  = self.G.R[r],
                          V0 = self.G.C[v0[0] - 1].B,
                          V1 = self.G.C[v1[0] - 1].B)
                self.G.C[v0[0] - 1].B.R = self.G.R[r]                   # C[v0].B.R = R
                self.G.C[v1[0] - 1].B.R = self.G.R[r]                   # C[v1].B.R = R
                self.G.E.append(ed)                                     # G.E[add E]
                self.G.R[r].E.append(self.G.E[-1])                      # R.E[add E]
                self.G.C[v0[0] - 1].B.E_a = self.G.E[-1]                # B.E_a
                self.G.C[v1[0] - 1].B.E_b = self.G.E[-1]                # B.E_b
        else:
            w = True
        
        if w:
            self.W += ('Conflict ' + 
                       '"{} {}"'.format(v0[0], v0[1]) + 
                       ' --> ' + 
                       '"{} {}"'.format(v1[0], v1[1]) + 
                       ' when following rope path ' +
                       '"{}()" '.format(inspect.stack()[0][3]) + 
                       '- {}.\n'.format(datetime.now()))
            self.w = True
        
    ####################################################################
    
    def setChirality(self, ve):
        #Every crossing will divide the plane into 4 regions,
        #they will be numbered as 0, 1, 2, 3 (a, b, A, B)
        #Until the numbers will be determined, every crossing will have 
        #4 regions, called F1, F2, F3, F4 and will always be positioned 
        #as the figure shows then will be assigned to a Face object
        #
        #              b                                     a
        #              1                                     0
        #          F1  ||  F2                           F4   ||   F1
        #              ||                                    ||
        #      a 0 ========> 2 A                     B 3 ====||===> 1 b
        #              ||                                    ||
        #          F4  ||  F3                                ||
        #              V                                F3   V    F2
        # Ch -1        3                       Ch +1         2
        # Left handed  B                       Right handed  A
        
        if ve[1] == 1:
            pass
        elif ve[1] == 2:
            self.G.C[ve[0] - 1].ch = -1
            self.G.C[ve[0] - 1].B.V_b = self.G.C[ve[0] - 1].b           # Vertex before of vertex 3 is 1
            self.G.C[ve[0] - 1].b.V_a = self.G.C[ve[0] - 1].B           # Vertex after  of vertex 1 is 3
        elif ve[1] == 4:
            self.G.C[ve[0] - 1].ch = 1
            self.G.C[ve[0] - 1].B.V_a = self.G.C[ve[0] - 1].b           # Vertex before of vertex 3 is 1
            self.G.C[ve[0] - 1].b.V_b = self.G.C[ve[0] - 1].B           # Vertex after  of vertex 1 is 3
        else:
            self.W += ('Value error ' + 
                       '"{}"'.format(ve[1]) + 
                       ', it should be 1, 2 or 4, ' + 
                       '"{}()" '.format(inspect.stack()[0][3]) + 
                       '- {}.\n'.format(datetime.now()))
            self.w = True
    
    ####################################################################
    
    def getAM(self, symp = True):                                       # Generate the adjacency matrices
        AM = np.zeros((self.c[0]*4, self.c[0]*4), dtype=int)            # without 
        dAM = np.zeros((self.c[0]*4, self.c[0]*4), dtype=int)           # and with symbols
        
        if symp:
            sAM = syp.Matrix(AM)
            A, B = syp.symbols('A B')
        
        for cr in self.G.C:
            AM[cr.id * 4 + cr.a.id, 
               cr.a.V_b.C.id * 4 + cr.a.V_b.id] =-1
            AM[cr.id * 4 + cr.A.id, 
               cr.A.V_a.C.id * 4 + cr.A.V_a.id] = 1
            dAM[cr.id * 4 + cr.A.id, 
                cr.A.V_a.C.id * 4 + cr.A.V_a.id] = 1
            
            if symp:
                sAM[cr.id * 4 + cr.a.id, 
                    cr.a.V_b.C.id * 4 + cr.a.V_b.id] =-A
                sAM[cr.id * 4 + cr.A.id, 
                    cr.A.V_a.C.id * 4 + cr.A.V_a.id] = A
            
            if cr.ch == 1:
                AM[cr.id * 4 + cr.b.id, 
                   cr.b.V_a.C.id * 4 + cr.b.V_a.id] = 1
                AM[cr.id * 4 + cr.B.id, 
                   cr.B.V_b.C.id * 4 + cr.B.V_b.id] =-1
                dAM[cr.id * 4 + cr.b.id, 
                    cr.b.V_a.C.id * 4 + cr.b.V_a.id] = 1
                
                if symp:
                    sAM[cr.id * 4 + cr.b.id, 
                        cr.b.V_a.C.id * 4 + cr.b.V_a.id]= B
                    sAM[cr.id * 4 + cr.B.id, 
                        cr.B.V_b.C.id * 4 + cr.B.V_b.id]=-B
                
            elif cr.ch == -1:
                AM[cr.id * 4 + cr.b.id, 
                   cr.b.V_b.C.id * 4 + cr.b.V_b.id] =-1
                AM[cr.id * 4 + cr.B.id, 
                   cr.B.V_a.C.id * 4 + cr.B.V_a.id] = 1
                dAM[cr.id * 4 + cr.B.id, 
                    cr.B.V_a.C.id * 4 + cr.B.V_a.id] = 1
                
                if symp:
                    sAM[cr.id * 4 + cr.b.id, 
                        cr.b.V_b.C.id * 4 + cr.b.V_b.id]=-B
                    sAM[cr.id * 4 + cr.B.id, 
                        cr.B.V_a.C.id * 4 + cr.B.V_a.id]= B
                
        self.G.A = AM
        self.G.dA = dAM
        
        if symp:
            self.G.sA = sAM

    ####################################################################
    
    def getsAM(self, symp = True):                                      # Generate the small adjacency matrices
        am = np.zeros((self.c[0], self.c[0]), dtype=int)                # without
        dam = np.zeros((self.c[0], self.c[0]), dtype=int)               # and with symbols
        
        if symp:
            sam = syp.Matrix(am)
        
        for i in range(self.G.c[0]):
            for j in range(self.G.c[0]):
                val = (abs(self.G.A[i * 4    ][j * 4    ]) + 
                       abs(self.G.A[i * 4 + 1][j * 4    ]) +
                       abs(self.G.A[i * 4 + 2][j * 4    ]) + 
                       abs(self.G.A[i * 4 + 3][j * 4    ]) +
                       abs(self.G.A[i * 4    ][j * 4 + 1]) + 
                       abs(self.G.A[i * 4 + 1][j * 4 + 1]) +
                       abs(self.G.A[i * 4 + 2][j * 4 + 1]) + 
                       abs(self.G.A[i * 4 + 3][j * 4 + 1]) +
                       abs(self.G.A[i * 4    ][j * 4 + 2]) + 
                       abs(self.G.A[i * 4 + 1][j * 4 + 2]) +
                       abs(self.G.A[i * 4 + 2][j * 4 + 2]) + 
                       abs(self.G.A[i * 4 + 3][j * 4 + 2]) +
                       abs(self.G.A[i * 4    ][j * 4 + 3]) + 
                       abs(self.G.A[i * 4 + 1][j * 4 + 3]) +
                       abs(self.G.A[i * 4 + 2][j * 4 + 3]) + 
                       abs(self.G.A[i * 4 + 3][j * 4 + 3]))
                
                dval = (self.G.dA[i * 4    ][j * 4    ] + 
                        self.G.dA[i * 4 + 1][j * 4    ] +
                        self.G.dA[i * 4 + 2][j * 4    ] + 
                        self.G.dA[i * 4 + 3][j * 4    ] +
                        self.G.dA[i * 4    ][j * 4 + 1] + 
                        self.G.dA[i * 4 + 1][j * 4 + 1] +
                        self.G.dA[i * 4 + 2][j * 4 + 1] + 
                        self.G.dA[i * 4 + 3][j * 4 + 1] +
                        self.G.dA[i * 4    ][j * 4 + 2] + 
                        self.G.dA[i * 4 + 1][j * 4 + 2] +
                        self.G.dA[i * 4 + 2][j * 4 + 2] + 
                        self.G.dA[i * 4 + 3][j * 4 + 2] +
                        self.G.dA[i * 4    ][j * 4 + 3] + 
                        self.G.dA[i * 4 + 1][j * 4 + 3] +
                        self.G.dA[i * 4 + 2][j * 4 + 3] + 
                        self.G.dA[i * 4 + 3][j * 4 + 3])
             
                if symp:
                    sval = (self.G.sA.row(i * 4    )[j * 4    ] + 
                            self.G.sA.row(i * 4 + 1)[j * 4    ] +
                            self.G.sA.row(i * 4 + 2)[j * 4    ] + 
                            self.G.sA.row(i * 4 + 3)[j * 4    ] +
                            self.G.sA.row(i * 4    )[j * 4 + 1] + 
                            self.G.sA.row(i * 4 + 1)[j * 4 + 1] +
                            self.G.sA.row(i * 4 + 2)[j * 4 + 1] + 
                            self.G.sA.row(i * 4 + 3)[j * 4 + 1] +
                            self.G.sA.row(i * 4    )[j * 4 + 2] + 
                            self.G.sA.row(i * 4 + 1)[j * 4 + 2] +
                            self.G.sA.row(i * 4 + 2)[j * 4 + 2] + 
                            self.G.sA.row(i * 4 + 3)[j * 4 + 2] +
                            self.G.sA.row(i * 4    )[j * 4 + 3] + 
                            self.G.sA.row(i * 4 + 1)[j * 4 + 3] +
                            self.G.sA.row(i * 4 + 2)[j * 4 + 3] + 
                            self.G.sA.row(i * 4 + 3)[j * 4 + 3])
                 
                am[i,j] = val
                dam[i,j] = dval
                
                if symp:
                    sam[i,j] = sval
        
        self.g.A = am
        self.g.dA = dam
        
        if symp:
            self.g.sA = sam
    
    ####################################################################
    
    def getPolynomial(self):                                            # Generate the polynomial
        mat = self.g.sA
        
        for cr in self.G.C:
            mat[cr.id, cr.id] = cr.ch
        
        self.P = mat.det().simplify()
    
    ####################################################################
    
    def findFaces(self):                                                # Function to find the faces of the graph
        f_count = 0
        for rope in self.G.R:
            for edge in rope.E:
                if edge.V0.C.id != edge.V1.C.id:
                    if edge.Fl == None:
                        self.findFaceL(edge, edge.V0, edge.V1, f_count)
                        f_count += 1
                    if edge.Fr == None:
                        self.findFaceR(edge, edge.V0, edge.V1, f_count)
                        f_count += 1
                else:
                    if not ((edge.V0.id == 0 and edge.V1.id == 2) or
                            (edge.V0.id == 1 and edge.V1.id == 3) or
                            (edge.V0.id == 2 and edge.V1.id == 0) or
                            (edge.V0.id == 3 and edge.V1.id == 1)):     # Including R1 (Reidemeister 1) case
                        
                        if edge.Fl == None:
                            self.findFaceL(edge, edge.V0, edge.V1, f_count)
                            f_count += 1
                        if edge.Fr == None:
                            self.findFaceR(edge, edge.V0, edge.V1, f_count)
                            f_count += 1
        
        self.getFaceAdjacencyM()
    
    ####################################################################
    
    def findFaceL(self, edge, V0, V1, f_count):                         # Function to find a face going left
        first_vertex = V0
        notsame = True
        
        while notsame:                                                  # While loop
            if V0.id == 0: ############################################## V0.id == 0
                if V1.id == 0: ########################################## V1.id == 0 - Not possible
                    self.W += ('0 <--> 0 conflict at ' + 
                               '"{} {}" - '.format(V0.C.id, V0.id) + 
                               '"{} {}" in '.format(V1.C.id, V1.id) + 
                               '"{}()" '.format(inspect.stack()[0][3]) + 
                               '- {}.\n'.format(datetime.now()))
                    self.w = True
                    notsame = False
                
                elif V1.id == 1: ######################################## V1.id == 1
                    V1.C.F2 = V0.C.F4                                   # F2 => F4
                    V0.C.F4.V.append(V0)
                    V0.C.F4.V.append(V1)
                    V0.C.F4.E.append(edge)
                    V0.C.F4.C.append(V0.C)
                                                                        # Opposite direction only
                    V0.C.F4.Ed.append(-1)
                    edge.Fr = V0.C.F4                                   # Set right faces
                    V0.Fr   = V0.C.F4
                    V1.Fr   = V0.C.F4
                    
                    V0   = V1.C.A                                       # Set new V0
                    edge = V1.C.A.E_a                                   # Set new Edge
                    V1   = edge.V1                                      # Set new V1
                
                elif V1.id == 2: ######################################## V1.id == 2
                    V1.C.F3 = V0.C.F4                                   # F3 => F4
                    V0.C.F4.V.append(V0)
                    V0.C.F4.V.append(V1)
                    V0.C.F4.E.append(edge)
                    V0.C.F4.C.append(V0.C)
                                                                        # Opposite direction only
                    V0.C.F4.Ed.append(-1)
                    edge.Fr = V0.C.F4                                   # Set right faces
                    V0.Fr   = V0.C.F4
                    V1.Fr   = V0.C.F4
                    
                    V0   = V1.C.B                                       # Set new V0
                    if V1.C.ch == 1:                                    # Set new Edge
                        edge = V1.C.B.E_b
                        V1   = edge.V0
                    else:                                               # Set new V1
                        edge = V1.C.B.E_a
                        V1   = edge.V1
                
                elif V1.id == 3: ######################################## V1.id == 3
                    V1.C.F4 = V0.C.F4                                   # F4 => F4
                    V0.C.F4.V.append(V0)
                    V0.C.F4.V.append(V1)
                    V0.C.F4.E.append(edge)
                    V0.C.F4.C.append(V0.C)
                                                                        # Opposite direction only
                    V0.C.F4.Ed.append(-1)
                    edge.Fr = V0.C.F4                                   # Set right faces
                    V0.Fr   = V0.C.F4
                    V1.Fr   = V0.C.F4
                    
                    V0   = V1.C.a                                       # Set new V0
                    edge = V1.C.a.E_b                                   # Set new Edge
                    V1   = edge.V0                                      # Set new V1   
                
                else:
                    self.W += ('Value error! ' + 
                               'Vertex id has to be 0, 1, 2 or 3. ' + 
                               '"{}"'.format(V1.id) + 
                               'was given instead. ' + 
                               '"{}()" '.format(inspect.stack()[0][3]) + 
                               '- {}.\n'.format(datetime.now()))
                    self.w = True
                    notsame = False
                
            elif V0.id == 1: ############################################ V0.id == 1
                if V1.id == 0: ########################################## V1.id == 0
                    V1.C.F1 = V0.C.F1                                   # F1 => F1
                    V0.C.F1.V.append(V0)
                    V0.C.F1.V.append(V1)
                    V0.C.F1.E.append(edge)
                    V0.C.F1.C.append(V0.C)
                                                                        # Same direction only
                    V0.C.F1.Ed.append(1)
                    edge.Fl = V0.C.F1                                   # Set left faces
                    V0.Fl   = V0.C.F1
                    V1.Fl   = V0.C.F1
                    
                    V0   = V1.C.b                                       # Set new V0
                    
                    if V1.C.ch == 1:                                    # Set new Edge
                        edge = V1.C.b.E_a
                        V1   = edge.V1
                    else:                                               # Set new V1
                        edge = V1.C.b.E_b
                        V1   = edge.V0 
                    
                elif V1.id == 1: ######################################## V1.id == 1
                    V1.C.F2 = V0.C.F1                                   # F2 => F1
                    V0.C.F1.V.append(V0)
                    V0.C.F1.V.append(V1)
                    V0.C.F1.E.append(edge)
                    V0.C.F1.C.append(V0.C)
                
                    if edge.V0.C.id == V0.C.id:                         # Same direction
                        V0.C.F1.Ed.append(1)
                        edge.Fl = V0.C.F1                               # Set left faces
                        V0.Fl   = V0.C.F1
                        V1.Fl   = V0.C.F1
                    
                    else:                                               # Opposite direction
                        V0.C.F1.Ed.append(-1)
                        edge.Fr = V0.C.F1                               # Set right faces
                        V0.Fr   = V0.C.F1
                        V1.Fr   = V0.C.F1
                    
                    V0   = V1.C.A                                       # Set new V0
                    edge = V1.C.A.E_a                                   # Set new Edge
                    V1   = edge.V1                                      # Set new V1
                    
                elif V1.id == 2: ######################################## V1.id == 2
                    V1.C.F3 = V0.C.F1                                   # F3 => F1
                    V0.C.F1.V.append(V0)
                    V0.C.F1.V.append(V1)
                    V0.C.F1.E.append(edge)
                    V0.C.F1.C.append(V0.C)
                                                                        # Opposite direction only
                    V0.C.F1.Ed.append(-1)
                    edge.Fr = V0.C.F1                                   # Set right faces
                    V0.Fr   = V0.C.F1
                    V1.Fr   = V0.C.F1
                    
                    V0   = V1.C.B                                       # Set new V0
                    
                    if V1.C.ch == 1:                                    # Set new Edge
                        edge = V1.C.B.E_b
                        V1   = edge.V0
                    else:                                               # Set new V1
                        edge = V1.C.B.E_a
                        V1   = edge.V1
                         
                elif V1.id == 3: ######################################## V1.id == 3
                    V1.C.F4 = V0.C.F1                                   # F4 => F1
                    V0.C.F1.V.append(V0)
                    V0.C.F1.V.append(V1)
                    V0.C.F1.E.append(edge)
                    V0.C.F1.C.append(V0.C)
                
                    if edge.V0.C.id == V0.C.id:                         # Same direction
                        V0.C.F1.Ed.append(1)
                        edge.Fl = V0.C.F1                               # Set left faces
                        V0.Fl   = V0.C.F1
                        V1.Fl   = V0.C.F1
                    
                    else:                                               # Opposite direction
                        V0.C.F1.Ed.append(-1)
                        edge.Fr = V0.C.F1                               # Set right faces
                        V0.Fr   = V0.C.F1
                        V1.Fr   = V0.C.F1
                    
                    V0   = V1.C.a                                       # Set new V0
                    edge = V1.C.a.E_b                                   # Set new Edge
                    V1   = edge.V0                                      # Set new V1
                    
                else:
                    self.W += ('Value error! ' + 
                               'Vertex id has to be 0, 1, 2 or 3. ' + 
                               '"{}"'.format(V1.id) + 
                               'was given instead. ' + 
                               '"{}()" '.format(inspect.stack()[0][3]) + 
                               '- {}.\n'.format(datetime.now()))
                    self.w = True
                    notsame = False
                
            elif V0.id == 2: ############################################ V0.id == 2
                if V1.id == 0: ########################################## V1.id == 0
                    V1.C.F1 = V0.C.F2                                   # F1 => F2
                    V0.C.F2.V.append(V0)
                    V0.C.F2.V.append(V1)
                    V0.C.F2.E.append(edge)
                    V0.C.F2.C.append(V0.C)
                                                                        # Same direction only
                    V0.C.F2.Ed.append(1)
                    edge.Fl = V0.C.F2                                   # Set left faces
                    V0.Fl   = V0.C.F2
                    V1.Fl   = V0.C.F2

                    V0   = V1.C.b                                       # Set new V0
                    
                    if V1.C.ch == 1:                                    # Set new Edge
                        edge = V1.C.b.E_a
                        V1   = edge.V1
                    else:                                               # Set new V1
                        edge = V1.C.b.E_b
                        V1   = edge.V0
                
                elif V1.id == 1: ######################################## V1.id == 1
                    V1.C.F2 = V0.C.F2                                   # F2 => F2
                    V0.C.F2.V.append(V0)
                    V0.C.F2.V.append(V1)
                    V0.C.F2.E.append(edge)
                    V0.C.F2.C.append(V0.C)
                                                                        # Same direction only
                    V0.C.F2.Ed.append(1)
                    edge.Fl = V0.C.F2                                   # Set left faces
                    V0.Fl   = V0.C.F2
                    V1.Fl   = V0.C.F2

                    V0   = V1.C.A                                       # Set new V0
                    edge = V1.C.A.E_a                                   # Set new Edge
                    V1   = edge.V1                                      # Set new V1
                    
                elif V1.id == 2: ######################################## V1.id == 2 - Not possible
                    self.W += ('2 -><- 2 conflict at ' + 
                               '"{} {}" - '.format(V0.C.id, V0.id) + 
                               '"{} {}" in '.format(V1.C.id, V1.id) + 
                               '"{}()" '.format(inspect.stack()[0][3]) + 
                               '- {}.\n'.format(datetime.now()))
                    self.w = True
                    notsame = False
                
                elif V1.id == 3: ######################################## V1.id == 3
                    V1.C.F4 = V0.C.F2                                   # F4 => F2
                    V0.C.F2.V.append(V0)
                    V0.C.F2.V.append(V1)
                    V0.C.F2.E.append(edge)
                    V0.C.F2.C.append(V0.C)
                                                                        # Same direction only
                    V0.C.F2.Ed.append(1)
                    edge.Fl = V0.C.F2                                   # Set left faces
                    V0.Fl   = V0.C.F2
                    V1.Fl   = V0.C.F2

                    V0   = V1.C.a                                       # Set new V0
                    edge = V1.C.a.E_b                                   # Set new Edge
                    V1   = edge.V0                                      # Set new V1
                    
                else:
                    self.W += ('Value error! ' + 
                               'Vertex id has to be 0, 1, 2 or 3. ' + 
                               '"{}"'.format(V1.id) + 
                               'was given instead. ' + 
                               '"{}()" '.format(inspect.stack()[0][3]) + 
                               '- {}.\n'.format(datetime.now()))
                    self.w = True
                    notsame = False
                
            elif V0.id == 3: ############################################ V0.id == 3
                if V1.id == 0: ########################################## V1.id == 0
                    V1.C.F1 = V0.C.F3                                   # F1 => F3
                    V0.C.F3.V.append(V0)
                    V0.C.F3.V.append(V1)
                    V0.C.F3.E.append(edge)
                    V0.C.F3.C.append(V0.C)
                                                                        # Same direction only
                    V0.C.F3.Ed.append(1)
                    edge.Fl = V0.C.F3                                   # Set left faces
                    V0.Fl   = V0.C.F3
                    V1.Fl   = V0.C.F3

                    V0   = V1.C.b                                       # Set new V0
                    
                    if V1.C.ch == 1:                                    # Set new Edge
                        edge = V1.C.b.E_a
                        V1   = edge.V1
                    else:                                               # Set new V1
                        edge = V1.C.b.E_b
                        V1   = edge.V0
                    
                elif V1.id == 1: ######################################## V1.id == 1
                    V1.C.F2 = V0.C.F3                                   # F2 => F3
                    V0.C.F3.V.append(V0)
                    V0.C.F3.V.append(V1)
                    V0.C.F3.E.append(edge)
                    V0.C.F3.C.append(V0.C)
                
                    if edge.V0.C.id == V0.C.id:                         # Same direction
                        V0.C.F3.Ed.append(1)
                        edge.Fl = V0.C.F3                               # Set left faces
                        V0.Fl   = V0.C.F3
                        V1.Fl   = V0.C.F3
                    
                    else:                                               # Opposite direction
                        V0.C.F3.Ed.append(-1)
                        edge.Fr = V0.C.F3                               # Set right faces
                        V0.Fr   = V0.C.F3
                        V1.Fr   = V0.C.F3
                    
                    V0   = V1.C.A                                       # Set new V0
                    edge = V1.C.A.E_a                                   # Set new Edge
                    V1   = edge.V1                                      # Set new V1
                    
                elif V1.id == 2: ######################################## V1.id == 2
                    V1.C.F3 = V0.C.F3                                   # F3 => F3
                    V0.C.F3.V.append(V0)
                    V0.C.F3.V.append(V1)
                    V0.C.F3.E.append(edge)
                    V0.C.F3.C.append(V0.C)
                                                                        # Opposite direction only
                    V0.C.F3.Ed.append(-1)
                    edge.Fr = V0.C.F3                                   # Set right faces
                    V0.Fr   = V0.C.F3
                    V1.Fr   = V0.C.F3

                    V0   = V1.C.B                                       # Set new V0
                    
                    if V1.C.ch == 1:                                    # Set new Edge
                        edge = V1.C.B.E_b
                        V1   = edge.V0
                    else:                                               # Set new V1
                        edge = V1.C.B.E_a
                        V1   = edge.V1
                    
                elif V1.id == 3: ######################################## V1.id == 3
                    V1.C.F4 = V0.C.F3                                   # F4 => F3
                    V0.C.F3.V.append(V0)
                    V0.C.F3.V.append(V1)
                    V0.C.F3.E.append(edge)
                    V0.C.F3.C.append(V0.C)
                
                    if edge.V0.C.id == V0.C.id:                         # Same direction
                        V0.C.F3.Ed.append(1)
                        edge.Fl = V0.C.F3                               # Set left faces
                        V0.Fl   = V0.C.F3
                        V1.Fl   = V0.C.F3
                    
                    else:                                               # Opposite direction
                        V0.C.F3.Ed.append(-1)
                        edge.Fr = V0.C.F3                               # Set right faces
                        V0.Fr   = V0.C.F3
                        V1.Fr   = V0.C.F3
                    
                    V0   = V1.C.a                                       # Set new V0
                    edge = V1.C.a.E_b                                   # Set new Edge
                    V1   = edge.V0                                      # Set new V1
                    
                else:
                    self.W += ('Value error! ' + 
                               'Vertex id has to be 0, 1, 2 or 3. ' + 
                               '"{}"'.format(V1.id) + 
                               'was given instead. ' + 
                               '"{}()" '.format(inspect.stack()[0][3]) + 
                               '- {}.\n'.format(datetime.now()))
                    self.w = True
                    notsame = False
                
            else:
                self.W += ('Value error! ' + 
                           'Vertex id has to be 0, 1, 2 or 3. ' + 
                           '"{}"'.format(V0.id) + 
                           'was given instead. ' + 
                           '"{}()" '.format(inspect.stack()[0][3]) + 
                           '- {}.\n'.format(datetime.now()))
                self.w = True
                notsame = False
            
            if first_vertex.C.id == V0.C.id:
                first_vertex.Fl.id = f_count
                self.G.F.append(first_vertex.Fl)
                self.G.f = f_count + 1
                
                summ = 0
                for i in first_vertex.Fl.Ed:
                    summ += i
                if summ == len(first_vertex.Fl.Ed) or summ == - len(first_vertex.Fl.Ed):
                    first_vertex.Fl.cy = True
                
                first_vertex.Fl.e = len(first_vertex.Fl.E)
                
                if self.G.n < first_vertex.Fl.e:
                    self.G.n  = first_vertex.Fl.e
                    self.G.f0 = first_vertex.Fl
                elif self.G.n == first_vertex.Fl.e:
                    if (self.G.f0.cy != True) and (first_vertex.Fl.cy == True):
                        self.G.f0 = first_vertex.Fl
                notsame = False
    
    ####################################################################
    
    def findFaceR(self, edge, V0, V1, f_count):                         # Function to find a face going right
        first_vertex = V0
        notsame = True
        
        while notsame:                                                  # While loop
            if V0.id == 0: ############################################## V0.id == 0
                if V1.id == 0: ########################################## V1.id == 0 - Not possible
                    self.W += ('0 <--> 0 conflict at ' + 
                               '"{} {}" - '.format(V0.C.id, V0.id) + 
                               '"{} {}" in '.format(V1.C.id, V1.id) + 
                               '"{}()" '.format(inspect.stack()[0][3]) + 
                               '- {}.\n'.format(datetime.now()))
                    self.w = True
                    notsame = False
                
                elif V1.id == 1: ######################################## V1.id == 1
                    V1.C.F1 = V0.C.F1                                   # F1 => F1
                    V0.C.F1.V.append(V0)
                    V0.C.F1.V.append(V1)
                    V0.C.F1.E.append(edge)
                    V0.C.F1.C.append(V0.C)
                                                                        # Opposite direction only
                    V0.C.F1.Ed.append(-1)
                    edge.Fl = V0.C.F1                                   # Set left faces
                    V0.Fl   = V0.C.F1
                    V1.Fl   = V0.C.F1
                    
                    V0   = V1.C.a                                       # Set new V0
                    edge = V1.C.a.E_b                                   # Set new Edge
                    V1   = edge.V0                                      # Set new V1
                
                elif V1.id == 2: ######################################## V1.id == 2
                    V1.C.F2 = V0.C.F1                                   # F2 => F1
                    V0.C.F1.V.append(V0)
                    V0.C.F1.V.append(V1)
                    V0.C.F1.E.append(edge)
                    V0.C.F1.C.append(V0.C)
                                                                        # Opposite direction only
                    V0.C.F1.Ed.append(-1)
                    edge.Fl = V0.C.F1                                   # Set left faces
                    V0.Fl   = V0.C.F1
                    V1.Fl   = V0.C.F1
                    
                    V0   = V1.C.b                                       # Set new V0
                    if V1.C.ch == 1:                                    # Set new Edge
                        edge = V1.C.b.E_a
                        V1   = edge.V1
                    else:                                               # Set new V1
                        edge = V1.C.b.E_b
                        V1   = edge.V0
                
                elif V1.id == 3: ######################################## V1.id == 3
                    V1.C.F3 = V0.C.F1                                   # F3 => F1
                    V0.C.F1.V.append(V0)
                    V0.C.F1.V.append(V1)
                    V0.C.F1.E.append(edge)
                    V0.C.F1.C.append(V0.C)
                                                                        # Opposite direction only
                    V0.C.F1.Ed.append(-1)
                    edge.Fl = V0.C.F1                                   # Set left faces
                    V0.Fl   = V0.C.F1
                    V1.Fl   = V0.C.F1
                    
                    V0   = V1.C.A                                       # Set new V0
                    edge = V1.C.A.E_a                                   # Set new Edge
                    V1   = edge.V1                                      # Set new V1   
                
                else:
                    self.W += ('Value error! ' + 
                               'Vertex id has to be 0, 1, 2 or 3. ' + 
                               '"{}"'.format(V1.id) + 
                               'was given instead. ' + 
                               '"{}()" '.format(inspect.stack()[0][3]) + 
                               '- {}.\n'.format(datetime.now()))
                    self.w = True
                    notsame = False
                
            elif V0.id == 1: ############################################ V0.id == 1
                if V1.id == 0: ########################################## V1.id == 0
                    V1.C.F4 = V0.C.F2                                   # F4 => F2
                    V0.C.F2.V.append(V0)
                    V0.C.F2.V.append(V1)
                    V0.C.F2.E.append(edge)
                    V0.C.F2.C.append(V0.C)
                                                                        # Same direction only
                    V0.C.F2.Ed.append(1)
                    edge.Fr = V0.C.F2                                   # Set right faces
                    V0.Fr   = V0.C.F2
                    V1.Fr   = V0.C.F2
                    
                    V0   = V1.C.B                                       # Set new V0
                    
                    if V1.C.ch == 1:                                    # Set new Edge
                        edge = V1.C.B.E_b
                        V1   = edge.V0
                    else:                                               # Set new V1
                        edge = V1.C.B.E_a
                        V1   = edge.V1 
                    
                elif V1.id == 1: ######################################## V1.id == 1
                    V1.C.F1 = V0.C.F2                                   # F1 => F2
                    V0.C.F2.V.append(V0)
                    V0.C.F2.V.append(V1)
                    V0.C.F2.E.append(edge)
                    V0.C.F2.C.append(V0.C)
                
                    if edge.V0.C.id == V0.C.id:                         # Same direction
                        V0.C.F2.Ed.append(1)
                        edge.Fr = V0.C.F2                               # Set right faces
                        V0.Fr   = V0.C.F2
                        V1.Fr   = V0.C.F2
                    
                    else:                                               # Opposite direction
                        V0.C.F2.Ed.append(-1)
                        edge.Fl = V0.C.F2                               # Set left faces
                        V0.Fl   = V0.C.F2
                        V1.Fl   = V0.C.F2
                    
                    V0   = V1.C.a                                       # Set new V0
                    edge = V1.C.a.E_b                                   # Set new Edge
                    V1   = edge.V0                                      # Set new V1
                    
                elif V1.id == 2: ######################################## V1.id == 2
                    V1.C.F2 = V0.C.F2                                   # F2 => F2
                    V0.C.F2.V.append(V0)
                    V0.C.F2.V.append(V1)
                    V0.C.F2.E.append(edge)
                    V0.C.F2.C.append(V0.C)
                                                                        # Opposite direction only
                    V0.C.F2.Ed.append(-1)
                    edge.Fl = V0.C.F2                                   # Set left faces
                    V0.Fl   = V0.C.F2
                    V1.Fl   = V0.C.F2
                    
                    V0   = V1.C.b                                       # Set new V0
                    
                    if V1.C.ch == 1:                                    # Set new Edge
                        edge = V1.C.b.E_a
                        V1   = edge.V1
                    else:                                               # Set new V1
                        edge = V1.C.b.E_b
                        V1   = edge.V0
                         
                elif V1.id == 3: ######################################## V1.id == 3
                    V1.C.F3 = V0.C.F2                                   # F3 => F2
                    V0.C.F2.V.append(V0)
                    V0.C.F2.V.append(V1)
                    V0.C.F2.E.append(edge)
                    V0.C.F2.C.append(V0.C)
                
                    if edge.V0.C.id == V0.C.id:                         # Same direction
                        V0.C.F2.Ed.append(1)
                        edge.Fr = V0.C.F2                               # Set right faces
                        V0.Fr   = V0.C.F2
                        V1.Fr   = V0.C.F2
                    
                    else:                                               # Opposite direction
                        V0.C.F2.Ed.append(-1)
                        edge.Fl = V0.C.F2                               # Set left faces
                        V0.Fl   = V0.C.F2
                        V1.Fl   = V0.C.F2
                    
                    V0   = V1.C.A                                       # Set new V0
                    edge = V1.C.A.E_a                                   # Set new Edge
                    V1   = edge.V1                                      # Set new V1
                    
                else:
                    self.W += ('Value error! ' + 
                               'Vertex id has to be 0, 1, 2 or 3. ' + 
                               '"{}"'.format(V1.id) + 
                               'was given instead. ' + 
                               '"{}()" '.format(inspect.stack()[0][3]) + 
                               '- {}.\n'.format(datetime.now()))
                    self.w = True
                    notsame = False
                
            elif V0.id == 2: ############################################ V0.id == 2
                if V1.id == 0: ########################################## V1.id == 0
                    V1.C.F4 = V0.C.F3                                   # F4 => F3
                    V0.C.F3.V.append(V0)
                    V0.C.F3.V.append(V1)
                    V0.C.F3.E.append(edge)
                    V0.C.F3.C.append(V0.C)
                                                                        # Same direction only
                    V0.C.F3.Ed.append(1)
                    edge.Fr = V0.C.F3                                   # Set right faces
                    V0.Fr   = V0.C.F3
                    V1.Fr   = V0.C.F3

                    V0   = V1.C.B                                       # Set new V0
                    
                    if V1.C.ch == 1:                                    # Set new Edge
                        edge = V1.C.B.E_b
                        V1   = edge.V0
                    else:                                               # Set new V1
                        edge = V1.C.B.E_a
                        V1   = edge.V1
                
                elif V1.id == 1: ######################################## V1.id == 1
                    V1.C.F1 = V0.C.F3                                   # F1 => F3
                    V0.C.F3.V.append(V0)
                    V0.C.F3.V.append(V1)
                    V0.C.F3.E.append(edge)
                    V0.C.F3.C.append(V0.C)
                                                                        # Same direction only
                    V0.C.F3.Ed.append(1)
                    edge.Fr = V0.C.F3                                   # Set right faces
                    V0.Fr   = V0.C.F3
                    V1.Fr   = V0.C.F3

                    V0   = V1.C.a                                       # Set new V0
                    edge = V1.C.a.E_b                                   # Set new Edge
                    V1   = edge.V0                                      # Set new V1
                    
                elif V1.id == 2: ######################################## V1.id == 2 - Not possible
                    self.W += ('2 -><- 2 conflict at ' + 
                               '"{} {}" - '.format(V0.C.id, V0.id) + 
                               '"{} {}" in '.format(V1.C.id, V1.id) + 
                               '"{}()" '.format(inspect.stack()[0][3]) + 
                               '- {}.\n'.format(datetime.now()))
                    self.w = True
                    notsame = False
                
                elif V1.id == 3: ######################################## V1.id == 3
                    V1.C.F3 = V0.C.F3                                   # F3 => F3
                    V0.C.F3.V.append(V0)
                    V0.C.F3.V.append(V1)
                    V0.C.F3.E.append(edge)
                    V0.C.F3.C.append(V0.C)
                                                                        # Same direction only
                    V0.C.F3.Ed.append(1)
                    edge.Fr = V0.C.F3                                   # Set right faces
                    V0.Fr   = V0.C.F3
                    V1.Fr   = V0.C.F3

                    V0   = V1.C.A                                       # Set new V0
                    edge = V1.C.A.E_a                                   # Set new Edge
                    V1   = edge.V1                                      # Set new V1
                    
                else:
                    self.W += ('Value error! ' + 
                               'Vertex id has to be 0, 1, 2 or 3. ' + 
                               '"{}"'.format(V1.id) + 
                               'was given instead. ' + 
                               '"{}()" '.format(inspect.stack()[0][3]) + 
                               '- {}.\n'.format(datetime.now()))
                    self.w = True
                    notsame = False
                
            elif V0.id == 3: ############################################ V0.id == 3
                if V1.id == 0: ########################################## V1.id == 0
                    V1.C.F4 = V0.C.F4                                   # F4 => F4
                    V0.C.F4.V.append(V0)
                    V0.C.F4.V.append(V1)
                    V0.C.F4.E.append(edge)
                    V0.C.F4.C.append(V0.C)
                                                                        # Same direction only
                    V0.C.F4.Ed.append(1)
                    edge.Fr = V0.C.F4                                   # Set right faces
                    V0.Fr   = V0.C.F4
                    V1.Fr   = V0.C.F4

                    V0   = V1.C.B                                       # Set new V0
                    
                    if V1.C.ch == 1:                                    # Set new Edge
                        edge = V1.C.B.E_b
                        V1   = edge.V0
                    else:                                               # Set new V1
                        edge = V1.C.B.E_a
                        V1   = edge.V1
                    
                elif V1.id == 1: ######################################## V1.id == 1
                    V1.C.F1 = V0.C.F4                                   # F1 => F4
                    V0.C.F4.V.append(V0)
                    V0.C.F4.V.append(V1)
                    V0.C.F4.E.append(edge)
                    V0.C.F4.C.append(V0.C)
                
                    if edge.V0.C.id == V0.C.id:                         # Same direction
                        V0.C.F4.Ed.append(1)
                        edge.Fr = V0.C.F4                               # Set right faces
                        V0.Fr   = V0.C.F4
                        V1.Fr   = V0.C.F4
                    
                    else:                                               # Opposite direction
                        V0.C.F4.Ed.append(-1)
                        edge.Fl = V0.C.F4                               # Set left faces
                        V0.Fl   = V0.C.F4
                        V1.Fl   = V0.C.F4
                    
                    V0   = V1.C.a                                       # Set new V0
                    edge = V1.C.a.E_b                                   # Set new Edge
                    V1   = edge.V0                                      # Set new V1
                    
                elif V1.id == 2: ######################################## V1.id == 2
                    V1.C.F2 = V0.C.F4                                   # F2 => F4
                    V0.C.F4.V.append(V0)
                    V0.C.F4.V.append(V1)
                    V0.C.F4.E.append(edge)
                    V0.C.F4.C.append(V0.C)
                                                                        # Opposite direction only
                    V0.C.F4.Ed.append(-1)
                    edge.Fl = V0.C.F4                                   # Set left faces
                    V0.Fl   = V0.C.F4
                    V1.Fl   = V0.C.F4

                    V0   = V1.C.b                                       # Set new V0
                    
                    if V1.C.ch == 1:                                    # Set new Edge
                        edge = V1.C.b.E_a
                        V1   = edge.V1
                    else:                                               # Set new V1
                        edge = V1.C.b.E_b
                        V1   = edge.V0
                    
                elif V1.id == 3: ######################################## V1.id == 3
                    V1.C.F3 = V0.C.F4                                   # F3 => F4
                    V0.C.F4.V.append(V0)
                    V0.C.F4.V.append(V1)
                    V0.C.F4.E.append(edge)
                    V0.C.F4.C.append(V0.C)
                
                    if edge.V0.C.id == V0.C.id:                         # Same direction
                        V0.C.F4.Ed.append(1)
                        edge.Fr = V0.C.F4                               # Set right faces
                        V0.Fr   = V0.C.F4
                        V1.Fr   = V0.C.F4
                    
                    else:                                               # Opposite direction
                        V0.C.F4.Ed.append(-1)
                        edge.Fl = V0.C.F4                               # Set left faces
                        V0.Fl   = V0.C.F4
                        V1.Fl   = V0.C.F4
                    
                    V0   = V1.C.A                                       # Set new V0
                    edge = V1.C.A.E_a                                   # Set new Edge
                    V1   = edge.V1                                      # Set new V1
                    
                else:
                    self.W += ('Value error! ' + 
                               'Vertex id has to be 0, 1, 2 or 3. ' + 
                               '"{}"'.format(V1.id) + 
                               'was given instead. ' + 
                               '"{}()" '.format(inspect.stack()[0][3]) + 
                               '- {}.\n'.format(datetime.now()))
                    self.w = True
                    notsame = False
                
            else:
                self.W += ('Value error! ' + 
                           'Vertex id has to be 0, 1, 2 or 3. ' + 
                           '"{}"'.format(V0.id) + 
                           'was given instead. ' + 
                           '"{}()" '.format(inspect.stack()[0][3]) + 
                           '- {}.\n'.format(datetime.now()))
                self.w = True
                notsame = False
            
            if first_vertex.C.id == V0.C.id:
                first_vertex.Fr.id = f_count
                self.G.F.append(first_vertex.Fr)
                self.G.f = f_count + 1
                
                summ = 0
                for i in first_vertex.Fr.Ed:
                    summ += i
                if summ == len(first_vertex.Fr.Ed) or summ == - len(first_vertex.Fr.Ed):
                    first_vertex.Fr.cy = True
                
                first_vertex.Fr.e = len(first_vertex.Fr.E)
                
                if self.G.n < first_vertex.Fr.e:
                    self.G.n  = first_vertex.Fr.e
                    self.G.f0 = first_vertex.Fr
                elif self.G.n == first_vertex.Fr.e:
                    if (self.G.f0.cy != True) and (first_vertex.Fr.cy == True):
                        self.G.f0 = first_vertex.Fr
                notsame = False
    
    ####################################################################
    
    def getFaceAdjacencyM(self):                                        # Generating 'fA' adj. matrix of face dual
        fAM = np.zeros((self.G.f, self.G.f), dtype=int)
        for n, edge in enumerate(self.G.E):
            if edge.Fl != None and edge.Fr != None:
                fAM[edge.Fl.id][edge.Fr.id] += 1
                fAM[edge.Fr.id][edge.Fl.id] += 1
        self.G.fA = fAM
        self.getFaceDistances()
        
    ####################################################################
    
    def getFaceDistances(self):                                         # Generating distances of faces from F0
        self.G.f0.df0 = 0
        fAMns = []
        fAMns_sum = 0
        any_zero = True
        n = 0
        
        while any_zero:
            fAMns.append(np.linalg.matrix_power(self.G.fA, n))
            fAMns_sum = sum(fAMns)
            any_zero = np.any(fAMns_sum == 0)
            n += 1
        
        n = 1
        for F in self.G.F:
            while F.df0 == None:
                if fAMns[n][self.G.f0.id][F.id] > 0:
                    F.df0 = n
                n += 1
            n = 1
        
        fAM_d = {}
        fAM_d[0] = [self.G.f0]
        
        for F in self.G.F:
            if F.df0 in fAM_d:
                fAM_d[F.df0].append(F)
            else:
                fAM_d[F.df0] = [F]
        
        for n, d in enumerate(list(reversed(sorted(fAM_d.keys())))):
            for F in fAM_d[d]:
                for E in F.E:
                    if E.LX['curve'] == None:
                        E.LX['curve'] = n * 10 + 5
        
        for C in self.G.C:
            C.E_a.LX['curve'] = (C.E_a.V0.E_b.LX['curve'] + 
                                 C.E_a.V1.E_a.LX['curve']) / 2
            C.E_b.LX['curve'] = (C.E_b.V0.E_b.LX['curve'] + 
                                 C.E_b.V1.E_a.LX['curve']) / 2
            
    ####################################################################
    
    def getGVE(self, rho = 100, phi0 = 0, rho_0 = None, phi_0 = None):
        self.Gr = gt.Graph()
        pin = self.Gr.new_vertex_property("boolean")                    # Creating vertex properties
        v_label = self.Gr.new_vertex_property("string")
        
        for i in range(self.G.c[0]):
            for j in range(4):
                v = self.Gr.add_vertex()
                if j == 0:
                    string = str(i + 1) + '₁' + '1'
                elif j == 1:
                    string = str(i + 1) + '₂' + '2'
                elif j == 2:
                    string = str(i + 1) + '₃' + '3'
                elif j == 3:
                    string = str(i + 1) + '₄' + '4'
                v_label[self.Gr.vertex(i * 4 + j)] = string             # Labeling vertices
        
        n_count = 0
        
        for n, E in enumerate(self.G.E):
            if E.V0.C.id != E.V1.C.id:
                v = self.Gr.add_vertex()
                string = 'E' + str(n_count + 1)
                v_label[self.G.c[0] * 4 + n_count] = string
                E.LX['id'] = self.G.c[0] * 4 + n_count                  # ID label in graph-tool object
                n_count += 1
                
            else:
                if not ((E.V0.id == 0 and E.V1.id == 2) or
                        (E.V0.id == 1 and E.V1.id == 3) or
                        (E.V0.id == 2 and E.V1.id == 0) or
                        (E.V0.id == 3 and E.V1.id == 1)):               # Including R1 (Reidemeister 1) case
                    v = self.Gr.add_vertex()
                    string = 'E' + str(n_count + 1)
                    v_label[self.G.c[0] * 4 + n_count] = string
                    E.LX['id'] = self.G.c[0] * 4 + n_count              # ID label in graph-tool object
                    n_count += 1
        
        n_count = 0
        
        for n, E in enumerate(self.G.E):
            if E.V0.C.id != E.V1.C.id:
                e = self.Gr.add_edge(int(E.V0.C.id) * 4 + 
                                     int(E.V0.id),
                                     self.G.c[0] * 4 + n_count)         # Adding between corssings edges V0
                e = self.Gr.add_edge(self.G.c[0] * 4 + n_count, 
                                     int(E.V1.C.id) * 4 + 
                                     int(E.V1.id))                      # Adding between corssings edges V1
                n_count += 1
                
            else:
                if not ((E.V0.id == 0 and E.V1.id == 2) or
                        (E.V0.id == 1 and E.V1.id == 3) or
                        (E.V0.id == 2 and E.V1.id == 0) or
                        (E.V0.id == 3 and E.V1.id == 1)):               # Including R1 (Reidemeister 1) case

                    e = self.Gr.add_edge(int(E.V0.C.id) * 4 + 
                                         int(E.V0.id),
                                         self.G.c[0] * 4 + n_count)     # Adding between corssings edges V0
                    e = self.Gr.add_edge(self.G.c[0] * 4 + n_count, 
                                         int(E.V1.C.id) * 4 + 
                                         int(E.V1.id))                  # Adding between corssings edges V1
                    n_count += 1
                    
        for n, C in enumerate(self.G.C):
            e = self.Gr.add_edge(int(C.id) * 4 + 0,
                                 int(C.id) * 4 + 1)
            e = self.Gr.add_edge(int(C.id) * 4 + 1,
                                 int(C.id) * 4 + 2)
            e = self.Gr.add_edge(int(C.id) * 4 + 2,
                                 int(C.id) * 4 + 3)
            e = self.Gr.add_edge(int(C.id) * 4 + 3,
                                 int(C.id) * 4 + 0)                     # Adding within crossings edges
        
        pos = gt.graph_tool.draw.sfdp_layout(self.Gr)                   # Initial positions
        
        if phi_0 == None:
            phi_0 = np.zeros(self.c[0] * 2, dtype=int)
        if rho_0 == None:
            rho_0 = np.zeros(self.c[0] * 2, dtype=int)
        
        angle = np.pi / len(self.G.f0.C)                                # Originally 2Pi / 2*nCr
        alternate = - angle / 2
        
        if phi0 != 0:
            phi = - angle + (2 * np.pi / phi0)
        else:
            phi = - angle
        
        for n, Phi in enumerate(phi_0):
            if Phi != 0:
                phi_0[n] = angle / Phi
        
        for n, v in enumerate(self.G.f0.V):
            pin[self.Gr.vertex(v.C.id * 4 + v.id)] = 1
            pos[v.C.id * 4 + v.id][0] = ((rho + rho_0[n]) * 
                                         np.cos(phi + phi_0[n]))
            pos[v.C.id * 4 + v.id][1] = ((rho + rho_0[n]) * 
                                         -np.sin(phi + phi_0[n]))
            
            alternate = -alternate
            phi += angle + alternate            

        pos = gt.graph_tool.draw.sfdp_layout(self.Gr, 
                                             pos = pos, 
                                             pin = pin,
                                             K = 20)                    # Creating positions
        
        gt.graph_draw(self.Gr,
                      pos = pos,
                      vertex_text = v_label, 
                      vertex_font_size=10, 
                      output_size=((self.g.c[0]*4)*10+500, 
                                   (self.g.c[0]*4)*10+500),
                      output = 'output_G.png')                          # Drawing graph
        
        self.Lpos = pos
    
    ####################################################################
        
    def genLatex(self, lw = 1.0, color = 'black', curvature = None, 
                 scale = 20, vx = None, vy = None, radii = None, 
                 vc = True, vs = 3):
        pos = self.Lpos
        
        if vx == None:
            vx = np.zeros(self.c[0] * 4, dtype=int)
        if vy == None:
            vy = np.zeros(self.c[0] * 4, dtype=int)
        if radii == None:
            radii = np.zeros(self.c[0], dtype=int)
            
        hlw = str(lw/2) + 'mm'                                          # hlw - half line width
        lw = str(lw) + 'mm'                                             # lw - line width
        
        head  = ('\\documentclass[12pt]{article}\n' +
                 '\\usepackage{amsmath}\n' + 
                 '\\usepackage{tikz}\n' +
                 '\\usepackage{bm}\n' +
                 '\\title{' + '{}'.format(self.name.split('.')[0]) + '}\n\n' + 
                 '\\begin{document}\n\n' + 
                 '\\begin{tikzpicture}\n')

        coord = ''
        lines = ''
        lines_ = ''
                    
        for n, V in enumerate(self.G.V):
            px = str((pos[V.C.id * 4 + V.id][0] / scale) + 
                     vx[V.C.id * 4 + V.id])                             # x
            py = str((pos[V.C.id * 4 + V.id][1] / scale) + 
                     vy[V.C.id * 4 + V.id])                             # y
            name = 'V_' + str(V.C.id + 1) + "_" + str(V.id + 1)         # V_n_m
            coord += ('  \\coordinate ({}) at '.format(name) + 
                      '({},{},{});\n'.format(px, py, 0))
            
            if V.id == 0:
                lines_ += ('  \\draw[line width={},'.format(hlw) + 
                           'color={},opacity=1,fill=white] '.format(color) + 
                           '({}) node [color=black] '.format(name) +
                           '{$\\bm{' + str(V.C.id + 1) + '_a}$} circle ' + 
                           '({}mm);\n'.format(vs))
            elif V.id == 1:
                if V.C.ch == 1:
                    lines_ += ('  \\draw[line width={},'.format(hlw) + 
                               'color={},opacity=1,fill={}] '.format(color, color) + 
                               '({}) node [color=white] '.format(name) + 
                               '{$\\bm{' + str(V.C.id + 1) + '_b}$} circle ' + 
                               '({}mm);\n'.format(vs))
                else:
                    lines_ += ('  \\draw[line width={},'.format(hlw) + 
                               'color={},opacity=1,fill=white] '.format(color) + 
                               '({}) node [color=black] '.format(name) +
                               '{$\\bm{' + str(V.C.id + 1) + '_b}$} circle ' + 
                               '({}mm);\n'.format(vs))
            elif V.id == 2:
                lines_ += ('  \\draw[line width={},'.format(hlw) + 
                           'color={},opacity=1,fill={}] '.format(color, color) + 
                           '({}) node [color=white] '.format(name) + 
                           '{$\\bm{' + str(V.C.id + 1) + '_A}$} circle ' + 
                           '({}mm);\n'.format(vs))
            elif V.id == 3:
                if V.C.ch == -1:
                    lines_ += ('  \\draw[line width={},'.format(hlw) + 
                               'color={},opacity=1,fill={}] '.format(color, color) + 
                               '({}) node [color=white] '.format(name) + 
                               '{$\\bm{' + str(V.C.id + 1) + '_B}$} circle ' + 
                               '({}mm);\n'.format(vs))
                else:
                    lines_ += ('  \\draw[line width={},'.format(hlw) + 
                               'color={},opacity=1,fill=white] '.format(color) + 
                               '({}) node [color=black] '.format(name) +
                               '{$\\bm{' + str(V.C.id + 1) + '_B}$} circle ' + 
                               '({}mm);\n'.format(vs))
            else:
                print('V.id == {}').format(V.id) #################################
        
        for n, E in enumerate(self.G.E):
            v0 = 'V_' + str(E.V0.C.id + 1) + '_' + str(E.V0.id + 1)
            v1 = 'V_' + str(E.V1.C.id + 1) + '_' + str(E.V1.id + 1)

            if curvature != None:
                bend = 'left=' + str(curvature[n])
            else:
                bend = 'left=' + str(E.LX['curve'])

            if E.V0.C.id != E.V1.C.id:
                lines += ('  \\draw[-,' + 
                          'line width={},'.format(lw) + 
                          'color={}] '.format(color) + 
                          '({}) to '.format(v0) + 
                          '[bend {}] '.format(bend) + 
                          '({});\n').format(v1)                         # Adding between crossings edges
            else:
                if not ((E.V0.id == 0 and E.V1.id == 2) or
                        (E.V0.id == 1 and E.V1.id == 3) or
                        (E.V0.id == 2 and E.V1.id == 0) or
                        (E.V0.id == 3 and E.V1.id == 1)):               # Including R1 (Reidemeister 1) case

                    lines += ('  \\draw[-,' + 
                              'line width={},'.format(lw) + 
                              'color={}] '.format(color) + 
                              '({}) to '.format(v0) + 
                              '[bend {}] '.format(bend) + 
                              '({});\n').format(v1)

        for n, C in enumerate(self.G.C):
            if curvature != None:
                bend_a = 'left=' + str(curvature[C.E_a.id])
                bend_b = 'left=' + str(curvature[C.E_b.id])
                
            else:
                bend_a = 'left=' + str(C.E_a.LX['curve'])
                bend_b = 'left=' + str(C.E_b.LX['curve'])

            va = 'V_' + str(C.id + 1) + '_1'
            vb = 'V_' + str(C.id + 1) + '_2'
            vA = 'V_' + str(C.id + 1) + '_3'
            vB = 'V_' + str(C.id + 1) + '_4'
            
            Cpx = ((pos[C.id * 4 + 0][0] / scale) + 
                   vx[C.id * 4 + 0] + 
                   (pos[C.id * 4 + 1][0] / scale) + 
                   vx[C.id * 4 + 1] + 
                   (pos[C.id * 4 + 2][0] / scale) + 
                   vx[C.id * 4 + 2] + 
                   (pos[C.id * 4 + 3][0] / scale) + 
                   vx[C.id * 4 + 3]) / 4
            
            Cpy = ((pos[C.id * 4 + 0][1] / scale) + 
                   vy[C.id * 4 + 0] + 
                   (pos[C.id * 4 + 1][1] / scale) + 
                   vy[C.id * 4 + 1] + 
                   (pos[C.id * 4 + 2][1] / scale) + 
                   vy[C.id * 4 + 2] + 
                   (pos[C.id * 4 + 3][1] / scale) + 
                   vy[C.id * 4 + 3]) / 4

            lines += ('  \\draw[->,' + 
                      'line width={},'.format(lw) + 
                      'color={}] '.format(color) + 
                      '({}) to '.format(va) + 
                      '[bend {}] '.format(bend_a) + 
                      '({});\n').format(vA)

            if C.ch == 1:
                lines += ('  \\draw[dotted, ->,' + 
                          'line width={},'.format(lw) + 
                          'color={}, opacity=0.25] '.format(color) + 
                          '({}) to '.format(vB) + 
                          '[bend {}] '.format(bend_b) + 
                          '({});\n').format(vb)
            elif C.ch == -1:
                lines += ('  \\draw[dotted, ->,' + 
                          'line width={},'.format(lw) + 
                          'color={}, opacity=0.25] '.format(color) + 
                          '({}) to '.format(vb) + 
                          '[bend {}] '.format(bend_b) + 
                          '({});\n').format(vB)
            else:
                print('chirality problem!!!')###################################################################

        tail  = ('\\end{tikzpicture}\n\n' +
                 '\\end{document}\n')
        
        if vc:
            lines += lines_
        
        text = head + coord + lines + tail
        self.L = text
    
    ####################################################################
########################################################################

