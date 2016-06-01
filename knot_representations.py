#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
from codecs import open as uopen
import codecs
import numpy
import re
import sympy
import scipy
import matplotlib.pyplot as pyplot
#import shapely
from scipy.sparse import csc_matrix
from graph_tool.all import *
from copy import deepcopy
import math

first_arg = sys.argv[1]
fn_in = open(first_arg, "r")
sys.stderr.write("Working with input file: {}\n".format(first_arg))
fn = []
for line in fn_in.readlines():
    fn.append(line)
filename = first_arg.rpartition('.')[0]
filename_out = filename + ".out"
fn_out = uopen(filename_out, encoding='utf-8', mode='w+')
#filename_all_rep = filename + "_all.rep"
#fn_all_rep = uopen(filename_all_rep, encoding='utf-8', mode='w+')
filename_png = filename + ".png"
filename_2png = filename + "__2.png"
filename_3png = filename + "__3.png"
filename_4png = filename + "__4.png"
filename_5png = filename + "__5.png"
filename_6png = filename + "__6.png"
filename_7png = filename + "__7.png"
filename_8png = filename + "__8.png"
filename_9png = filename + "__9.png"

A = sympy.Symbol('A')
B = sympy.Symbol('B')

##################################################################################################################
#Start of functions

##########################################################################################################
################################################################################ RemoveEmptyLines
def RemoveEmptyLines(fn_in):#                                               #### Takes -fn_in-
    '''Removing empty lines from inputfile; '''
    re_blankLine = re.compile(r"^\s*\n?$")
    in_put = []
    infinite = len(fn)
    for line in range(infinite):
        string = str(fn[line])
        if re_blankLine.match(string):
            continue
        else:
            in_put.append(string)
    return in_put############################################################### Returns -in_put-
##########################################################################################################


##########################################################################################################
################################################################################ NrOfRopes
def NrOfRopes(First_line):#                                                 #### Takes -First_line-
    '''Defining the number of the ropes from first line; '''
    Nr_of_ropes = int(len(First_line) / 2)
    return Nr_of_ropes########################################################## Returns -Nr_of_ropes-
##########################################################################################################


##########################################################################################################
################################################################################ StartingPoints
def StartingPoints(First_line, Nr_of_ropes):#                               #### Takes -First_line-, -Nr_of_ropes-
    '''Creating starting points list; '''
    #Creating starting points list
    Starting_points = numpy.empty((Nr_of_ropes, 2), numpy.int32)
    Rope_nr = 1
    while Rope_nr <= Nr_of_ropes:
        Starting_points[Rope_nr -1, 0] = int(First_line[(2 * Rope_nr) -2])
        Starting_points[Rope_nr -1, 1] = int(First_line[(2 * Rope_nr) -1])
        Rope_nr += 1
    return Starting_points###################################################### Returns -Starting_Points-
##########################################################################################################


##########################################################################################################
################################################################################ SeparateFirstLine
def SeparateFirstLine(in_put):#                                             #### Takes -in_put-
    '''Separating first line; '''
    First_line = in_put[0]
    First_line = First_line.split()
    del in_put[0]
    return in_put, First_line################################################### Returnes -in_put-, -First_line-
##########################################################################################################


##########################################################################################################
################################################################################ InputToMatrix
def InputToMatrix(in_put, Nr_of_crossings):#                                #### Takes -in_put-, -Nr_of_crossings-
    '''Creating connections matrix based on inputfile; '''
    #Transforming input string into list of list of lists
    Connections_matrix = numpy.empty((Nr_of_crossings,4,2), numpy.int32)
    a = 0
    while a < Nr_of_crossings:
        b = in_put[a].split()
        Connections_matrix[a,0,0] = int(b[0])
        Connections_matrix[a,0,1] = int(b[1])
        Connections_matrix[a,1,0] = int(b[2])
        Connections_matrix[a,1,1] = int(b[3])
        Connections_matrix[a,2,0] = int(b[4])
        Connections_matrix[a,2,1] = int(b[5])
        Connections_matrix[a,3,0] = int(b[6])
        Connections_matrix[a,3,1] = int(b[7])
        a += 1

    #Creating the connection points position matrix
    Position_matrix = numpy.empty((Nr_of_crossings,4,2), numpy.int32)
    for i in range(Nr_of_crossings):
        for j in range(4):
            Position_matrix[i,j,0] = i + 1
            Position_matrix[i,j,1] = j + 1                                  #### Returns -Connections_matrix-
    return Connections_matrix, Position_matrix##################################         -Position_matrix-
##########################################################################################################


##########################################################################################################
################################################################################ ReadRopePaths
def ReadRopePaths(Nr_of_ropes, Nr_of_crossings, Starting_points, Connections_matrix):
    '''Reading rope path; '''
    #Reading the path based on the connections
    Crossing = {1:3, 2:4, 3:1, 4:2}
    Moving_pt = numpy.empty((1,2,2), numpy.int32)
    Ending_pt = numpy.empty((1,1,2), numpy.int32)
    Starting_pt = numpy.empty((1,1,2), numpy.int32)
    This_path = numpy.empty((1,2,2), numpy.int32)
    Rope_paths = {}
    for i in range(Nr_of_ropes):
        This_path = numpy.empty((1,2,2), numpy.int32)
        This_path[0,0,0] = Starting_points[i,0]
        This_path[0,0,1] = Starting_points[i,1]
        This_path[0,1,0] = Starting_points[i,0]
        This_path[0,1,1] = Crossing[Starting_points[i,1]]
        Starting_pt[0,0,0] = Starting_points[i,0]
        Starting_pt[0,0,1] = Starting_points[i,1]
        Moving_pt[0,0,0] = Starting_points[i,0]
        Moving_pt[0,0,1] = Starting_points[i,1]
        Moving_pt[0,1,0] = Starting_points[i,0]
        Moving_pt[0,1,1] = Crossing[Starting_points[i,1]]
        Ending_pt[0,0,0] = Starting_points[i,0]
        Ending_pt[0,0,1] = Crossing[Starting_points[i,1]]
        EndEQStart = (Ending_pt == Starting_pt)
        if EndEQStart.all() == True:
            whileloop = 2
        elif EndEQStart.any() == True:
            whileloop = 1
        elif EndEQStart.all() == False:
            whileloop = 0
        while whileloop < 2:
            j = Moving_pt[0,1,0]
            k = Moving_pt[0,1,1]
            Moving_pt[0,0,0] = Connections_matrix[j-1,k-1,0]
            Moving_pt[0,0,1] = Connections_matrix[j-1,k-1,1]
            Moving_pt[0,1,0] = Moving_pt[0,0,0]
            Moving_pt[0,1,1] = Crossing[Moving_pt[0,0,1]]
            Ending_pt[0,0,0] = Moving_pt[0,0,0]
            Ending_pt[0,0,1] = Moving_pt[0,0,1]
            This_path = numpy.append(This_path, Moving_pt, axis=1)
            EndEQStart = (Ending_pt == Starting_pt)
            if EndEQStart.all() == True:
                whileloop = 2
            elif EndEQStart.any() == True:
                whileloop = 1
            elif EndEQStart.all() == False:
                whileloop = 0
        This_path = numpy.delete(This_path, -1,1)
        This_path = numpy.delete(This_path, -1,1)
        Rope_paths[i] = This_path
    return Rope_paths########################################################### Returns -Rope_paths-
##########################################################################################################


##########################################################################################################
################################################################################ ConnectionsDict
def ConnectionsDict(Nr_of_ropes, Rope_paths):#                              #### Takes -Nr_of_ropes-, -Rope_paths-
    '''Creating connections dictionary; '''
    Connections_dict = {}
    Rope_path_length = {}

    #First we need the length of the rope paths (how many connection points are involved from start to end)
    for i in range(Nr_of_ropes):
        j = Rope_paths[i]
        k = list(j.shape)
        Rope_path_length[i] = k[1]

    #Now we can generate a dictionary for the connections of the connection points (based on rope paths)
    for i in range(Nr_of_ropes):
        for j in range(Rope_path_length[i]):
            k = Rope_paths[i]
            a = k[0,j,0]
            b = k[0,j,1]
            c = [a, b]
            if j == (Rope_path_length[i] - 1):
                m = k[0,0,0]
                n = k[0,0,1]
                p = [m, n]
                Connections_dict[tuple(c)] = tuple(p)
            elif j < (Rope_path_length[i] - 1):
                e = k[0,j+1,0]
                f = k[0,j+1,1]
                g = [e, f]
                Connections_dict[tuple(c)] = tuple(g)#                      #### Returns -Rope_path_length-
    return Rope_path_length, Connections_dict###################################         -Connections_dict-
##########################################################################################################


##########################################################################################################
################################################################################ CrossingsChirality
def CrossingsChirality(Nr_of_crossings, Connections_dict):#                 #### Takes -Nr_of_crossings-, -Connections_dict-
    '''Creating crossings chirality; '''
    #Determine the chirality of the crossings (Clockwise (-1) or anticlockwise (1))
    #Now follow the path and create a dictionary for the chirality of the crossings
    Crossings_chirality_dict = {}

    for i in range(Nr_of_crossings):
        if Connections_dict[tuple([i+1, 2])] == tuple([i+1, 4]):
            if Connections_dict[tuple([i+1, 4])] != tuple([i+1, 2]):
                Crossings_chirality_dict[i] = -1
            elif Connections_dict[tuple([i+1, 4])] == tuple([i+1, 2]):
                Crossings_chirality_dict[i] = 2
        elif Connections_dict[tuple([i+1, 2])] != tuple([i+1, 4]):
            if Connections_dict[tuple([i+1, 4])] == tuple([i+1, 2]):
                Crossings_chirality_dict[i] = 1
            elif Connections_dict[tuple([i+1, 4])] != tuple([i+1, 2]):
                Crossings_chirality_dict[i] = 2

    #Now create array for the crossings chirality
    Crossings_chirality = numpy.zeros((Nr_of_crossings), numpy.int32)
    for i in range(Nr_of_crossings):
        Crossings_chirality[i] = Crossings_chirality_dict[i]#               #### Returns -Crossings_chirality-
    return Crossings_chirality, Crossings_chirality_dict########################       -Crossings_chirality_dict-
##########################################################################################################


##########################################################################################################
################################################################################ ConnectionpointsConnectivity
def ConnectionpointsConnectivity(Nr_of_crossings, Crossings_chirality, Connections_dict_inv, Connections_dict):
    '''Creating connectionpoints connectivity matrix; '''
    #Now generate the connectionpoints connectivity matrix
    row = scipy.array([])
    col = scipy.array([])
    data = scipy.array([])

    for i in range(Nr_of_crossings):
        c = Crossings_chirality[i]
        for j in range(4):
            row = scipy.append(row,[i*4 + j])
            if c == 1:
                if j == 0:
                    a = Connections_dict_inv[tuple([i+1,j+1])]
                    col = scipy.append(col, [(a[0]-1)*4 + a[1] - 1])
                    data = scipy.append(data,[-1])
                elif j == 1:
                    a = Connections_dict[tuple([i+1,j+1])]
                    col = scipy.append(col, [(a[0]-1)*4 + a[1] - 1])
                    data = scipy.append(data,[1])
                elif j == 2:
                    a = Connections_dict[tuple([i+1,j+1])]
                    col = scipy.append(col, [(a[0]-1)*4 + a[1] - 1])
                    data = scipy.append(data,[1])
                elif j == 3:
                    a = Connections_dict_inv[tuple([i+1,j+1])]
                    col = scipy.append(col, [(a[0]-1)*4 + a[1] - 1])
                    data = scipy.append(data,[-1])
            elif c == -1:
                if j == 0:
                    a = Connections_dict_inv[tuple([i+1,j+1])]
                    col = scipy.append(col, [(a[0]-1)*4 + a[1] - 1])
                    data = scipy.append(data,[-1])
                elif j == 1:
                    a = Connections_dict_inv[tuple([i+1,j+1])]
                    col = scipy.append(col, [(a[0]-1)*4 + a[1] - 1])
                    data = scipy.append(data,[-1])
                elif j == 2:
                    a = Connections_dict[tuple([i+1,j+1])]
                    col = scipy.append(col, [(a[0]-1)*4 + a[1] - 1])
                    data = scipy.append(data,[1])
                elif j == 3:
                    a = Connections_dict[tuple([i+1,j+1])]
                    col = scipy.append(col, [(a[0]-1)*4 + a[1] - 1])
                    data = scipy.append(data,[1])

    Connectionpoints_connectivity_matrix = csc_matrix( (data,(row,col)), shape=(len(row),len(col)) )#.todense()
    return Connectionpoints_connectivity_matrix################################# Returns -Connectionpoints_connectivity_matrix-
##########################################################################################################


##########################################################################################################
################################################################################ ConnectionpointsConnectivityAB
def ConnectionpointsConnectivityAB(Nr_of_crossings, Crossings_chirality, Connections_dict_inv, Connections_dict, A, B):
    '''Creating connectionpoints connectivity matrix with symbols; '''
    #Adding the u and o variables to the matrix
    CCM = sympy.SparseMatrix(Nr_of_crossings*4, Nr_of_crossings*4, {(3, 3): A})
    CCM = CCM.add(-CCM)

    for i in range(Nr_of_crossings):
        c = Crossings_chirality[i]
        for j in range(4):
            m = i*4 + j
            if c == 1:
                if j == 0:
                    a = Connections_dict_inv[tuple([i+1,j+1])]
                    n = (a[0]-1)*4 + a[1] - 1
                    CCM = CCM.add(sympy.SparseMatrix(Nr_of_crossings*4, Nr_of_crossings*4, {(m,n): -A}))
                elif j == 1:
                    a = Connections_dict[tuple([i+1,j+1])]
                    n = (a[0]-1)*4 + a[1] - 1
                    CCM = CCM.add(sympy.SparseMatrix(Nr_of_crossings*4, Nr_of_crossings*4, {(m,n): B}))
                elif j == 2:
                    a = Connections_dict[tuple([i+1,j+1])]
                    n = (a[0]-1)*4 + a[1] - 1
                    CCM = CCM.add(sympy.SparseMatrix(Nr_of_crossings*4, Nr_of_crossings*4, {(m,n): A}))
                elif j == 3:
                    a = Connections_dict_inv[tuple([i+1,j+1])]
                    n = (a[0]-1)*4 + a[1] - 1
                    CCM = CCM.add(sympy.SparseMatrix(Nr_of_crossings*4, Nr_of_crossings*4, {(m,n): -B}))
            elif c == -1:
                if j == 0:
                    a = Connections_dict_inv[tuple([i+1,j+1])]
                    n = (a[0]-1)*4 + a[1] - 1
                    CCM = CCM.add(sympy.SparseMatrix(Nr_of_crossings*4, Nr_of_crossings*4, {(m,n): -A}))
                elif j == 1:
                    a = Connections_dict_inv[tuple([i+1,j+1])]
                    n = (a[0]-1)*4 + a[1] - 1
                    CCM = CCM.add(sympy.SparseMatrix(Nr_of_crossings*4, Nr_of_crossings*4, {(m,n): -B}))
                elif j == 2:
                    a = Connections_dict[tuple([i+1,j+1])]
                    n = (a[0]-1)*4 + a[1] - 1
                    CCM = CCM.add(sympy.SparseMatrix(Nr_of_crossings*4, Nr_of_crossings*4, {(m,n): A}))
                elif j == 3:
                    a = Connections_dict[tuple([i+1,j+1])]
                    n = (a[0]-1)*4 + a[1] - 1
                    CCM = CCM.add(sympy.SparseMatrix(Nr_of_crossings*4, Nr_of_crossings*4, {(m,n): B}))
    Connectionpoints_connectivity_matrix_AB = sympy.SparseMatrix(CCM)
    return Connectionpoints_connectivity_matrix_AB############## Returns -Connectionpoints_connectivity_matrix_AB-
##########################################################################################################


##########################################################################################################
############################################################################ ConnectionpointsConnectivityABSmall
def ConnectionpointsConnectivityABSmall(Nr_of_crossings, Connectionpoints_connectivity_matrix_AB, Crossings_chirality, A, B):
    '''Creating small matrix with symbols; '''
    #Make the small matrix
    Connectionpoints_connectivity_matrix_small_AB = sympy.zeros(Nr_of_crossings, Nr_of_crossings)
    Ccms = sympy.SparseMatrix(Connectionpoints_connectivity_matrix_small_AB)

    for i in range(Nr_of_crossings):
        for j in range(Nr_of_crossings):
            m = Connectionpoints_connectivity_matrix_AB
            a = m[4*i,(4*j)] + m[4*i,(4*j+1)] + m[4*i,(4*j+2)] + m[4*i,(4*j+3)]
            b = m[(4*i+1),(4*j)] + m[(4*i+1),(4*j+1)] + m[(4*i+1),(4*j+2)] + m[(4*i+1),(4*j+3)]
            c = m[(4*i+2),(4*j)] + m[(4*i+2),(4*j+1)] + m[(4*i+2),(4*j+2)] + m[(4*i+2),(4*j+3)]
            d = m[(4*i+3),(4*j)] + m[(4*i+3),(4*j+1)] + m[(4*i+3),(4*j+2)] + m[(4*i+3),(4*j+3)]
            Ccms = Ccms.add(sympy.SparseMatrix(Nr_of_crossings, Nr_of_crossings, {(i,j): a + b + c + d}))
    for i in range(Nr_of_crossings):
        Ccms = Ccms.add(sympy.SparseMatrix(Nr_of_crossings, Nr_of_crossings, {(i,i): int(Crossings_chirality[i])}))

    Connectionpoints_connectivity_matrix_small_AB = sympy.SparseMatrix(Ccms)
    return Connectionpoints_connectivity_matrix_small_AB#### Returns -Connectionpoints_connectivity_matrix_small_AB- 
##########################################################################################################

#End of functions
##################################################################################################################
##################################################################################################################
#Start of printing functions

##########################################################################################################
################################################################################ FormatMatrix
def FormatMatrix(m, i="auto", sep=","):#                             #### Takes a -matrix- an integer and separator
    '''Pretty formatting given matrix before printing; '''
    #Making numpy pretty printing
    #numpy.set_printoptions(linewidth=140)
    j = m.shape[0]
    if i.lower() == "auto":
        i = 1
        for k in range(j):
            for l in range(j):
                if len(str(m[k, l]).replace(' ','')) > i:
                    i= len(str(m[k, l]).replace(' ',''))
    elif not type(i) == type(1):
        print "ERROR: i needs to be 'auto' or integer"
        return None
        
    s = "{\n"
    myFormat = " {:>" + str(i) + "s}"+ sep  # " {:2s}"
    for k in range(j):
        s += "{"
        for l in range(j):
            s += myFormat.format(str(m[k, l]).replace(' ',''))
        if s.endswith(sep):
            s = s[:s.rfind(sep)]
        s += "},\n"
    s += "}\n"
    return s#                                               #### Returns pretty formatted matrix as string
##########################################################################################################

##########################################################################################################
################################################################################ PrintInput
def PrintInput(Nr_of_ropes, Nr_of_crossings, Starting_points, Connections_matrix, Position_matrix, out_put):
    '''Printing input to output; '''
    #Writing to standard output file for the first time
    out_put.write("Starting points:\n\n[\n")
    for i in range(Nr_of_ropes):
        out_put.write("[{} {}] ".format(Starting_points[i,0], Starting_points[i,1]))
    out_put.write("\n]\n\nKnot representation connections read from inputfile :\n\n[\n")
    for i in range(Nr_of_crossings):
        out_put.write("[[{} {}] [".format(Connections_matrix[i,0,0], Connections_matrix[i,0,1]))
        out_put.write("{} {}] [".format(Connections_matrix[i,1,0], Connections_matrix[i,1,1]))
        out_put.write("{} {}] [".format(Connections_matrix[i,2,0], Connections_matrix[i,2,1]))
        out_put.write("{} {}]]\n".format(Connections_matrix[i,3,0], Connections_matrix[i,3,1]))
    out_put.write("]\n\nConnection points position matrix for {} crossings representation:\n\n[\n".format(Nr_of_crossings))
    for i in range(Nr_of_crossings):
        out_put.write("[[{} {}] [".format(Position_matrix[i,0,0], Position_matrix[i,0,1]))
        out_put.write("{} {}] [".format(Position_matrix[i,1,0], Position_matrix[i,1,1]))
        out_put.write("{} {}] [".format(Position_matrix[i,2,0], Position_matrix[i,2,1]))
        out_put.write("{} {}]]\n".format(Position_matrix[i,3,0], Position_matrix[i,3,1]))
    out_put.write("]\n\n")
################################################################################ No returns, print to output
##########################################################################################################


##########################################################################################################
################################################################################ PrintRopePaths
def PrintRopePaths(Nr_of_ropes, Rope_paths):#                               #### Takes -Nr_of_ropes-, -Rope_paths-
    '''Printing rope paths; '''
    #Writing out paths to the output file
    out_put.write("Rope paths:\n\n[\n")
    for i in range(Nr_of_ropes):
        j = list((Rope_paths[i]).shape)
        l = Rope_paths[i]
        out_put.write("[ ")
        for k in range(j[1]):
            out_put.write("[{} {}] ".format(l[0,k,0], l[0,k,1]))
        out_put.write("]\n")
    out_put.write("]\n\n")
################################################################################ No returns, print to output
##########################################################################################################


##########################################################################################################
################################################################################ PrintCrossingsChirality
def PrintCrossingsChirality(Crossings_chirality):#                          #### Takes -Crossings_chirality-
    '''Printing crossings chirality; '''
    out_put.write("Crossings chirality array:\n\n{}\n".format(Crossings_chirality))
################################################################################ No returns, print to output-
##########################################################################################################


##########################################################################################################
################################################################################ PrintConnectionpointsConnectivity
def PrintConnectionpointsConnectivity(Connectionpoints_connectivity_matrix):#### Takes -Connectionpoints_connectivity_matrix-
    '''Printing connectionpoints connectivity matrix; '''
    out_put.write("\nConnectionpoints connectivity matrix:\n\n{}\n".format(Connectionpoints_connectivity_matrix))
################################################################################ No returns, print to output
##########################################################################################################


##########################################################################################################
################################################################################ PrintConnectionpointsConnectivityAB
def PrintConnectionpointsConnectivityAB(Connectionpoints_connectivity_matrix_AB):
    '''Printing connectionpoints connectivity matrix with symbols; '''
    AsList = FormatMatrix(Connectionpoints_connectivity_matrix_AB)
    out_put.write("\nConnectionpoints connectivity matrix with the symbols A - B:\n\n{}".format(AsList))
################################################################################ No returns, print to output
##########################################################################################################


##########################################################################################################
################################################################################ PrintConnectionpointsConnectivityABSmall
def PrintConnectionpointsConnectivityABSmall(Connectionpoints_connectivity_matrix_small_AB):
    '''Printing small matrix with symbols; '''
    AsList = FormatMatrix(Connectionpoints_connectivity_matrix_small_AB, sep=",")
    out_put.write("\nConnections connectivity matrix:\n\n{}".format(AsList))
################################################################################ No returns, print to output
##########################################################################################################



#End of printing functions
##################################################################################################################
##################################################################################################################
#Start of knot representation generating functions



##########################################################################################################
################################################################################ DataFromInput
def DataFromInput(std_err, fn_in, out_put):#                       #### Takes -std_err-, -fn_in-, -out_put-
    '''Get data from inputfile, turn it into arrays, matrices; '''
    #First read the knot representation matrix input data from the standard input file
    #Removing empty lines, that contain just whitespace caracters and new line
    std_err.write( "{}".format(RemoveEmptyLines.__doc__))
    in_put = RemoveEmptyLines(fn_in)
    std_err.write( "{}".format(SeparateFirstLine.__doc__))
    in_put, First_line = SeparateFirstLine(in_put)
    #Defining the number of the ropes from first line
    std_err.write( "{}".format(NrOfRopes.__doc__))
    Nr_of_ropes = NrOfRopes(First_line)

    #Adding the starting point values to Starting points list
    std_err.write( "{}".format(StartingPoints.__doc__))
    Starting_points = StartingPoints(First_line, Nr_of_ropes)

    #Define the number of crossings by reading how many rows the file has
    std_err.write( "Defining the number of crossings; ")
    Nr_of_crossings = int(len(in_put))

    std_err.write( "{}".format(InputToMatrix.__doc__))
    Connections_matrix, Position_matrix = InputToMatrix(in_put, Nr_of_crossings)

    #std_err.write( "{}".format(PrintInput.__doc__))
    PrintInput(Nr_of_ropes, Nr_of_crossings, Starting_points, Connections_matrix, Position_matrix, out_put)
    return in_put, Nr_of_ropes, Starting_points, Nr_of_crossings, Connections_matrix, Position_matrix
##########################################################################################################
##########################################################################################################

##########################################################################################################
################################################################################ ConnectionpointsConnectivityMatrix
def ConnectionpointsConnectivityMatrix(std_err, Nr_of_ropes, Nr_of_crossings, Starting_points, Connections_matrix):
    '''Create the connectionpoints connectivity matrix; '''
    #Reading rope paths
    std_err.write( "{}".format(ReadRopePaths.__doc__))
    Rope_paths = ReadRopePaths(Nr_of_ropes, Nr_of_crossings, Starting_points, Connections_matrix)

    #std_err.write( "{}".format(PrintRopePaths.__doc__))
    PrintRopePaths(Nr_of_ropes, Rope_paths)
    #Crate connections dictionary
    std_err.write( "{}".format(ConnectionsDict.__doc__))
    Rope_path_length, Connections_dict = ConnectionsDict(Nr_of_ropes, Rope_paths)

    std_err.write( "{}".format(CrossingsChirality.__doc__))
    Crossings_chirality, Crossings_chirality_dict = CrossingsChirality(Nr_of_crossings, Connections_dict)

    #std_err.write( "{}".format(PrintCrossingsChirality.__doc__))
    PrintCrossingsChirality(Crossings_chirality)

    #Now generate the invers dictionary of connections
    std_err.write( "Creating the invers of connections dictionary; ")
    Connections_dict_inv = {Connections_dict[key] : key for key in Connections_dict}

    std_err.write( "{}".format(ConnectionpointsConnectivity.__doc__))
    Connectionpoints_connectivity_matrix = ConnectionpointsConnectivity(Nr_of_crossings, Crossings_chirality, Connections_dict_inv, Connections_dict)

    #std_err.write( "{}".format(PrintConnectionpointsConnectivity.__doc__))
    PrintConnectionpointsConnectivity(Connectionpoints_connectivity_matrix)
    return Rope_paths, Rope_path_length, Connections_dict, Crossings_chirality, Crossings_chirality_dict, Connections_dict_inv, Connectionpoints_connectivity_matrix
##########################################################################################################
##########################################################################################################

##########################################################################################################
################################################################################ GetPolynomial
def GetPolynomial(Nr_of_crossings, Crossings_chirality, Connections_dict_inv, Connections_dict, A, B):
    '''Adding the symbols A-bove and B-elow to the connections matrix, and calculating the determinant to get the polynomial; '''
    std_err.write( "{}".format(ConnectionpointsConnectivityAB.__doc__))
    Connectionpoints_connectivity_matrix_AB = ConnectionpointsConnectivityAB(Nr_of_crossings, Crossings_chirality, Connections_dict_inv, Connections_dict, A, B)

    #std_err.write( "{}".format(PrintConnectionpointsConnectivityAB.__doc__))
    PrintConnectionpointsConnectivityAB(Connectionpoints_connectivity_matrix_AB)

    std_err.write( "{}".format(ConnectionpointsConnectivityABSmall.__doc__))
    Connectionpoints_connectivity_matrix_small_AB = ConnectionpointsConnectivityABSmall(Nr_of_crossings, Connectionpoints_connectivity_matrix_AB, Crossings_chirality, A, B)

    #std_err.write( "{}".format(PrintConnectionpointsConnectivityABSmall.__doc__))
    PrintConnectionpointsConnectivityABSmall(Connectionpoints_connectivity_matrix_small_AB)

    #Calculate the determinant to get the polynomial
    std_err.write( "Calculating the determinant/polynomial. ")
    Determinant = sympy.SparseMatrix.det(Connectionpoints_connectivity_matrix_small_AB)

    PlusOrMinusOne = str(Determinant)
#    print(PlusOrMinusOne)
    if PlusOrMinusOne[-1] == str("1"):
        if PlusOrMinusOne[-3] == str("-"):
            Correction = 1
        elif PlusOrMinusOne[-3] == str("+"):
            Correction = -1
        else:
            print("\nSign is not + or - in determinant.\n")
    else:
        print("Determinant doesn't end with 1.")

    FactoredDet = Determinant + Correction
#    FactoredDet = sympy.powsimp(FactoredDet)
#    print(sympy.factor(FactoredDet))
    FactoredDet = FactoredDet - Correction
#    print(FactoredDet)

    #std_err.write( "Printing the determinant/polynomial.\n")
    out_put.write("\nDeterminant in SymPy format:\n\n{}\n\n".format(Determinant))

#    out_put.write("Determinant in unicode format:\n\n{}\n\n".format(sympy.pretty(Determinant, use_unicode = True)))
#    This is only ASCII, unicode doesn't work for some reason

    out_put.write("\nDeterminant in LATEX format:\n\n{}\n".format(sympy.latex(Determinant)))

    return Connectionpoints_connectivity_matrix_AB, Determinant
##########################################################################################################
##########################################################################################################

##########################################################################################################
################################################################################ GetPolynomial
def CrossingsAdjacencyMatrix(Nr_of_crossings, Nr_of_ropes, Starting_points, Rope_path_length, Crossings_chirality, Connections_dict_inv, Connections_dict):
    '''Creating the simple adjacency matrix of the crossings, '''
    std_err.write( "{}".format(CrossingsAdjacencyMatrix.__doc__))
    row = scipy.array([])
    col = scipy.array([])
    data = scipy.array([])

    for i in range(Nr_of_crossings):
        c = Crossings_chirality[i]
        for j in range(4):
            row = scipy.append(row,[i*4 + j])
            if c == 1:
                if j == 0:
                    a = Connections_dict_inv[tuple([i+1,j+1])]
                    col = scipy.append(col, [(a[0]-1)*4 + a[1] - 1])
                    data = scipy.append(data,[int(1)])
                elif j == 1:
                    a = Connections_dict[tuple([i+1,j+1])]
                    col = scipy.append(col, [(a[0]-1)*4 + a[1] - 1])
                    data = scipy.append(data,[int(1)])
                elif j == 2:
                    a = Connections_dict[tuple([i+1,j+1])]
                    col = scipy.append(col, [(a[0]-1)*4 + a[1] - 1])
                    data = scipy.append(data,[int(1)])
                elif j == 3:
                    a = Connections_dict_inv[tuple([i+1,j+1])]
                    col = scipy.append(col, [(a[0]-1)*4 + a[1] - 1])
                    data = scipy.append(data,[int(1)])
            elif c == -1:
                if j == 0:
                    a = Connections_dict_inv[tuple([i+1,j+1])]
                    col = scipy.append(col, [(a[0]-1)*4 + a[1] - 1])
                    data = scipy.append(data,[int(1)])
                elif j == 1:
                    a = Connections_dict_inv[tuple([i+1,j+1])]
                    col = scipy.append(col, [(a[0]-1)*4 + a[1] - 1])
                    data = scipy.append(data,[int(1)])
                elif j == 2:
                    a = Connections_dict[tuple([i+1,j+1])]
                    col = scipy.append(col, [(a[0]-1)*4 + a[1] - 1])
                    data = scipy.append(data,[int(1)])
                elif j == 3:
                    a = Connections_dict[tuple([i+1,j+1])]
                    col = scipy.append(col, [(a[0]-1)*4 + a[1] - 1])
                    data = scipy.append(data,[int(1)])

    Crossings_connectivity_matrix = csc_matrix( (data,(row,col)), shape=(len(row),len(col)) )

    AsList = FormatMatrix(Crossings_connectivity_matrix.astype(int), sep=",")
    out_put.write("\nCrossings connectivity matrix:\n\n{}".format(AsList))

    Crossings_adjacency_matrix = numpy.empty((Nr_of_crossings, Nr_of_crossings), numpy.int32)
    Ccms = Crossings_adjacency_matrix

    for i in range(Nr_of_crossings):
        for j in range(Nr_of_crossings):
            m = Crossings_connectivity_matrix
            a = m[4*i,(4*j)] + m[4*i,(4*j+1)] + m[4*i,(4*j+2)] + m[4*i,(4*j+3)]
            b = m[(4*i+1),(4*j)] + m[(4*i+1),(4*j+1)] + m[(4*i+1),(4*j+2)] + m[(4*i+1),(4*j+3)]
            c = m[(4*i+2),(4*j)] + m[(4*i+2),(4*j+1)] + m[(4*i+2),(4*j+2)] + m[(4*i+2),(4*j+3)]
            d = m[(4*i+3),(4*j)] + m[(4*i+3),(4*j+1)] + m[(4*i+3),(4*j+2)] + m[(4*i+3),(4*j+3)]
            Ccms[i, j] = a + b + c + d
    
    Crossings_adjacency_matrix = Ccms

    AsList = FormatMatrix(Crossings_adjacency_matrix.astype(int), sep=",")
    out_put.write("\nCrossings adjacency matrix:\n\n{}".format(AsList))

###########################Generate faces and connectivity (regions), Alexander polynomial approach
#These fesces will be used to generate the coordinates of the vertices to draw a planar graph with graphtool
#Every crossing will divide the plane into 4 regions (Alexander polynomial approach), they will be numbered (easier than text)
#Until the numbers will be determined, every crossing will have 4 regions, called R1, R2, R3, R4 and will always be positioned as the figure shows
#
#             2                                                    1
#         R1  ||  R2                                          R4   ||   R1
#             ||                                                   ||
#       1 ========> 3                                        4 ====||===> 2
#             ||                                                   ||
#         R4  ||  R3                                               ||
#             V                                               R3   V    R2
#   -1        4                                            +1       3

    def Regioning(Nr_of_crossings, Nr_of_ropes, Crossings_chirality, Connections_dict_inv, Connections_dict, ConnectionP, RCode):
        if ConnectionP[1] == 1:
            NextConnectionP = Connections_dict_inv[ConnectionP]
            if NextConnectionP[1] == 3:
                if RCode == 1:
                    AfterNextConnectionP = tuple((NextConnectionP[0], 2))
                    RCode = 2
                elif RCode == 4:
                    AfterNextConnectionP = tuple((NextConnectionP[0], 4))
                    RCode = 3
                else:
                    print("Error : Wrong RCode in Regioning for ConnectionP, should be 1 or 4 (3 --> 1)")
            elif NextConnectionP[1] == 2:
                if Crossings_chirality[NextConnectionP[0] - 1] == 1:
                    if RCode == 1:
                        AfterNextConnectionP = tuple((NextConnectionP[0], 1))
                        RCode = 1
                    elif RCode == 4:
                        AfterNextConnectionP = tuple((NextConnectionP[0], 3))
                        RCode = 2
                    else:
                        print("Error : Wrong RCode in Regioning for ConnectionP, should be 1 or 4 (2 --> 1)")
                else:
                    print("Error : Wrong chirality or connectionpoint given, cannot be Ch = -1 and 2 --> 1")
            elif NextConnectionP[1] == 4:
                if Crossings_chirality[NextConnectionP[0] - 1] == -1:
                    if RCode == 1:
                        AfterNextConnectionP = tuple((NextConnectionP[0], 3))
                        RCode = 3
                    elif RCode == 4:
                        AfterNextConnectionP = tuple((NextConnectionP[0], 1))
                        RCode = 4
                    else:
                        print("Error : Wrong RCode in Regioning for ConnectionP, should be 1 or 4 (4 --> 1)")
                else:
                    print("Error : Wrong chirality or connectionpoint given, cannot be Ch = 1 and 4 --> 1")
        ##############################################################################################################1///////////2
        elif ConnectionP[1] == 2:
            if Crossings_chirality[ConnectionP[0] - 1] == 1:
                NextConnectionP = Connections_dict[ConnectionP]
                if NextConnectionP[1] == 1:
                    if RCode == 1:
                        AfterNextConnectionP = tuple((NextConnectionP[0], 2))
                        RCode = 1
                    elif RCode == 2:
                        AfterNextConnectionP = tuple((NextConnectionP[0], 4))
                        RCode = 4
                    else:
                        print("Error : Wrong RCode in Regioning for ConnectionP, should be 1 or 2 (2 --> 1)")
                elif NextConnectionP[1] == 2:
                    if Crossings_chirality[NextConnectionP[0] - 1] == -1:
                        if RCode == 1:
                            AfterNextConnectionP = tuple((NextConnectionP[0], 3))
                            RCode = 2
                        elif RCode == 2:
                            AfterNextConnectionP = tuple((NextConnectionP[0], 1))
                            RCode = 1
                        else:
                            print("Error : Wrong RCode in Regioning for ConnectionP, should be 1 or 2 (2 --> 2)")
                    elif Crossings_chirality[NextConnectionP[0] - 1] == 1:
                        print("Error : Wrong chirality or connectionpoint given, cannot be Ch = 1 and 2 --> 2")
                elif NextConnectionP[1] == 4:
                    if Crossings_chirality[NextConnectionP[0] - 1] == 1:
                        if RCode == 1:
                            AfterNextConnectionP = tuple((NextConnectionP[0], 1))
                            RCode = 4
                        elif RCode == 2:
                            AfterNextConnectionP = tuple((NextConnectionP[0], 3))
                            RCode = 3
                        else:
                            print("Error : Wrong RCode in Regioning for ConnectionP, should be 1 or 2 (2 --> 4)")
                    elif Crossings_chirality[NextConnectionP[0] - 1] == -1:
                        print("Error : Wrong chirality or connectionpoint given, cannot be Ch = -1 and 2 --> 4")
                else:
                    print("Error : Direction problem in Regioning, 2 --> 3 ?")
            ##################################################################Ch = +1  /////////// Ch = -1
            elif Crossings_chirality[ConnectionP[0] - 1] == -1:
                NextConnectionP = Connections_dict_inv[ConnectionP]
                if NextConnectionP[1] == 3:
                    if RCode == 1:
                        AfterNextConnectionP = tuple((NextConnectionP[0], 4))
                        RCode = 3
                    elif RCode == 2:
                        AfterNextConnectionP = tuple((NextConnectionP[0], 2))
                        RCode = 2
                    else:
                        print("Error : Wrong RCode in Regioning for ConnectionP, should be 1 or 2 (3 --> 2)")
                elif NextConnectionP[1] == 2:
                    if Crossings_chirality[NextConnectionP[0] - 1] == 1:
                        if RCode == 1:
                            AfterNextConnectionP = tuple((NextConnectionP[0], 3))
                            RCode = 2
                        elif RCode == 2:
                            AfterNextConnectionP = tuple((NextConnectionP[0], 1))
                            RCode = 1
                        else:
                            print("Error : Wrong RCode in Regioning for ConnectionP, should be 1 or 2 (2 --> 2)")
                    elif Crossings_chirality[NextConnectionP[0] - 1] == -1:
                        print("Error : Wrong chirality or connectionpoint given, cannot be Ch = -1 and 2 --> 2")
                elif NextConnectionP[1] == 4:
                    if Crossings_chirality[NextConnectionP[0] - 1] == -1:
                        if RCode == 1:
                            AfterNextConnectionP = tuple((NextConnectionP[0], 1))
                            RCode = 4
                        elif RCode == 2:
                            AfterNextConnectionP = tuple((NextConnectionP[0], 3))
                            RCode = 3
                        else:
                            print("Error : Wrong RCode in Regioning for ConnectionP, should be 1 or 2 (4 --> 2)")
                    elif Crossings_chirality[NextConnectionP[0] - 1] == 1:
                        print("Error : Wrong chirality or connectionpoint given, cannot be Ch = 1 and 4 --> 2")
                else:
                    print("Error : Direction problem in Regioning, 1 --> 2 ?")
            else:
                print("Invalid value given in Regioning for the Crossings chirality of input ConnectionP")
        ################################################################################################################2//////////3
        elif ConnectionP[1] == 3:
            NextConnectionP = Connections_dict[ConnectionP]
            if NextConnectionP[1] == 1:
                if RCode == 2:
                    AfterNextConnectionP = tuple((NextConnectionP[0], 2))
                    RCode = 1
                elif RCode == 3:
                    AfterNextConnectionP = tuple((NextConnectionP[0], 4))
                    RCode = 4
                else:
                    print("Error : Wrong RCode in Regioning for ConnectionP, should be 2 or 3 (3 --> 1)")
            elif NextConnectionP[1] == 2:
                if Crossings_chirality[NextConnectionP[0] - 1] == -1:
                    if RCode == 2:
                        AfterNextConnectionP = tuple((NextConnectionP[0], 3))
                        RCode = 2
                    elif RCode == 3:
                        AfterNextConnectionP = tuple((NextConnectionP[0], 1))
                        RCode = 1
                    else:
                        print("Error : Wrong RCode in Regioning for ConnectionP, should be 2 or 3 (3 --> 2)")
                else:
                    print("Error : Wrong chirality or connectionpoint given, cannot be Ch = 1 and 3 --> 2")
            elif NextConnectionP[1] == 4:
                if Crossings_chirality[NextConnectionP[0] - 1] == 1:
                    if RCode == 2:
                        AfterNextConnectionP = tuple((NextConnectionP[0], 1))
                        RCode = 4
                    elif RCode == 3:
                        AfterNextConnectionP = tuple((NextConnectionP[0], 3))
                        RCode = 3
                    else:
                        print("Error : Wrong RCode in Regioning for ConnectionP, should be 2 or 3 (3 --> 4)")
                else:
                    print("Error : Wrong chirality or connectionpoint given, cannot be Ch = -1 and 3 --> 4")
        ################################################################################################################3//////////4
        elif ConnectionP[1] == 4:
            if Crossings_chirality[ConnectionP[0] - 1] == 1:
                NextConnectionP = Connections_dict_inv[ConnectionP]
                if NextConnectionP[1] == 3:
                    if RCode == 3:
                        AfterNextConnectionP = tuple((NextConnectionP[0], 4))
                        RCode = 3
                    elif RCode == 4:
                        AfterNextConnectionP = tuple((NextConnectionP[0], 2))
                        RCode = 2
                    else:
                        print("Error : Wrong RCode in Regioning for ConnectionP, should be 3 or 4 (3 --> 4)")
                elif NextConnectionP[1] == 2:
                    if Crossings_chirality[NextConnectionP[0] - 1] == 1:
                        if RCode == 3:
                            AfterNextConnectionP = tuple((NextConnectionP[0], 3))
                            RCode = 2
                        elif RCode == 4:
                            AfterNextConnectionP = tuple((NextConnectionP[0], 1))
                            RCode = 1
                        else:
                            print("Error : Wrong RCode in Regioning for ConnectionP, should be 3 or 4 (2 --> 4)")
                    elif Crossings_chirality[NextConnectionP[0] - 1] == -1:
                        print("Error : Wrong chirality or connectionpoint given, cannot be Ch = -1 and 2 --> 4")
                elif NextConnectionP[1] == 4:
                    if Crossings_chirality[NextConnectionP[0] - 1] == -1:
                        if RCode == 3:
                            AfterNextConnectionP = tuple((NextConnectionP[0], 1))
                            RCode = 4
                        elif RCode == 4:
                            AfterNextConnectionP = tuple((NextConnectionP[0], 3))
                            RCode = 3
                        else:
                            print("Error : Wrong RCode in Regioning for ConnectionP, should be 3 or 4 (4 --> 4)")
                    elif Crossings_chirality[NextConnectionP[0] - 1] == 1:
                        print("Error : Wrong chirality or connectionpoint given, cannot be Ch = 1 and 4 --> 4")
                else:
                    print("Error : Direction problem in Regioning, 1 --> 4 ?")
            ##################################################################Ch = +1  /////////// Ch = -1
            elif Crossings_chirality[ConnectionP[0] - 1] == -1:
                NextConnectionP = Connections_dict[ConnectionP]
                if NextConnectionP[1] == 1:
                    if RCode == 3:
                        AfterNextConnectionP = tuple((NextConnectionP[0], 2))
                        RCode = 1
                    elif RCode == 4:
                        AfterNextConnectionP = tuple((NextConnectionP[0], 4))
                        RCode = 4
                    else:
                        print("Error : Wrong RCode in Regioning for ConnectionP, should be 3 or 4 (4 --> 1)")
                elif NextConnectionP[1] == 2:
                    if Crossings_chirality[NextConnectionP[0] - 1] == -1:
                        if RCode == 3:
                            AfterNextConnectionP = tuple((NextConnectionP[0], 3))
                            RCode = 2
                        elif RCode == 4:
                            AfterNextConnectionP = tuple((NextConnectionP[0], 1))
                            RCode = 1
                        else:
                            print("Error : Wrong RCode in Regioning for ConnectionP, should be 3 or 4 (4 --> 2)")
                    elif Crossings_chirality[NextConnectionP[0] - 1] == 1:
                        print("Error : Wrong chirality or connectionpoint given, cannot be Ch = 1 and 4 --> 2")
                elif NextConnectionP[1] == 4:
                    if Crossings_chirality[NextConnectionP[0] - 1] == 1:
                        if RCode == 3:
                            AfterNextConnectionP = tuple((NextConnectionP[0], 1))
                            RCode = 4
                        elif RCode == 4:
                            AfterNextConnectionP = tuple((NextConnectionP[0], 3))
                            RCode = 3
                        else:
                            print("Error : Wrong RCode in Regioning for ConnectionP, should be 3 or 4 (4 --> 4)")
                    elif Crossings_chirality[NextConnectionP[0] - 1] == -1:
                        print("Error : Wrong chirality or connectionpoint given, cannot be Ch = -1 and 4 --> 4")
                else:
                    print("Error : Direction problem in Regioning, 4 --> 3 ?")
            else:
                print("Invalid value given in Regioning for the Crossings chirality of input ConnectionP")
        else:
            print("Error : ConnectionP[1] must be 1, 2, 3 or 4 in Regioning")

        return NextConnectionP, AfterNextConnectionP, RCode

####################################################################################################################################
    def GetEdgesDict(Nr_of_crossings, Nr_of_ropes, Rope_paths, Rope_path_length, Connections_dict, Starting_points, Nr_of_edges):
        '''Calculating Connectionpoint_Which_Rope_dict; '''
        Connectionpoint_Which_Rope_dict = {}
        Edges_dict = {}
        Edges_list = []
        Edges_cross_dict = {}
        Edges_cross_list = []

        for i in range(Nr_of_ropes):
            a = Rope_paths[i]
            b = a[0]
            for j in range(Rope_path_length[i]):
                c = b[j]
                Connectionpoint_Which_Rope_dict[tuple(c)] = i
        
        for i in range(Nr_of_ropes):
            StartP = Connections_dict[tuple(Starting_points[i])]
            for j in range(Rope_path_length[i] / 2):
                MovingP = Connections_dict[StartP]
                c = tuple((StartP, MovingP))
                Edges_list.append(c)
                StartP = MovingP
        
        for i in range(Nr_of_edges):
            Edges_dict[i] = Edges_list[i]

        for i in range(Nr_of_ropes):
            StartP = tuple(Starting_points[i])
            for j in range(Rope_path_length[i] / 2):
                MovingP = Connections_dict[StartP]
                c = tuple((StartP, MovingP))
                Edges_cross_list.append(c)
                StartP = MovingP
        
        for i in range(Nr_of_edges):
            Edges_cross_dict[i] = Edges_cross_list[i]

        return Connectionpoint_Which_Rope_dict, Edges_dict, Edges_list, Edges_cross_list, Edges_cross_dict
####################################################################################################################################

    Connectionpoint_Which_Rope_dict, Edges_dict, Edges_list, Edges_cross_list, Edges_cross_dict = GetEdgesDict(Nr_of_crossings, Nr_of_ropes, Rope_paths, Rope_path_length, Connections_dict, Starting_points, Nr_of_edges)

    Crossing_Regions = numpy.empty((Nr_of_crossings,4), numpy.int32)
    for i in range(Nr_of_crossings):
        for j in range(4):
            Crossing_Regions[i, j] = j + 1

    Crossing_Faces = numpy.empty((Nr_of_crossings,4), numpy.int32)
    for i in range(Nr_of_crossings):
        for j in range(4):
            Crossing_Faces[i, j] = 0

#    ConnectionP = tuple((1, 1))
#    RCode = 1

#    NextConnectionP, AfterNextConnectionP, RC = Regioning(Nr_of_crossings, Nr_of_ropes, Crossings_chirality, Connections_dict_inv, Connections_dict, ConnectionP, RCode)

    V = Nr_of_crossings
    E = Nr_of_edges
    F = Nr_of_faces
    R = Nr_of_ropes
    RC = 1
    CP = tuple((1, 1))
    NC = 0 #NextConnectionP
    ANC = 0 #AfterNextConnectionP
        
    Cycles_list = []
    f = 1
    for i in range(Nr_of_crossings):
        for j in range(4):
            if Crossing_Faces[i, j] == 0:
                CP = tuple((i + 1, j + 1))
                RC = j + 1
                Cycle_list = []
                Cycle_list.append(CP)
                Crossing_Faces[i, j] = f

                NC, ANC, RC = Regioning(V, R, Crossings_chirality, Connections_dict_inv, Connections_dict, CP, RC)
                Cycle_list.append(NC)
                Cycle_list.append(ANC)
                Crossing_Faces[ANC[0] - 1, RC - 1] = f 
                while ANC != CP:
                    NC, ANC, RC = Regioning(V, R, Crossings_chirality, Connections_dict_inv, Connections_dict, ANC, RC)
                    Cycle_list.append(NC)
                    Cycle_list.append(ANC)
                    Crossing_Faces[ANC[0] - 1, RC - 1] = f
                del Cycle_list[-1]
                Cycles_list.append(Cycle_list)
                f = f + 1
    out_put.write("\nCycles in the graph, borders of the faces (In an arbitrary order):\n")
    for i in range(len(Cycles_list)):
        pr = str(Cycles_list[i])
        out_put.write("\n{}".format(pr))
    out_put.write("\n\nFace codes of the crossing regions:\n\n")
    out_put.write("{}".format(Crossing_Faces))

    def GetDualGraph(V, E, F, R, Crossing_Faces):
        '''Calculating the dual of the graph; '''
        DualGraphAdjacency = numpy.zeros((F, F), numpy.int32)
        for i in range(V):
            for j in range(4):
                if j < 3:
                    a = Crossing_Faces[i, j] - 1
                    b = Crossing_Faces[i, j + 1] - 1
                    if DualGraphAdjacency[b, a] == 1:
                        a = 0                                       #Do nothing
                    elif DualGraphAdjacency[a, b] == 1:
                        a = 0
                    else:
                        if a < b:
                            DualGraphAdjacency[a, b] = 1
                        elif b < a:
                            DualGraphAdjacency[b, a] = 1
                if j == 3:
                    a = Crossing_Faces[i, j] - 1
                    b = Crossing_Faces[i, 0] - 1
                    if DualGraphAdjacency[b, a] == 1:
                        a = 0
                    elif DualGraphAdjacency[a, b] == 1:
                        a = 0
                    else:
                        if a < b:
                            DualGraphAdjacency[a, b] = 1
                        elif b < a:
                            DualGraphAdjacency[b, a] = 1
        return DualGraphAdjacency

    DualGraphAdjacency = GetDualGraph(V, E, F, R, Crossing_Faces)
    out_put.write("\n\nAdjacency matrix of the dual graph:\n\n")
    out_put.write("{}\n\n".format(DualGraphAdjacency))

    GCrossingsDual = Graph(directed=False)
    duver = GCrossingsDual.add_vertex(Nr_of_faces)
    fprop_string = GCrossingsDual.new_vertex_property("string")
    for i in range(Nr_of_faces):
        fprop_string[GCrossingsDual.vertex(i)] = i + 1

    ed = Crossing_Faces[0, 0]
    ge = Crossing_Faces[0, 1]
    GCrossingsDual.add_edge(ed, ge) #I have to add this so it has at least one edge, otherwise it will crash

    for i in range(Nr_of_faces):
        for j in range(Nr_of_faces):
            if DualGraphAdjacency[i, j] == 1:
                exists = 0
                
                for e in GCrossingsDual.edges():
                    a = tuple((i, j))
                    b = tuple((e.source(), e.target()))
                    if a == b:
                        exists += 1
                for e in GCrossingsDual.edges():
                    a = tuple((i, j))
                    b = tuple((e.target(), e.source()))
                    if a == b:
                        exists += 1
                if exists == 0:
                    GCrossingsDual.add_edge(i, j)

    Nr_of_cycles = len(Cycles_list)     #pretty obvious, it gives the number of cycles in the crossings graph
    Max_cycle_length = 0
    for i in range(Nr_of_cycles):
        if (len(Cycles_list[i]) / 2) > Max_cycle_length:
            Max_cycle_length = (len(Cycles_list[i]) / 2)        #It will be set to the length of the longest cycle(s)
    
    Nr_of_max_length_cycles = 0
    for i in range(Nr_of_cycles):
        if (len(Cycles_list[i]) / 2) == Max_cycle_length:
            Nr_of_max_length_cycles += 1                        #This will count how many cycles we have with the maximum length

    Remove_outside_face = 1                                     #This will tell how many faces you want to remove (if zero, nothing is removed)
    Skip = 0                 #This will tell how many faces you want to skip before removing one, for the case when the first would not be the right one
    for i in range(Nr_of_faces):
        if (len(Cycles_list[i]) / 2) == Max_cycle_length:
            if Remove_outside_face > 0:
                if Skip > 0:
                    Skip -= 1
                else:
                    for j in range(Nr_of_faces):
                        continue
                        GCrossingsDual.remove_edge(GCrossingsDual.vertex(i), GCrossingsDual.vertex(j))
                        GCrossingsDual.remove_edge(GCrossingsDual.vertex(j), GCrossingsDual.vertex(i))    #This will remove the edges of the first vertex in dual that is in the middle of one of the longest cycles
                        
                    out_put.write("Removing edges of vertex '{}' from the Dual Graph for drawing purposes.\nIt is assumed that this is the outside face of the original graph.\n\n".format(i + 1))
                    Outside_face = i + 1
                    Remove_outside_face -= 1                        #but will leave the rest as is

    #graph_tool.draw.graphviz_draw(GCrossingsDual, sep = 1.0, overlap = "scalexy")          #this is good as well, if it can give us the coordinates
    #graph_draw(GCrossingsDual, vertex_text=fprop_string, vertex_font_size=60, output_size=(1500, 1500), output=filename_2png)

    def BiColorGraph(DualGraphAdjacency, F, Outside_face): ######################################################################Not done yet!!!!!!!!!!!!!!!!!
        '''Bicoloring the graph; '''
        Seifert_coloring = {}
        Seifert_coloring[Outside_face - 1] = 0
        for i in range(F):
            if Seifert_coloring[i] == 0:
                for j in range(F):
                    if DualGraphAdjacency[i, j] == 1:
                        Seifert_coloring[j] = 1
            if Seifert_coloring[i] == 1:
                for j in range(F):
                    if DualGraphAdjacency[i, j] == 1:
                        Seifert_coloring[j] = 0

        return Seifert_coloring

#    Seifert_coloring = BiColorGraph(DualGraphAdjacency, F, Outside_face)   #Has to be finished, problem, when outside face is not the first face


###########################Using graphtool to generate the graph and the image of the graph
#Connectionpoints graph

    GConnectionpoints = Graph(directed=True)                        # A directed graph is defined for the connectionpoints
    vlist = GConnectionpoints.add_vertex((4*Nr_of_crossings))       #Adding the number of vertices, it's 4* the number of crossings

    vprop_string = GConnectionpoints.new_vertex_property("string")  #Creating a property for the vertices that is a string
    for i in range(Nr_of_crossings):                                #this will be used to label the vertices
        for j in range(4):
            string1 = str(i+1)
            if j == 0:
                string2 = u'\u2081'     #unicode text is used to print subscript, for the connectionpoints
            elif j == 1:
                string2 = u'\u2082'     #"u208a" will give "a" as subscript, where "a" goes from 1-9
            elif j == 2:
                string2 = u'\u2083'
            elif j == 3:
                string2 = u'\u2084'
            else:
                print("Connectionpoint index is other than 1,2,3,4")    #just testing if evrything is ok
            string3 = string2.encode('utf-8')
            strings = "".join((string1, string3))
            vprop_string[GConnectionpoints.vertex(4*i + j)] = strings   #This line will actually add the property to that vertex
    
    for r in range(Nr_of_ropes):
        startP = tuple(Starting_points[r])
        moveP = Connections_dict[startP]                                            #Following the connections dictionary
        a = startP[0]                                                               #using the starting points in each rope
        b = startP[1]
        c = moveP[0]                                                                #a moving point is definedthat will follow the path
        d = moveP[1]
        e = GConnectionpoints.add_edge(4 * (a - 1) + b - 1, 4 * (c - 1) + d - 1)    #and add edges to the graph, first the very first one
        for i in range(Rope_path_length[r]-1):                                      #that is the starting point connected the the first moving point
            a = moveP[0]
            b = moveP[1]
            moveP = Connections_dict[moveP]                                         #than the moving point moves one step forward on the path
            c = moveP[0]
            d = moveP[1]
            e = GConnectionpoints.add_edge(4 * (a - 1) + b - 1, 4 * (c - 1) + d - 1)    #and the next edge is added, and so on

    Dense = (Connectionpoints_connectivity_matrix.todense()).tolist()
#    print Dense

    pos = graph_tool.draw.sfdp_layout(GConnectionpoints)
    for i in range(Nr_of_crossings * 4):
        for j in range(Nr_of_crossings * 4):
            if Dense[i][j] != 0.0:
                pos[i][0] = i * 100
                pos[i][1] = j * 100

#    graph_tool.draw.sfdp_layout(GConnectionpoints)
    graph_draw(GConnectionpoints, pos = pos, vertex_text=vprop_string, vertex_font_size=60, output_size=((Nr_of_crossings * 4 + 2) * 100, (Nr_of_crossings * 4 + 2) * 100), output=filename_3png)


#Crossings graph

    GCrossings = Graph(directed=False)                                   # A directed graph is defined for the crossings
    vlist = GCrossings.add_vertex(Nr_of_crossings)                      #Adding the number of vertices = the number of crossings

    vprop2_string = GCrossings.new_vertex_property("int32_t")           #since we don't need subscript the crossing number can be used to label
    for i in range(Nr_of_crossings):
        vprop2_string[GCrossings.vertex(i)] = i + 1             #we still have to create the property, because the vertices are numbered from 0 otherwise

    for r in range(Nr_of_ropes):                                #same idea as with the connectionpoints, adding the edges with starting and moving point
        startP = tuple(Starting_points[r])
        moveP = Connections_dict[startP]
        a = startP[0]
        c = moveP[0]
        if a != c:
            print("Error: Starting point is connected to a different crossings connectionpoint!")  #just a check
        for i in range(Rope_path_length[r]-1):
            a = moveP[0]
            moveP = Connections_dict[moveP]
            c = moveP[0]
            if a != c:                                          #if the crossing number of the two connectionpoints is not the same
                e = GCrossings.add_edge(a - 1, c - 1)           #the edge is added to connect those vertices (otherwise it would self connect)
    pin = GCrossings.new_vertex_property("boolean")
    for i in range(Nr_of_crossings):
        pin[GCrossings.vertex(i)] = 1

    #pos = graph_tool.draw.sfdp_layout(GCrossings, pos=pos, pin=pin)
                           #pos=pos,
#    graph_draw(GCrossings,           nodesfirst = True, vertex_text=vprop2_string, vertex_font_size=60, output_size=(2000, 1000), output=filename_png)
#    graph_tool.draw.interactive_window(GCrossings, vertex_text=vprop2_string, vertex_font_size=30)

    GCrossingSkeleton = Graph(directed=False)                                   # An undirected graph is defined for the crossings
    vlist = GCrossingSkeleton.add_vertex(Nr_of_crossings)                              #Adding the number of vertices = the number of crossings
    vprop3_string = GCrossingSkeleton.new_vertex_property("int32_t")           #since we don't need subscript the crossing number can be used to label
    for i in range(Nr_of_crossings):
        vprop3_string[GCrossingSkeleton.vertex(i)] = i + 1             #we still have to create the property, because the vertices are numbered from 0 otherwise
    edge_counter = {}
    for r in range(Nr_of_ropes):                                #same idea as with the connectionpoints, adding the edges with starting and moving point
        startP = tuple(Starting_points[r])
        moveP = Connections_dict[startP]
        a = startP[0]
        c = moveP[0]
        if a != c:
            print("Error: Starting point is connected to a different crossings connectionpoint!")  #just a check
        for i in range(Rope_path_length[r]-1):
            a = moveP[0]
            moveP = Connections_dict[moveP]
            c = moveP[0]
            if a != c:                                          #if the crossing number of the two connectionpoints is not the same
                if tuple((c - 1, a - 1)) in edge_counter:       #if the edge with the opposite direction already exists
                    continue                                    #do nothing
                if tuple((a - 1, c - 1)) in edge_counter:       #or if the edge with the already exists
                    continue                                    #do nothing
                else:
                    edge_counter[tuple((a - 1, c - 1))] = 1     #otherwise add the edge to the counter
                    e = GCrossingSkeleton.add_edge(a - 1, c - 1)       #and to the graph

    p, kur = is_planar(GCrossingSkeleton, kuratowski=True)
    GCrossingSkeleton.set_edge_filter(kur, True)
    make_maximal_planar(GCrossingSkeleton)
#    graph_draw(GCrossingSkeleton, vertex_text=vprop3_string, vertex_font_size=60, output_size=(2000, 2000), output=filename_9png)

#    AdjacencySkeleton = graph_tool.spectral.adjacency(GCrossingSkeleton)
#    ew_A, ev_A = scipy.linalg.eig(AdjacencySkeleton.todense())
#    pyplot.figure(figsize=(8, 2))
#    pyplot.scatter(numpy.real(ew_A), numpy.imag(ew_A), c=numpy.sqrt(abs(ew_A)), linewidths=0, alpha=0.6)
#    pyplot.xlabel(r"$\operatorname{Re}(\lambda)$")
#    pyplot.ylabel(r"$\operatorname{Im}(\lambda)$")
#    pyplot.tight_layout()
#    pyplot.savefig(filename_9png)

#    LaplacianSkeleton = graph_tool.spectral.laplacian(GCrossingSkeleton, deg='total', normalized=True, weight=None, index=None)
#    ew_L, ev_L = scipy.linalg.eig(LaplacianSkeleton.todense())
#    pyplot.figure(figsize=(8, 2))
#    pyplot.scatter(numpy.real(ew_L), numpy.imag(ew_L), c=numpy.sqrt(abs(ew_L)), linewidths=0, alpha=0.6)
#    pyplot.xlabel(r"$\operatorname{Re}(\lambda)$")
#    pyplot.ylabel(r"$\operatorname{Im}(\lambda)$")
#    pyplot.tight_layout()
#    pyplot.savefig(filename_8png)

#    TransitionSkeleton = graph_tool.spectral.transition(GCrossingSkeleton)
#    ew_T, ev_T = scipy.linalg.eig(TransitionSkeleton.todense())
#    pyplot.figure(figsize=(8, 2))
#    pyplot.scatter(numpy.real(ew_T), numpy.imag(ew_T), c=numpy.sqrt(abs(ew_T)), linewidths=0, alpha=0.6)
#    pyplot.xlabel(r"$\operatorname{Re}(\lambda)$")
#    pyplot.ylabel(r"$\operatorname{Im}(\lambda)$")
#    pyplot.tight_layout()
#    pyplot.savefig(filename_7png)

#    ModularitySkeleton = graph_tool.spectral.modularity_matrix(GCrossingSkeleton)
#    ModularitySkeleton = ModularitySkeleton * numpy.identity(ModularitySkeleton.shape[0])
#    ew_M, ev_M = scipy.linalg.eig(ModularitySkeleton)
#    pyplot.figure(figsize=(8, 2))
#    pyplot.scatter(numpy.real(ew_M), numpy.imag(ew_M), c=numpy.sqrt(abs(ew_M)), linewidths=0, alpha=0.6)
#    pyplot.xlabel(r"$\operatorname{Re}(\lambda)$")
#    pyplot.ylabel(r"$\operatorname{Im}(\lambda)$")
#    pyplot.tight_layout()
#    pyplot.savefig(filename_6png)

#    Ag = AdjacencySkeleton.todense()
#    NegAg = Ag * (-1)
#    Lg = NegAg
#    for i in range(Nr_of_crossings):
#        a = 0
#        for j in range(Nr_of_crossings):
#            a = a + Ag[i, j]
#        Lg[i, i] = a
    
#    Diff = Lg - LaplacianSkeleton
#    print Diff

#    ew_l, ev_l = scipy.linalg.eig(Lg)

#    pos = graph_tool.draw.sfdp_layout(GCrossings)
#    for i in range(Nr_of_crossings):
#        pos[i][0] = ev_L[1][i]
#        pos[i][1] = ev_L[2][i]

#    BuildUpGraph = Graph(directed=False)
#    vlist = BuildUpGraph.add_vertex(Nr_of_crossings)
#    vprop3_string = GCrossingSkeleton.new_vertex_property("int32_t")
#    for i in range(Nr_of_crossings):
#        vprop3_string[GCrossingSkeleton.vertex(i)] = i + 1

#    StepByStep = Graph(directed=False)
#    vlist = StepByStep.add_vertex(Nr_of_crossings)
#    vertex_name = StepByStep.new_vertex_property("int32_t")
#    for i in range(Nr_of_crossings):
#        vertex_name[StepByStep.vertex(i)] = i + 1
#    lengths = []
#    for array in Cycles_list:
#        length = len(array)
#        lengths.append(length)

#    maxlength = 0
#    for length in lengths:
#        if length > maxlength:
#            maxlength = length

#    maxlength_count = 0
#    longest_lengths = []
#    for i in range(Nr_of_cycles):
#        if lengths[i] == maxlength:
#            maxlength_count += 1
#            longest_lengths.append(i)
    
#    possible_cycles = []

#    for i in range(maxlength_count):
#        Cycles_list_copy = deepcopy(Cycles_list)
#        Cycles_list_copy.pop(longest_lengths[i])
#        possible_cycles.append(Cycles_list_copy)
    
#    out_put.write("\nPossible cycles after one face removed as the outside face:\n")
#    for i in range(maxlength_count):
#        out_put.write("\n{}\n".format(possible_cycles[i]))

#    edges = []
#    empty = []
#    groups = []
#    for i in range(Nr_of_cycles - 1):
#        edges.append(empty)

#    for i in range(maxlength_count):
#        for j in range(Nr_of_cycles - 1):
#            cycles_length = len(possible_cycles[i][j])
#            for k in range(cycles_length):
#                if k < cycles_length - 1:
#                    if possible_cycles[i][j][k][0] == possible_cycles[i][j][k + 1][0]:
#                        continue
#                    else:
#                        break_it = 0
#                        for l in range(j + 1):
#                            if tuple((possible_cycles[i][j][k][0] - 1, possible_cycles[i][j][k + 1][0] - 1)) in edges[l]:
#                                break_it = 1
#                                continue
#                            elif tuple((possible_cycles[i][j][k + 1][0] - 1, possible_cycles[i][j][k][0] - 1)) in edges[l]:
#                                break_it = 1
#                                continue
#                        if break_it == 0:
#                            if tuple((possible_cycles[i][j][k][0] - 1, possible_cycles[i][j][k + 1][0] - 1)) in edges[j]:
#                                break_it = 1
#                                continue
#                            elif tuple((possible_cycles[i][j][k + 1][0] - 1, possible_cycles[i][j][k][0] - 1)) in edges[j]:
#                                break_it = 1
#                                continue
#                            else:
#                                if break_it == 0:
#                                    edges[j].append(tuple((possible_cycles[i][j][k][0] - 1, possible_cycles[i][j][k + 1][0] - 1)))
#                                    groups.append(tuple((j, possible_cycles[i][j][k][0] - 1, possible_cycles[i][j][k + 1][0] - 1)))
#                else:
#                    if possible_cycles[i][j][k][0] == possible_cycles[i][j][0][0]:
#                        continue
#                    else:
#                        break_it = 0
#                        for l in range(j + 1):
#                            if tuple((possible_cycles[i][j][k][0] - 1, possible_cycles[i][j][0][0] - 1)) in edges[l]:
#                                break_it = 1
#                                continue
#                            elif tuple((possible_cycles[i][j][0][0] - 1, possible_cycles[i][j][k][0] - 1)) in edges[l]:
#                                break_it = 1
#                                continue
#                        if break_it == 0:
#                            if tuple((possible_cycles[i][j][k][0] - 1, possible_cycles[i][j][0][0] - 1)) in edges[j]:
#                                break_it = 1
#                                continue
#                            elif tuple((possible_cycles[i][j][0][0] - 1, possible_cycles[i][j][k][0] - 1)) in edges[j]:
#                                break_it = 1
#                                continue
#                            else:
#                                if break_it == 0:
#                                    edges[j].append(tuple((possible_cycles[i][j][k][0] - 1, possible_cycles[i][j][0][0] - 1)))
#                                    groups.append(tuple((j, possible_cycles[i][j][k][0] - 1, possible_cycles[i][j][0][0] - 1)))
    
#    pin = StepByStep.new_vertex_property("boolean")
#    for i in range(Nr_of_crossings):
#        pin[StepByStep.vertex(i)] = 0

#    for i in range(Nr_of_cycles - 1):
#        for edge in groups:
#            if edge[0] == i:
#                StepByStep.add_edge(edge[1], edge[2])
#                pin[StepByStep.vertex(edge[1])] = 1
#                pin[StepByStep.vertex(edge[1])] = 1
#        pos = graph_tool.draw.sfdp_layout(StepByStep, pin = None)
#        filename_cycle_png = filename + "_" + str(i) + ".png"
#        graph_draw(StepByStep, pos = pos, pin = pin, vertex_text=vertex_name, vertex_font_size=60, output_size=(2000, 2000), output=filename_cycle_png)

#    graph_draw(StepByStep, vertex_text=vertex_name, vertex_font_size=60, output_size=(2000, 2000), output=filename_png)

    Graph_of_the_knot = Graph(directed=True)
    vlist = Graph_of_the_knot.add_vertex(Nr_of_crossings * 4)
    vertex_name = Graph_of_the_knot.new_vertex_property("string")
    for i in range(Nr_of_crossings):
        for j in range(4):
            string1 = str(i+1)
            if j == 0:
                string2 = u'\u2081'         #unicode text is used to print subscript, for the connectionpoints
            elif j == 1:
                string2 = u'\u2082'         #"u208a" will give "a" as subscript, where "a" goes from 1-9
            elif j == 2:
                string2 = u'\u2083'         #for example "u2083" "gives underscore 3" can be printed in the output text too
            elif j == 3:
                string2 = u'\u2084'
            else:
                print("Connectionpoint index is other than 1,2,3,4")        #just a check
            string3 = string2.encode('utf-8')
            strings = "".join((string1, string3))
            vertex_name[Graph_of_the_knot.vertex(4*i + j)] = strings      #this line will add the property to that vertex, the name

    for i in range(Nr_of_crossings):
        e = Graph_of_the_knot.add_edge(4*i + 0, 4*i + 1)
        e = Graph_of_the_knot.add_edge(4*i + 1, 4*i + 2)
        e = Graph_of_the_knot.add_edge(4*i + 2, 4*i + 3)
        e = Graph_of_the_knot.add_edge(4*i + 3, 4*i + 0)

    lengths = []
    for array in Cycles_list:
        length = len(array)
        lengths.append(length)

    maxlength = 0
    for length in lengths:
        if length > maxlength:
            maxlength = length

    maxlength_count = 0
    longest_lengths = []
    for i in range(Nr_of_cycles):
        if lengths[i] == maxlength:
            maxlength_count += 1
            longest_lengths.append(i)
    
    possible_cycles = []

    for i in range(maxlength_count):
        Cycles_list_copy = deepcopy(Cycles_list)
        Cycles_list_copy.pop(longest_lengths[i])
        possible_cycles.append(Cycles_list_copy)
    
    out_put.write("\nPossible cycles after one face removed as the outside face:\n")
    for i in range(maxlength_count):
        out_put.write("\n{}\n".format(possible_cycles[i]))

    empty_list = []

    pin = Graph_of_the_knot.new_vertex_property("boolean")
    for i in range(4*Nr_of_crossings):
        pin[Graph_of_the_knot.vertex(i)] = 0

    for i in range(1):
        possible_edges = []
        possible_edges_dict = {}
        for cycle in possible_cycles[i]:
            edges_of_this_cycle = []
            for j in range(len(cycle)):
                if j != len(cycle) - 1:
                    if cycle[j][0] != cycle[j + 1][0]:
                        if Connections_dict[tuple((cycle[j][0], cycle[j][1]))] == tuple((cycle[j + 1][0], cycle[j + 1][1])):
                            if tuple((4*(cycle[j][0] - 1) + cycle[j][1] - 1, 4*(cycle[j+1][0] - 1) + cycle[j+1][1] - 1)) in possible_edges_dict:
                                continue
                            else:    
                                edges_of_this_cycle.append([4*(cycle[j][0] - 1) + cycle[j][1] - 1, 4*(cycle[j+1][0] - 1) + cycle[j+1][1] - 1])
                                possible_edges_dict[tuple((4*(cycle[j][0] - 1) + cycle[j][1] - 1, 4*(cycle[j+1][0] - 1) + cycle[j+1][1] - 1))] = 1
                        elif Connections_dict[tuple((cycle[j+1][0], cycle[j+1][1]))] == tuple((cycle[j][0], cycle[j][1])):
                            if tuple((4*(cycle[j+1][0] - 1) + cycle[j+1][1] - 1, 4*(cycle[j][0] - 1) + cycle[j][1] - 1)) in possible_edges_dict:
                                continue
                            else:
                                edges_of_this_cycle.append([4*(cycle[j+1][0] - 1) + cycle[j+1][1] - 1, 4*(cycle[j][0] - 1) + cycle[j][1] - 1])
                                possible_edges_dict[tuple((4*(cycle[j+1][0] - 1) + cycle[j+1][1] - 1, 4*(cycle[j][0] - 1) + cycle[j][1] - 1))] = 1
                elif j == len(cycle) -1:
                    if cycle[j][0] != cycle[0][0]:
                        if Connections_dict[tuple((cycle[j][0], cycle[j][1]))] == tuple((cycle[0][0], cycle[0][1])):
                            if tuple((4*(cycle[j][0] - 1) + cycle[j][1] - 1, 4*(cycle[0][0] - 1) + cycle[0][1] - 1)) in possible_edges_dict:
                                continue
                            else:
                                edges_of_this_cycle.append([4*(cycle[j][0] - 1) + cycle[j][1] - 1, 4*(cycle[0][0] - 1) + cycle[0][1] - 1])
                                possible_edges_dict[tuple((4*(cycle[j][0] - 1) + cycle[j][1] - 1, 4*(cycle[0][0] - 1) + cycle[0][1] - 1))] = 1
                        elif Connections_dict[tuple((cycle[0][0], cycle[0][1]))] == tuple((cycle[j][0], cycle[j][1])):
                            if tuple((4*(cycle[0][0] - 1) + cycle[0][1] - 1, 4*(cycle[j][0] - 1) + cycle[j][1] - 1)) in possible_edges_dict:
                                continue
                            else:
                                edges_of_this_cycle.append([4*(cycle[0][0] - 1) + cycle[0][1] - 1, 4*(cycle[j][0] - 1) + cycle[j][1] - 1])
                                possible_edges_dict[tuple((4*(cycle[0][0] - 1) + cycle[0][1] - 1, 4*(cycle[j][0] - 1) + cycle[j][1] - 1))] = 1
            if edges_of_this_cycle != empty_list:
                possible_edges.append(edges_of_this_cycle)
        for l in range(len(possible_edges)):
            for edge in possible_edges[l]:
                e = Graph_of_the_knot.add_edge(edge[0], edge[1])
                pin[Graph_of_the_knot.vertex(edge[0])] = 1
                pin[Graph_of_the_knot.vertex(edge[1])] = 1
            pos = graph_tool.draw.sfdp_layout(Graph_of_the_knot, pin = None)
            filename_pos_cyc_png = filename + "_" + str(l) + "_" + str(i) + ".png"
    ADJ = graph_tool.spectral.adjacency(Graph_of_the_knot)

    Crossings_set = set()

    for i in range(Nr_of_crossings * 4):
        Crossings_set.add(i)

    Edges_set = set()
    
    for i in range(Nr_of_ropes):
        for vertex in Rope_paths[i][0]:
            a = (vertex[0] - 1) * 4 + vertex[1] - 1
            vertex2 = Connections_dict[tuple(vertex)]
            b = (vertex2[0] - 1) * 4 + vertex2[1] - 1
            if vertex[0] == vertex2[0]:
                continue
            else:
                Edges_set.add((a, b))

    Graph_cycle_lengths = []
    Graph_cycles_max_length = 0

    for cycle in Cycles_list:
        Graph_cycle_lengths.append([len(cycle)])
        if len(cycle) < Graph_cycles_max_length:
            continue
        else:
            Graph_cycles_max_length = len(cycle)
    
    Longest_cycles = set()

    for i in range(Nr_of_faces):
        if Graph_cycle_lengths[i][0] < Graph_cycles_max_length:
            continue
        else:
            Longest_cycles.add(int(i))
    
    Connectionpoints_order = []

    for r in range(Nr_of_ropes):
        This_list = []
        for i in range(len(Rope_paths[r][0])):
            This_list.append((Rope_paths[r][0][i][0] - 1) * 4 + Rope_paths[r][0][i][1] - 1)
        Connectionpoints_order.append(This_list)

    for cycle in Longest_cycles:
        Vertices_in_cycle = set()
        Vertices_in_cycle_array = []
        Edges_in_cycle = set()
        Temporary_graph = Graph(directed=True)
        temp_v_list = Temporary_graph.add_vertex(Nr_of_crossings * 4)
        temp_vertex_name = Temporary_graph.new_vertex_property("string")

        for vertex in Cycles_list[cycle]:
            Vertices_in_cycle.add((vertex[0] - 1) * 4 + vertex[1] - 1)
            Vertices_in_cycle_array.append((vertex[0] - 1) * 4 + vertex[1] - 1)
        
        for i in range(Nr_of_crossings):
            if i * 4 in Vertices_in_cycle and i * 4 + 1 in Vertices_in_cycle:
                e = Temporary_graph.add_edge(i * 4, i * 4 + 1)
            elif i * 4 + 1 in Vertices_in_cycle and i * 4 + 2 in Vertices_in_cycle:
                e = Temporary_graph.add_edge(i * 4 + 1, i * 4 + 2)
            elif i * 4 + 2 in Vertices_in_cycle and i * 4 + 3 in Vertices_in_cycle:
                e = Temporary_graph.add_edge(i * 4 + 2, i * 4 + 3)
            elif i * 4 + 3 in Vertices_in_cycle and i * 4 in Vertices_in_cycle:
                e = Temporary_graph.add_edge(i * 4 + 3, i * 4)
            else:
                continue

        for edge in Edges_set:
            if edge[0] in Vertices_in_cycle and edge [1] in Vertices_in_cycle:
                Edges_in_cycle.add(edge)
        
        for i in range(Nr_of_crossings):
            for j in range(4):
                string1 = str(i+1)
                if j == 0:
                    string2 = u'\u2081'         #unicode text is used to print subscript, for the connectionpoints
                elif j == 1:
                    string2 = u'\u2082'         #"u208a" will give "a" as subscript, where "a" goes from 1-9
                elif j == 2:
                    string2 = u'\u2083'         #for example "u2083" "gives underscore 3" can be printed in the output text too
                elif j == 3:
                    string2 = u'\u2084'
                else:
                    print("Connectionpoint index is other than 1,2,3,4")        #just a check
                string3 = string2.encode('utf-8')
                strings = "".join((string1, string3))
                temp_vertex_name[Temporary_graph.vertex(4*i + j)] = strings      #this line will add the property to that vertex, the name
        
        for edge in Edges_in_cycle:
            e = Temporary_graph.add_edge(edge[0], edge[1])
        
        pos = graph_tool.draw.sfdp_layout(Temporary_graph)
        pin = Temporary_graph.new_vertex_property("boolean")

        for i in range(Nr_of_crossings * 4):
            if i in Vertices_in_cycle:
                pin[Temporary_graph.vertex(i)] = 1
            else:
                pin[Temporary_graph.vertex(i)] = 0

        for i in range(Nr_of_crossings * 4):
            if i in Vertices_in_cycle:
                continue
            else:
                pos[i][0] = 0
                pos[i][1] = 0
        shift = float(72.0 / float(Graph_cycles_max_length))
        if Vertices_in_cycle_array[0] % 4 == Vertices_in_cycle_array[1] % 4:
            sign = 1.0
        else:
            sign = -1.0

        for i in range(len(Vertices_in_cycle)):
            r = 1000.0
            theta = float((360.0 / float(Graph_cycles_max_length)) * float(i + 1) + 44.0 + (shift * sign))
            pos[Vertices_in_cycle_array[i]][0] = r * math.degrees(math.sin(math.radians(theta)))
            pos[Vertices_in_cycle_array[i]][1] = r * math.degrees(math.cos(math.radians(theta)))
            sign = sign * (-1.0)

        for edge in Edges_set:
            e = Temporary_graph.add_edge(edge[0], edge[1])
        for i in range(Nr_of_crossings):
            e = Temporary_graph.add_edge(i * 4, i * 4 + 1)
            e = Temporary_graph.add_edge(i * 4 + 1, i * 4 + 2)
            e = Temporary_graph.add_edge(i * 4 + 2, i * 4 + 3)
            e = Temporary_graph.add_edge(i * 4 + 3, i * 4)
        
        Crossing_groups = Temporary_graph.new_vertex_property("int32_t")
        for i in range(Nr_of_crossings):
            Crossing_groups[i * 4] = i
            Crossing_groups[i * 4 + 1] = i
            Crossing_groups[i * 4 + 2] = i
            Crossing_groups[i * 4 + 3] = i

        pos2 = graph_tool.draw.sfdp_layout(Temporary_graph, pos = pos, pin = pin, K = 100, groups = Crossing_groups)
        picture = filename + "__" + str(cycle) + ".png"

        Connectionpoint_coordinates = []

        for r in range(Nr_of_ropes):
            This_list = []
            for i in range(len(Rope_paths[r][0])):
                This_list.append([i, pos[Connectionpoints_order[r][i]][0], pos[Connectionpoints_order[r][i]][1]])
            Connectionpoint_coordinates.append(This_list)

        print Connectionpoint_coordinates

        graph_draw(Graph_of_the_knot, pos = pos2, vertex_text=temp_vertex_name, vertex_font_size=60, edge_pen_width = 25, output_size=(3000, 3000), output=picture)
#        graph_draw(GConnectionpoints, pos = pos2, vertex_text=temp_vertex_name, vertex_font_size=60, edge_pen_width = 25, output_size=(3000, 3000), output=picture)

        #print Cycles_list[cycle]
        #print Vertices_in_cycle
#innen
    it_is_planar, Graphs_embed_order = graph_tool.topology.is_planar(Graph_of_the_knot, embedding=True)
#    for i in range(Nr_of_crossings * 4):
#        print Graphs_embed_order[i][0], Graphs_embed_order[i][1], Graphs_embed_order[i][2]

#    planar_position = graph_tool.draw.planar_layout(Graph_of_the_knot, pos=None)
#    graph_draw(Graph_of_the_knot, vertex_text=vertex_name, vertex_font_size=60, output_size=(2000, 2000), output=filename_pos_cyc_png)
#            graph_draw(Graph_of_the_knot, pos = pos, pin = pin, vertex_text=vertex_name, vertex_font_size=60, output_size=(2000, 2000), output=filename_pos_cyc_png)
#            p, embed_order = is_planar(Graph_of_the_knot, embedding=True)
#            print(list(embed_order[Graph_of_the_knot.vertex(0)]))
#            graph_draw(Graph_of_the_knot, vertex_text=vertex_name, vertex_font_size=60, output_size=(2000, 2000), output=filename_pos_cyc_png)

#    graph_draw(Graph_of_the_knot, vertex_text=vertex_name, vertex_font_size=60, output_size=(2000, 2000), output=filename_png)
#    graph_tool.draw.interactive_window(Graph_of_the_knot, vertex_text=vertex_name, vertex_font_size=30)


#    graph_draw(GCrossings, pos = pos, nodesfirst = True, vertex_text=vprop2_string, vertex_font_size=60, output_size=(2000, 1000), output=filename_png)
#    graph_draw(GCrossingSkeleton, pos = pos, nodesfirst = True, vertex_text=vprop2_string, vertex_font_size=60, output_size=(2000, 1000), output=filename_2png)

#    graph_draw(GCrossingSkeleton, nodesfirst = True, vertex_text=vprop3_string, vertex_font_size=60, output_size=(2000, 1000), output=filename_png)

#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################
#########################################################################################################################################################


#    TheGraph = Graph(directed=True)                     #A directed graph is defined
#    vlist = TheGraph.add_vertex(Nr_of_crossings)        #Adding the number of vertices for the crossings, = the number of crossings

#    vprop_string = TheGraph.new_vertex_property("string")   #Creating a property for the vertices that is a string
#    for i in range(Nr_of_crossings):                        #this will be used to label the vertices
#        vprop_string[TheGraph.vertex(i)] = i + 1            #It's easy for the crossings, just take the vertex number and add 1 (starts with 0)

#    for r in range(Nr_of_ropes):
#        startP = tuple(Starting_points[r])
#        moveP = Connections_dict[startP]                    #Following the connections dictionary
#        a = startP[0]                                       #using the starting points in each rope
#        c = moveP[0]                                        #a moving point is defined that will follow the path
#        if a != c:
#            print("Error: Starting point is connected to a different crossings connectionpoint!") #just a check, if fails, there is big input problem
#        for i in range(Rope_path_length[r]-1):
#            a = moveP[0]
#            moveP = Connections_dict[moveP]
#            c = moveP[0]
#            if a != c:                                      #if the crossing number of the two connectionpoints is not the same
#                e = TheGraph.add_edge(a - 1, c - 1)         #the edge is added to connect those vertices (otherwise it would self connect)

#    vlist = TheGraph.add_vertex((4*Nr_of_crossings))    #Adding the number of vertices for connectionpoints, it's 4* the number of crossings

#    for i in range(Nr_of_crossings):
#        for j in range(4):
#            string1 = str(i+1)
#            if j == 0:
#                string2 = u'\u2081'         #unicode text is used to print subscript, for the connectionpoints
#            elif j == 1:
#                string2 = u'\u2082'         #"u208a" will give "a" as subscript, where "a" goes from 1-9
#            elif j == 2:
#                string2 = u'\u2083'         #for example "u2083" "gives underscore 3" can be printed in the output text too
#            elif j == 3:
#                string2 = u'\u2084'
#            else:
#                print("Connectionpoint index is other than 1,2,3,4")        #just a check
#            string3 = string2.encode('utf-8')
#            strings = "".join((string1, string3))
#            vprop_string[TheGraph.vertex(Nr_of_crossings + 4*i + j)] = strings      #this line will add the property to that vertex, the name

#    vlist = TheGraph.add_vertex(Nr_of_faces)    #Adding vertices for the dual graph

#    for i in range(Nr_of_faces):
#        vprop_string[TheGraph.vertex(Nr_of_crossings*5 + i)] = i + 1        #Labeling the vertices of the dual that is the faces in this case

#    TheGraph.remove_vertex(Nr_of_crossings * 5 + Outside_face - 1)        #Removing the outside face vertex of the dual

#    pos = graph_tool.draw.sfdp_layout(GCrossingsDual)       #This will calculate the position of the vertices of the dual graph, might need to scale it

#    pin = TheGraph.new_vertex_property("boolean")
#    for i in range(Nr_of_crossings * 5 + Nr_of_faces - 1):
#        if i < (Nr_of_crossings * 5):
#            pin[TheGraph.vertex(i)] = 0
#        else:
#            pin[TheGraph.vertex(i)] = 1

    #Adding edges

#    for i in range(Nr_of_faces):
#        if i == Outside_face - 1:
#            continue
#        elif i != Outside_face - 1:
#            cycle_length = len(Cycles_list[i])
#            for j in range(cycle_length):
#                a = Cycles_list[i][j][0]
#                b = Cycles_list[i][j][1]
#                TheGraph.add_edge(Nr_of_crossings * 5 + i, Nr_of_crossings + 4 * (a - 1) + b - 1)

#    for i in range(Nr_of_crossings):
#        e = TheGraph.add_edge(Nr_of_crossings - 1 + 4*i + 1, Nr_of_crossings - 1 + 4*i + 2)
#        e = TheGraph.add_edge(Nr_of_crossings - 1 + 4*i + 2, Nr_of_crossings - 1 + 4*i + 3)
#        e = TheGraph.add_edge(Nr_of_crossings - 1 + 4*i + 3, Nr_of_crossings - 1 + 4*i + 4)
#        e = TheGraph.add_edge(Nr_of_crossings - 1 + 4*i + 4, Nr_of_crossings - 1 + 4*i + 1)
    
    #Creating new vertex property, this will be used to group them so the printing software can put the vertices that are in the same group close to eachother
#    group_prop = TheGraph.new_vertex_property("int32_t")

#    for i in range(5 * Nr_of_crossings):
#        if i < Nr_of_crossings:
#            continue
#        else:
#            group_prop[TheGraph.vertex(i)] = int(vprop_string[i][0])    #Add property to the connectionpoints

#    position = graph_tool.draw.sfdp_layout(TheGraph, groups=group_prop)

#    for i in range(Nr_of_faces):
#        if i == Outside_face - 1:
#            continue
#        else:
#            position[5 * Nr_of_crossings + i][0] = pos[i][0]
#            position[5 * Nr_of_crossings + i][1] = pos[i][1]

#    position = graph_tool.draw.sfdp_layout(TheGraph, pos = position, pin = pin, groups=group_prop)

#    graph_draw(TheGraph, pos = position, nodesfirst = True, vertex_text=vprop_string, vertex_font_size=60, output_size=(1000, 1000), output=filename_2png)

    return Crossings_adjacency_matrix
##########################################################################################################
##########################################################################################################



##########################################################################################################
################################################################################ NextCrossing
def NextCrossing(Connections_dict, Crossings_chirality_dict):
    '''Defining the next crossings chirality; '''
    Next_crossings_point = Connections_dict[Connections_dict[(1, 1)]]
    Next_crossings_chirality = Crossings_chirality_dict[Next_crossings_point[0] - 1]
    #Next_crossings_alternation is 1 if the rope goes over, and -1 if it goes under (the meaning of it is that which part of the crossing is first in the rope path, the over or under part)
    if Next_crossings_point[1] == 1:
        Next_crossings_alternation = 1
    if Next_crossings_point[1] == 2:
        Next_crossings_alternation = -1
    if Next_crossings_point[1] == 4:
        Next_crossings_alternation = -1
    if Next_crossings_point[1] == 3:
        Next_crossings_alternation = 0
        print 'Invalid parameter "3" given as the next crossings point' #It is impossible that we go to 3 first, because 3 allways points towards the next crossing)
    return Next_crossings_point, Next_crossings_chirality, Next_crossings_alternation
##########################################################################################################
##########################################################################################################

##########################################################################################################
################################################################################ NextRepresentationsChirality
def NextRepresentationsChirality(Nr_of_crossings, Crossings_chirality_dict, Next_crossings_alternation):
    '''Create the next representations chirality array; '''
    New_chirality = numpy.zeros((Nr_of_crossings), numpy.int32)
    for a in range(Nr_of_crossings):
        New_chirality[a] = Crossings_chirality_dict[a] * Next_crossings_alternation #If the next crossings alternation is 1, we won't change the chirality, but if it's -1, the chirality will be changed for every crossing
    return New_chirality
##########################################################################################################
##########################################################################################################

##########################################################################################################
################################################################################ CrossingsChangeDict
def CrossingsChangeDict(Nr_of_crossings, Connections_dict, Starting_points):
    '''Creating dictionary of the crossings number change, and the inverse of that; '''
    Crossings_change_dict = {}
    for i in range(Nr_of_crossings):
        Crossings_change_dict[i + 1] = 0

    StartingP = Connections_dict[Connections_dict[tuple(Starting_points[0])]]
#    EndingP = tuple(Starting_points[0]) #No need for this, I'm not using it
    MovingP = Connections_dict[Connections_dict[tuple(Starting_points[0])]]

    Crossings_counter = 0
    Rope_counter = 0

    for i in range(Nr_of_crossings*2): #Nr_of_crossings *2, because we need to follow the path under (*1) and (+) over (*1) to get back to the starting point
        if Crossings_change_dict[MovingP[0]] == 0:
            Crossings_change_dict[MovingP[0]] = Crossings_counter + 1
            Crossings_counter += 1
            MovingP = Connections_dict[Connections_dict[MovingP]]
        else:
            MovingP = Connections_dict[Connections_dict[MovingP]]
        
        if MovingP == StartingP:
            if Rope_counter < Nr_of_ropes - 1:
                Rope_counter += 1
                MovingP = tuple(Starting_points[Rope_counter])
                StartingP = tuple(Starting_points[Rope_counter])

    Crossings_change_dict_inv = {Crossings_change_dict[key] : key for key in Crossings_change_dict}

    ConnectionP_change_clock = {2: 1, 1: 4, 3: 2, 4: 3} #For clockwise, the chirality is taken as -1
    ConnectionP_change_anti = {2: 3, 1: 2, 3: 4, 4: 1} #For anticlockwise, the chirality is taken as 1

    #How the numbers will change if the chirality is changed. 
    #(Double lines represent the rope, single lines show the transformation)
    #Note that this change is symmetric, so clocwise dict is the inverse of anticlockwise dict.
    #
    #       --------------------------
    #       |                        |
    #       |     2 ---------------------- 1
    #       |     ||                 |     ||
    #             ||                       ||
    #       1 ========> 3            4 ====||===> 2
    #             ||                       ||   ^
    #             ||  / |                  ||  /  |
    #             V  v  |                  V      |
    #             4 ---------------------- 3      |
    #                   |                         |
    #                   ---------------------------
    #
    #       Clockwise = -1   <-->   Anticlockwise = 1

    return Crossings_change_dict, Crossings_change_dict_inv, ConnectionP_change_clock, ConnectionP_change_anti
##########################################################################################################
##########################################################################################################

##########################################################################################################
################################################################################ 
def NextRepresentationsMatrix(Nr_of_crossings, Nr_of_ropes, Connections_matrix, Starting_points, New_chirality, Next_crossings_alternation, Crossings_change_dict, Crossings_change_dict_inv, ConnectionP_change_clock, ConnectionP_change_anti, Crossings_chirality_dict):
    '''Creating the next representations matrix; '''
    ConMatStart = Connections_matrix.copy()
    ConMatEnd = ConMatStart.copy()
    New_chiralities = New_chirality.copy()

    for i in range(Nr_of_crossings):
        New_chiralities[i] = New_chirality[Crossings_change_dict_inv[i + 1] - 1]

    New_chiralities_dict = {}
    for i in range(Nr_of_crossings):
        New_chiralities_dict[i] = New_chiralities[i]

    for i in range(Nr_of_crossings):
        for j in range(4):
            ConMatEnd[i,j,0] = Crossings_change_dict[ConMatEnd[i,j,0]]

    ConMatMid = ConMatEnd.copy()
    for i in range(Nr_of_crossings):
        for j in range(4):
            if Next_crossings_alternation == -1:
                if New_chiralities_dict[ConMatEnd[i,j,0] - 1] == 1:
                    ConMatEnd[i,j,1] = ConnectionP_change_clock[ConMatEnd[i,j,1]]
                if New_chiralities_dict[ConMatEnd[i,j,0] - 1] == -1:
                    ConMatEnd[i,j,1] = ConnectionP_change_anti[ConMatEnd[i,j,1]]

    ConMatMid = ConMatEnd.copy()
    for i in range(Nr_of_crossings):
        for j in range(4):
            if Next_crossings_alternation == -1:
                if Crossings_chirality_dict[i] == 1:
                    ConMatEnd[i,j,0] = ConMatMid[i, ConnectionP_change_clock[j + 1] - 1, 0]
                    ConMatEnd[i,j,1] = ConMatMid[i, ConnectionP_change_clock[j + 1] - 1, 1]
                if Crossings_chirality_dict[i] == -1:
                    ConMatEnd[i,j,0] = ConMatMid[i, ConnectionP_change_anti[j + 1] - 1, 0]
                    ConMatEnd[i,j,1] = ConMatMid[i, ConnectionP_change_anti[j + 1] - 1, 1]

    ConMatMid = ConMatEnd.copy()
    for i in range(Nr_of_crossings):
        for j in range(4):
            ConMatEnd[i,j,0] = ConMatMid[Crossings_change_dict_inv[i + 1] - 1, j, 0]
            ConMatEnd[i,j,1] = ConMatMid[Crossings_change_dict_inv[i + 1] - 1, j, 1]

    Next_representation_matrix = ConMatEnd.copy()
    
# This part is not working for 
#    Next_starting_points = Starting_points.copy()
#    for i in range(Nr_of_ropes):
#        print Next_starting_points[i, 0]
#        print Next_starting_points[i, 1]
#        something = Connections_dict[Connections_dict[(Next_starting_points[1, 0], Next_starting_points[i, 1])]]
#        print something
#        Next_starting_points[i, 0] = something[0]
#        Next_starting_points[i, 1] = something[1]
#        print Next_starting_points
        
    return Next_representation_matrix#, Next_starting_points
##########################################################################################################
##########################################################################################################



#End of knot representation generating functions
##################################################################################################################
##################################################################################################################
#Start of code

std_err = sys.stderr
std_err.write("\nSTART\nOpening the standard input file; ")
fn_in = sys.stdin
std_err.write("Opening the standard output file; ")
out_put = fn_out

#Get data from inputfile, turn it into arrays, matrices
in_put, Nr_of_ropes, Starting_points, Nr_of_crossings, Connections_matrix, Position_matrix = DataFromInput(std_err, fn_in, out_put)

def GetNumberOfEdges(Nr_of_crossings):
    '''Calculating Nr_of_edges and Nr_of_faces; '''
    Nr_of_degrees = Nr_of_crossings * 4 #Every crossing has 4 degrees, two in and two out, since the graph is 4-regular
    Nr_of_edges = Nr_of_degrees / 2     #Simple way of calculating the number of edges from the number of degrees, every edge has two ends
    Nr_of_faces = Nr_of_edges - Nr_of_crossings + 2 #Since knot projections are planar, the Euler characteristic is 2, ==> V-E+F=2 (F=E-V+2)
    return Nr_of_edges, Nr_of_faces

Nr_of_edges, Nr_of_faces = GetNumberOfEdges(Nr_of_crossings)

out_put.write("Based on the Euler formula (V - E + F = 2) the Crossings 4-regular planar graph has: (R - ropes)\n")
out_put.write("\nV = {}".format(Nr_of_crossings))
out_put.write("\nE = {}".format(Nr_of_edges))
out_put.write("\nF = {}".format(Nr_of_faces))
out_put.write("\nR = {}\n\n".format(Nr_of_ropes))

###################################################################################################################################

#Take the arrays and matrices, to create connectionpoints connectivity matrix
Rope_paths, Rope_path_length, Connections_dict, Crossings_chirality, Crossings_chirality_dict, Connections_dict_inv, Connectionpoints_connectivity_matrix = ConnectionpointsConnectivityMatrix(std_err, Nr_of_ropes, Nr_of_crossings, Starting_points, Connections_matrix)

#Adding the symbols A (for "Above") and B (for "Below"), and calculating the determinant/polynomial
Connectionpoints_connectivity_matrix_AB, Determinant = GetPolynomial(Nr_of_crossings, Crossings_chirality, Connections_dict_inv, Connections_dict, A, B)

An = []
Bn = []

for i in range(Nr_of_crossings):
    An.append ("A" + "{}".format(i + 1))
    An[i] = sympy.Symbol("A" + "{}".format(i + 1))
    Bn.append ("B" + "{}".format(i + 1))
    Bn[i] = sympy.Symbol("B" + "{}".format(i + 1))

CCMn = sympy.SparseMatrix(Nr_of_crossings*4, Nr_of_crossings*4, {(3, 3): A})
CCMn = CCMn.add(-CCMn)

for i in range(Nr_of_crossings):
    c = Crossings_chirality[i]
    for j in range(4):
        m = i*4 + j
        if c == 1:
            if j == 0:
                a = Connections_dict_inv[tuple([i+1,j+1])]
                n = (a[0]-1)*4 + a[1] - 1
                CCMn = CCMn.add(sympy.SparseMatrix(Nr_of_crossings*4, Nr_of_crossings*4, {(m,n): -An[i]}))
            elif j == 1:
                a = Connections_dict[tuple([i+1,j+1])]
                n = (a[0]-1)*4 + a[1] - 1
                CCMn = CCMn.add(sympy.SparseMatrix(Nr_of_crossings*4, Nr_of_crossings*4, {(m,n): Bn[i]}))
            elif j == 2:
                a = Connections_dict[tuple([i+1,j+1])]
                n = (a[0]-1)*4 + a[1] - 1
                CCMn = CCMn.add(sympy.SparseMatrix(Nr_of_crossings*4, Nr_of_crossings*4, {(m,n): An[i]}))
            elif j == 3:
                a = Connections_dict_inv[tuple([i+1,j+1])]
                n = (a[0]-1)*4 + a[1] - 1
                CCMn = CCMn.add(sympy.SparseMatrix(Nr_of_crossings*4, Nr_of_crossings*4, {(m,n): -Bn[i]}))
        elif c == -1:
            if j == 0:
                a = Connections_dict_inv[tuple([i+1,j+1])]
                n = (a[0]-1)*4 + a[1] - 1
                CCMn = CCMn.add(sympy.SparseMatrix(Nr_of_crossings*4, Nr_of_crossings*4, {(m,n): -An[i]}))
            elif j == 1:
                a = Connections_dict_inv[tuple([i+1,j+1])]
                n = (a[0]-1)*4 + a[1] - 1
                CCMn = CCMn.add(sympy.SparseMatrix(Nr_of_crossings*4, Nr_of_crossings*4, {(m,n): -Bn[i]}))
            elif j == 2:
                a = Connections_dict[tuple([i+1,j+1])]
                n = (a[0]-1)*4 + a[1] - 1
                CCMn = CCMn.add(sympy.SparseMatrix(Nr_of_crossings*4, Nr_of_crossings*4, {(m,n): An[i]}))
            elif j == 3:
                a = Connections_dict[tuple([i+1,j+1])]
                n = (a[0]-1)*4 + a[1] - 1
                CCMn = CCMn.add(sympy.SparseMatrix(Nr_of_crossings*4, Nr_of_crossings*4, {(m,n): Bn[i]}))
Connectionpoints_connectivity_matrix_AnBn = sympy.SparseMatrix(CCMn)

Connectionpoints_connectivity_matrix_small_AnBn = sympy.zeros(Nr_of_crossings, Nr_of_crossings)
Ccmsn = sympy.SparseMatrix(Connectionpoints_connectivity_matrix_small_AnBn)

for i in range(Nr_of_crossings):
    for j in range(Nr_of_crossings):
        m = Connectionpoints_connectivity_matrix_AnBn
        a = m[4*i,(4*j)] + m[4*i,(4*j+1)] + m[4*i,(4*j+2)] + m[4*i,(4*j+3)]
        b = m[(4*i+1),(4*j)] + m[(4*i+1),(4*j+1)] + m[(4*i+1),(4*j+2)] + m[(4*i+1),(4*j+3)]
        c = m[(4*i+2),(4*j)] + m[(4*i+2),(4*j+1)] + m[(4*i+2),(4*j+2)] + m[(4*i+2),(4*j+3)]
        d = m[(4*i+3),(4*j)] + m[(4*i+3),(4*j+1)] + m[(4*i+3),(4*j+2)] + m[(4*i+3),(4*j+3)]
        Ccmsn = Ccmsn.add(sympy.SparseMatrix(Nr_of_crossings, Nr_of_crossings, {(i,j): a + b + c + d}))
for i in range(Nr_of_crossings):
    Ccmsn = Ccmsn.add(sympy.SparseMatrix(Nr_of_crossings, Nr_of_crossings, {(i,i): int(Crossings_chirality[i])}))

Connectionpoints_connectivity_matrix_small_AnBn = sympy.SparseMatrix(Ccmsn)
out_put.write("\n\nAn - Bn Symbols matrix:\n")
out_put.write("\n{}\n".format(Connectionpoints_connectivity_matrix_small_AnBn))
#AnBnDeterminant = sympy.SparseMatrix.det(Connectionpoints_connectivity_matrix_small_AnBn)
#out_put.write("\nAn - Bn Symbols determinant:\n")
#out_put.write("\n{}\n".format(AnBnDeterminant))
#kikapcsolva

###################################################################################################################################

#Create dictionary for the important data, for every representation
CodRopCroStaConPatDicChiConSymPol = (Nr_of_ropes, Nr_of_crossings, Starting_points, Connections_matrix, Rope_paths, Connections_dict, Crossings_chirality, Crossings_chirality_dict, Connectionpoints_connectivity_matrix, Connectionpoints_connectivity_matrix_AB, Determinant)
#out_put.write("\nCodRopCroStaConPatDicChiConSymPol:\n\n{}\n".format(CodRopCroStaConPatDicChiConSymPol))
All_representations = {}

All_representations[0] = CodRopCroStaConPatDicChiConSymPol

###################################################################################################################################

#Create the next representation
#Follow the path and define what is the next crossings point going to be, based on the path. 
#Path has to be changed before this step if we want the rope to go in the opposite direction!!!
Next_crossings_point, Next_crossings_chirality, Next_crossings_alternation = NextCrossing(Connections_dict, Crossings_chirality_dict)
#Crating the new chirality
New_chirality = NextRepresentationsChirality(Nr_of_crossings, Crossings_chirality_dict, Next_crossings_alternation)

#Creating the dictionary of the crossings number change and the inverse of that
Crossings_change_dict, Crossings_change_dict_inv, ConnectionP_change_clock, ConnectionP_change_anti = CrossingsChangeDict(Nr_of_crossings, Connections_dict, Starting_points)

#Creating the next representations matrix
Next_representation_matrix = NextRepresentationsMatrix(Nr_of_crossings, Nr_of_ropes, Connections_matrix, Starting_points, New_chirality, Next_crossings_alternation, Crossings_change_dict, Crossings_change_dict_inv, ConnectionP_change_clock, ConnectionP_change_anti, Crossings_chirality_dict)


out_put.write("\nNext representations connectivity matrix:\n\n")
for i in range(Nr_of_crossings):
        out_put.write("[[{} {}] [".format(Next_representation_matrix[i,0,0], Next_representation_matrix[i,0,1]))
        out_put.write("{} {}] [".format(Next_representation_matrix[i,1,0], Next_representation_matrix[i,1,1]))
        out_put.write("{} {}] [".format(Next_representation_matrix[i,2,0], Next_representation_matrix[i,2,1]))
        out_put.write("{} {}]]\n".format(Next_representation_matrix[i,3,0], Next_representation_matrix[i,3,1]))
###################################################################################################################################

#Create a simple crossings connectivity (adjacency) matrix

Crossings_adjacency_matrix = CrossingsAdjacencyMatrix(Nr_of_crossings, Nr_of_ropes, Starting_points, Rope_path_length, Crossings_chirality, Connections_dict_inv, Connections_dict)

#Take the arrays and matrices, to create connectionpoints connectivity matrix
#Next_Rope_paths, Next_Connections_dict, Next_Crossings_chirality, Next_Crossings_chirality_dict, Next_Connections_dict_inv, Next_Connectionpoints_connectivity_matrix = ConnectionpointsConnectivityMatrix(std_err, Nr_of_ropes, Nr_of_crossings, Next_starting_points, Next_representation_matrix)

#Adding the symbols A (for "Above") and B (for "Below"), and calculating the determinant/polynomial
#Connectionpoints_connectivity_matrix_AB, Determinant = GetPolynomial(Nr_of_crossings, Crossings_chirality, Connections_dict_inv, Connections_dict, A, B)

###################################################################################################################################

#Create dictionary for the important data, for every representation
#CodRopCroStaConPatDicChiConSymPol = (Nr_of_ropes, Nr_of_crossings, Starting_points, Connections_matrix, Rope_paths, Connections_dict, Crossings_chirality_dict, Connectionpoints_connectivity_matrix, Connectionpoints_connectivity_matrix_AB, Determinant)

#All_representations = {}

#All_representations[0] = CodRopCroStaConPatDicChiConSymPol

std_err.write( "\nEND\n")

fn_in.close()
#fn_all_rep.close()
fn_out.close()

exit(0)