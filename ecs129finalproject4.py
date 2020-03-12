acidSequence = {'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,
                'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19}

import math
def getSeqFile(file_in):

    fp_in = open(file_in, "r")
    #skips first line of fasta file
    fp_in.readline()
    #stores rest of fasta file
    seq = fp_in.readlines()

    #stores sequence from file into an array
    i=0
    seqTemp = ""
    while i < len(seq):
        seqTemp = seqTemp+(seq[i])
        i = i+1

    seq = ''
    for char in seqTemp:
        if char != '\n':
            seq = seq+char


    return seq

def getMatrix(beta):

    file_in = input("What is the name of the input file to get the matrix\n")
    fp_in = open(file_in, "r")
    temp = fp_in.readline()

    #take out spaces between amino acids
    temp = temp.replace(" ","")
    temp = temp.replace("\n","")

    #create matrix
    size = len(temp)+1 #number of amino acids + 1 (for the labels)
    tempMatrix = [(['*'] * size) for row in range(size)]

    #fill top row with Amino Acids
    i = 1
    j = 0
    while j < len(temp):
        tempMatrix[0][i] = temp[j]
        j = j+1
        i = i+1



    # fill rest of matrix

    k= 1 #the row being filled
    while k < size:
        temp = fp_in.readline()
        temp = temp.split()


        i = 0
        j = 0
        while j < len(temp):
            if temp[j].isalpha()==True:
                tempMatrix[k][i] = temp[j]
                j = j + 1
                i = i + 1
            else:
                base = temp[j]
                base =float(base)
                tempMatrix[k][i] = pow(base,beta)
                j = j + 1
                i = i + 1

        k = k+1


    return tempMatrix


# conver to distance
def donkeyKong(matrix, sequenceA, sequenceB):
    seqA = []
    seqB = []

    for acidA in sequenceA:
        seqA.append(acidSequence[acidA])

    for acidB in sequenceB:
        seqB.append(acidSequence[acidB])
    dk = math.sqrt(2 * (1 - k3Tilda(matrix, seqA, seqB)))
    # print('dk = ', dk)
    return dk


# normalize for length
def k3Tilda(matrix, seqA, seqB):
    final = kay3(matrix, seqA, seqB) / math.sqrt((kay3(matrix, seqA, seqA) * kay3(matrix, seqB, seqB)))
    # print("final", final)

    return final


def kay1(matrix, acidA, acidB):
    return (matrix[acidA + 1][acidB + 1])


def kay2(matrix, sliceSeqA, sliceSeqB):
    k2 = 1

    for i in range(len(sliceSeqA)):
        temp = kay1(matrix, sliceSeqA[i], sliceSeqB[i])

        k2 = k2 * temp

    return k2


def kay3(matrix, seqA, seqB):
    if len(seqA) <= len(seqB):
        shorter = seqA
        longer = seqB
    else:
        shorter = seqB
        longer = seqA

    k3 = 0
    counter = 1
    i = 1
    s = 0
    while counter <= min(10, len(shorter)):

        s = 0
        i = counter

        while i <= len(shorter):

            window = shorter[s:i]
            j = 0

            while j + len(window) - 1 < len(longer):
                k3 = k3 + kay2(matrix, window, longer[j:j + len(window)])
                j += 1
            i += 1
            s = i - counter

        counter += 1

    return k3


#make NxN matrix of
def groupCompare(matrix, orgList):

    outer = []
    animalAindex = 0


    for animalA in orgList:
        inner = []

        for animalB in orgList:
            dkResult = donkeyKong(matrix, animalA,animalB)
            inner.append(dkResult)
        outer.append(inner)
        animalAindex+=1

    return outer


# this segent of code was based off of a geeks for gees page link
from collections import defaultdict


# Class to represent a graph
class Graph:

    def __init__(self, vertices):
        self.V = vertices  # No. of vertices
        self.graph = []  # default dictionary
        # to store graph

    # function to add an edge to graph
    def addEdge(self, u, v, w):
        self.graph.append([u, v, w])

        # A utility function to find set of an element i

    # (uses path compression technique)
    def find(self, parent, i):
        if parent[i] == i:
            return i
        return self.find(parent, parent[i])

        # A function that does union of two sets of x and y

    # (uses union by rank)
    def union(self, parent, rank, x, y):
        xroot = self.find(parent, x)
        yroot = self.find(parent, y)

        # Attach smaller rank tree under root of
        # high rank tree (Union by Rank)
        if rank[xroot] < rank[yroot]:
            parent[xroot] = yroot
        elif rank[xroot] > rank[yroot]:
            parent[yroot] = xroot

            # If ranks are same, then make one as root
        # and increment its rank by one
        else:
            parent[yroot] = xroot
            rank[xroot] += 1

    # The main function to construct MST using Kruskal's
    # algorithm
    def KruskalMST(self):

        result = []  # This will store the resultant MST

        i = 0  # An index variable, used for sorted edges
        e = 0  # An index variable, used for result[]

        # Step 1:  Sort all the edges in non-decreasing
        # order of their
        # weight.  If we are not allowed to change the
        # given graph, we can create a copy of graph
        self.graph = sorted(self.graph, key=lambda item: item[2])

        parent = [];
        rank = []

        # Create V subsets with single elements
        for node in range(self.V):
            parent.append(node)
            rank.append(0)

            # Number of edges to be taken is equal to V-1
        while e < self.V - 1:

            # Step 2: Pick the smallest edge and increment
            # the index for next iteration
            u, v, w = self.graph[i]
            i = i + 1
            x = self.find(parent, u)
            y = self.find(parent, v)

            # If including this edge does't cause cycle,
            # include it in result and increment the index
            # of result for next edge
            if x != y:
                e = e + 1
                result.append([u, v, w])
                self.union(parent, rank, x, y)
                # Else discard the edge

        NIE = []
        for result in result:
            # print str(u) + " -- " + str(v) + " == " + str(weight)
            NIE.append(result)

        return NIE


def main():

    blosum62 = [[]]

    seqNum = int(input("Please enter the number of sequences to find the distance between: "))
    orgSeqs = []
    orgStings = []

    #enters filename to get sequence into list
    for i in range(seqNum):
        orgName = input("Please enter the name of the sequence ID: ")
        orgSeq = input("Please enter sequence file name: ")
        orgSeqs.append(getSeqFile(orgSeq))
        orgStings.append(orgName)

    beta = input("please enter value for beta\n")
    beta = float(beta)

    #takes in matrix file and raises all values by beta
    blosum62 = getMatrix(beta)

    print("start watching the star wars or LOTR trillogy")

    #Actual function to compare sequences/calculate distances --> puts into a matrix
    dkGraph = groupCompare(blosum62, orgSeqs)

    print(dkGraph)

    #Generate MST from distance matrix
    answer = input("Do you want to cluster results?")
    answer.capitalize()
    if answer == 'YES':
        clusterNum = int(input("how many clusters do you want to make?"))
        #Deletes K-1 max edges
        clusterNum = clusterNum - 1

        # PutDkgraph into a graph
        # Creates graph object of distance results with organims as nodes and distances as edges
        mstTemp = Graph(len(orgSeqs))


        for orgNodeA in range(len(orgSeqs)):
            for orgNodeB in range(len(orgSeqs)):
                mstTemp.addEdge(orgNodeA, orgNodeB, dkGraph[orgNodeA][orgNodeB])

        # make MST off of graph
        # Kruskal's implementation sourced from: _____, with modifications
        mst = mstTemp.KruskalMST()

        # remove K-1 max edges
        count = 0
        while count < clusterNum:
            del mst[-1]
            count += 1

        return mst

    else:
        return






if __name__ == '__main__':
    main()