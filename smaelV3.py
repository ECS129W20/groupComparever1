acidSequence = {'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,
                'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19}

import math


def getSeqFile():

    file_in = input("What is the name of the input file to get the sequence\n")
    fp_in = open(file_in, "r")
    fp_in.readline()
    seq = fp_in.readlines()

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

    #print matrix --> optional
    i = 0
    for row in tempMatrix:
        print('row = ', i, row)
        i = i +1


    return tempMatrix


def donkeyKong(matrix, sequenceA, sequenceB):
    seqA = []
    seqB = []

    #Create Dictionary for Sequence A, where
    for acidA in sequenceA:
        seqA.append(acidSequence[acidA])

    for acidB in sequenceB:
        seqB.append(acidSequence[acidB])
    print("k3Tild", k3Tilda(matrix, seqA, seqB))
    dk = math.sqrt(2 * (1 - k3Tilda(matrix, seqA, seqB)))
    print('dk = ', dk)
    return dk


def k3Tilda(matrix, seqA, seqB):
    final = kay3(matrix, seqA, seqB) / math.sqrt((kay3(matrix, seqA, seqA) * kay3(matrix, seqB, seqB)))
    print("final", final)

    return final


def kay1(blosum62, acidA, acidB):

    return (blosum62[acidA+1][acidB+1])


def kay2(matrix, sliceSeqA, sliceSeqB):

    k2 = 1
    for i in range(len(sliceSeqA)):
        temp = kay1(matrix, sliceSeqA[i], sliceSeqB[i])

        k2 = k2 * temp

    return k2


def kay3(matrix, seqA, seqB):
    # print('sequences: ', seqA)
    # print('sequences: ', seqB)
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
    #while counter <= len(shorter):
    while counter <= min(10, len(shorter)):

        s = 0
        i = counter
        while i <= len(shorter):
            window = shorter[s:i]
            # print('shorter:', shorter)
            # print("window",window)
            j = 0
            while j + len(window) - 1 < len(longer):
                # print("j",j)
                # print('longer window: ', longer[j:j+len(window)])
                k3 = k3 + kay2(matrix, window, longer[j:j + len(window)])
                j += 1
                # print('j has been incremented',j)

            i += 1
            s = i - counter
            # print('i has been incremented: ', i)
            # print("i",i)
        counter += 1

    # print("k3",k3)
    return k3





def main():

    seq1 = []
    seq2 = []
    blosum62 = [[]]

    seq1 = getSeqFile()
    print(' ')
    seq2 = getSeqFile()
    beta = input("please enter value for beta\n")
    beta = float(beta)
    blosum62 = getMatrix(beta)


    x = donkeyKong(blosum62, seq1, seq2)
    print(x)





if __name__ == '__main__':
    main()






