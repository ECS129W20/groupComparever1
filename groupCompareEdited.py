import math
acidSequence = {'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,
                'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19}

blosum62 = [[3.90294070015052, 0.612698600061261, 0.58830764005965, 0.544605274637711, 0.867987663625351, 0.756803942764855, 0.741264113108527, 1.05686960775367, 0.569364849349247, 0.632481034524687, 0.601945974581018, 0.775390239493712, 0.723150342301527, 0.464893827242731, 0.754121369072143, 1.47210398545007, 0.984401955935144, 0.416548781469996, 0.542611868915233,0.936458396120263],
            [0.612698600061261, 6.66557706993898, 0.858630477662975, 0.573200023500285, 0.308939296217125, 1.40579606250098, 0.960797602466529, 0.449983999048674, 0.917048020652714, 0.354751311390641, 0.473919278116426, 2.07680866910689, 0.622623369170503, 0.380726330360237, 0.48153490494176, 0.767165632779031, 0.677754679194329, 0.395102105602751, 0.555965424563688, 0.42007231624752],
            [0.58830764005965, 0.858630477662975, 7.09409487818815, 1.55385280556386, 0.397802620047207, 1.00058441852805, 0.911298183018242, 0.86371140576969, 1.22200066958752, 0.327934751806163, 0.310043275665557, 0.939841128716106, 0.474529654685916, 0.354288952229033, 0.499932835964896, 1.23152924484831, 0.984152634507759, 0.277782895565922, 0.486030805784542, 0.369033853043791],
            [0.544605274637711, 0.573200023500285, 1.55385280556386, 7.39792738079911, 0.301454344668722, 0.89708112869296, 1.68781074564091, 0.634301018724725, 0.678558838858087, 0.339015407074848, 0.286613045833504, 0.784090405759266, 0.346454633856285, 0.298969081268819, 0.598716825883862, 0.913504624354301, 0.694789868062839, 0.232102315145553, 0.345683565218847, 0.336500141524566],
            [0.867987663625351, 0.308939296217125, 0.397802620047207, 0.301454344668722, 19.5765856868537, 0.3657815306054, 0.285934574128461, 0.420387869540714, 0.355049504965318, 0.653458800603899, 0.642275633431893, 0.349128464920366, 0.61135434012195, 0.438990117695358, 0.379562691207627, 0.738415701073377, 0.740551692220627, 0.449983902793404, 0.434203398141003, 0.75584405464154],
            [0.756803942764855, 1.40579606250098, 1.00058441852805, 0.89708112869296, 0.3657815306054, 6.24442175356205, 1.90173784203935, 0.538649627426744, 1.16798103533111, 0.382937802207239, 0.477325586336375, 1.5543230772441, 0.864250292645428, 0.333972401843889, 0.641280588751714, 0.96555522798098, 0.79132074056634, 0.509360271558287, 0.611094097106686, 0.466777931405709],
            [0.741264113108527, 0.960797602466529, 0.911298183018242, 1.68781074564091, 0.285934574128461, 1.90173784203935, 5.46952607963445, 0.481267654658343, 0.960040718354581, 0.330522558376655, 0.372873704285776, 1.30827885329714, 0.50034228947388, 0.33074399059478, 0.679202586642317, 0.950357185031325, 0.741425610477113, 0.374300211820363, 0.496467353893267, 0.428943129877398],
            [1.05686960775367, 0.449983999048674, 0.86371140576969, 0.634301018724725, 0.420387869540714, 0.538649627426744, 0.481267654658343, 6.87630690865387, 0.492966575788069, 0.275009721763455, 0.284504011912594, 0.588871736039716, 0.395486600257494, 0.340640908478402, 0.477385507184256, 0.90359652515418, 0.579271581711225, 0.421690355206204, 0.348714366361603, 0.3369549123791],
            [0.569364849349247, 0.917048020652714, 1.22200066958752, 0.678558838858087, 0.355049504965318, 1.16798103533111, 0.960040718354581, 0.492966575788069, 13.5059996886779, 0.326288124625136, 0.380675485808673, 0.778887489609194, 0.584132623334439, 0.651990520809943, 0.472879830723747, 0.736731739892316, 0.55750325361248, 0.44408895489718, 1.79790413031311, 0.339447441760233],
            [0.632481034524687, 0.354751311390641, 0.327934751806163, 0.339015407074848, 0.653458800603899, 0.382937802207239, 0.330522558376655, 0.275009721763455, 0.326288124625136, 3.997929939961, 1.69443475437089, 0.396372934422445, 1.4777445015865, 0.945769882931625, 0.384662859733384, 0.443163582314639, 0.779816109586742, 0.408874390481926, 0.630388930627921, 2.41751209060932],
            [0.601945974581018, 0.473919278116426, 0.310043275665557, 0.286613045833504, 0.642275633431893, 0.477325586336375, 0.372873704285776, 0.284504011912594, 0.380675485808673, 1.69443475437089, 3.79662136919197, 0.428270363066123, 1.99429556770288, 1.15459749441297, 0.37112172360338, 0.428893742551906, 0.660328974513946, 0.568037074030476, 0.692059423023677, 1.31423572845207],
            [0.775390239493712, 2.07680866910689, 0.939841128716106, 0.784090405759266, 0.349128464920366, 1.5543230772441, 1.30827885329714, 0.588871736039716, 0.778887489609194, 0.396372934422445, 0.428270363066123, 4.76433717338922, 0.625302816237689, 0.344043118911871, 0.703774478956202, 0.931919140710646, 0.792905802702891, 0.358930070683472, 0.532179332619096, 0.456542719723346],
            [0.723150342301527, 0.622623369170503, 0.474529654685916, 0.346454633856285, 0.61135434012195, 0.864250292645428, 0.50034228947388, 0.395486600257494, 0.584132623334439, 1.4777445015865, 1.99429556770288, 0.625302816237689, 6.48145120779836, 1.00437163122058, 0.423898023856973, 0.598558924100088, 0.793801615982561, 0.61029621403986, 0.708364627674993, 1.26893679116311],
            [0.464893827242731, 0.380726330360237, 0.354288952229033, 0.298969081268819, 0.438990117695358, 0.333972401843889, 0.33074399059478, 0.340640908478402, 0.651990520809943, 0.945769882931625, 1.15459749441297, 0.344043118911871, 1.00437163122058, 8.1287970162524, 0.287444757613266, 0.439973596778216, 0.481693682890345, 1.37437942379832, 2.76938062915766, 0.745089737790822],
            [0.754121369072143, 0.48153490494176, 0.499932835964896, 0.598716825883862, 0.379562691207627, 0.641280588751714, 0.679202586642317, 0.477385507184256, 0.472879830723747, 0.384662859733384, 0.37112172360338, 0.703774478956202, 0.423898023856973, 0.287444757613266, 12.8375437364914, 0.755503259406695, 0.688897122172827, 0.281833163504864, 0.363521118919126, 0.443082983953355],
            [1.47210398545007, 0.767165632779031, 1.23152924484831, 0.913504624354301, 0.738415701073377, 0.96555522798098, 0.950357185031325, 0.90359652515418, 0.736731739892316, 0.443163582314639, 0.428893742551906, 0.931919140710646, 0.598558924100088, 0.439973596778216, 0.755503259406695, 3.8428474099213, 1.61392097340711, 0.385303034789668, 0.557520051020311, 0.565223766047939],
            [0.984401955935144, 0.677754679194329, 0.984152634507759, 0.694789868062839, 0.740551692220627, 0.79132074056634, 0.741425610477113, 0.579271581711225, 0.55750325361248, 0.779816109586742, 0.660328974513946, 0.792905802702891, 0.793801615982561, 0.481693682890345, 0.688897122172827, 1.61392097340711, 4.83210516236962, 0.430934143757138, 0.573156574120723, 0.980943004996173],
            [0.416548781469996, 0.395102105602751, 0.277782895565922, 0.232102315145553, 0.449983902793404, 0.509360271558287, 0.374300211820363, 0.421690355206204, 0.44408895489718, 0.408874390481926, 0.568037074030476, 0.358930070683472, 0.61029621403986, 1.37437942379832, 0.281833163504864, 0.385303034789668, 0.430934143757138, 38.1077832575802, 2.10980811550359, 0.374456331815776],
            [0.542611868915233, 0.555965424563688, 0.486030805784542, 0.345683565218847, 0.434203398141003, 0.611094097106686, 0.496467353893267, 0.348714366361603, 1.79790413031311, 0.630388930627921, 0.692059423023677, 0.532179332619096, 0.708364627674993, 2.76938062915766, 0.363521118919126, 0.557520051020311, 0.573156574120723, 2.10980811550359, 9.83220341258545, 0.658038692870502],
            [0.936458396120263, 0.42007231624752, 0.369033853043791, 0.336500141524566, 0.75584405464154, 0.466777931405709, 0.428943129877398, 0.3369549123791, 0.339447441760233, 2.41751209060932, 1.31423572845207, 0.456542719723346, 1.26893679116311, 0.745089737790822, 0.443082983953355, 0.565223766047939, 0.980943004996173, 0.374456331815776, 0.658038692870502, 3.69215640428348]]

def getSeqFile(file_in):

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
    print(seq)

    return seq

human = getSeqFile("human.py")
#oran = getSeqFile("oran.py")
# bat = getSeqFile("bat.py")
mouse = getSeqFile("mouse.py")
# rat = getSeqFile("rat.py")
# cat = getSeqFile("cat.py")
# fish = getSeqFile("fish.py")
# thaliana = getSeqFile("thaliana.py")
# corn = getSeqFile("corn.py")
#bact = getSeqFile("bact.py")
#fus = getSeqFile("fus.py")
yeast = getSeqFile("yeast.py")

seqList = [human,mouse,yeast]



def groupCompare(seqList):
    size = len(seqList)
    i =0
    seqMatrix = [([] * size) for row in range(size)]
    for orgA in seqList:

            innerList = []
            for orgB in seqList:

                innerList.append(donkeyKong(orgA,orgB))

            print(innerList, '<--- inner list')
            #seqMatrix.append(innerList)
            seqMatrix[i] = innerList
            i+=1


    print(seqMatrix)
    return seqMatrix


def donkeyKong(sequenceA, sequenceB):
    seqA = []
    seqB = []


    for acidA in sequenceA:
        seqA.append(acidSequence[acidA])

    for acidB in sequenceB:
        seqB.append(acidSequence[acidB])
    print("k3Tild", k3Tilda(seqA, seqB))
    dk = math.sqrt(2 * (1 - k3Tilda(seqA, seqB)))
    #print('dk = ', dk)
    return dk


def k3Tilda(seqA, seqB):
    final = kay3(seqA, seqB) / math.sqrt((kay3(seqA, seqA) * kay3(seqB, seqB)))
    #print("final", final)

    return final


def kay1(acidA, acidB):
    beta = 0.01
    # print('acidA:', acidA)
    # print('acidB:', acidB)
    return pow(blosum62[acidA][acidB], beta)


def kay2(sliceSeqA, sliceSeqB):
    # print('slice:', sliceSeqA, sliceSeqB)
    k2 = 1

    for i in range(len(sliceSeqA)):
        temp = kay1(sliceSeqA[i], sliceSeqB[i])

        k2 = k2 * temp

    # print("k2",k2)
    return k2


def kay3(seqA, seqB):
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
    while counter <= 100:
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
                k3 = k3 + kay2(window, longer[j:j + len(window)])
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
    groupCompare(seqList)



    
if __name__ == '__main__':
    main()