from bb_algorithms import *
from bb_merge_alg import *
exec(open("./Code/rna_poly.py").read())


seq_list = "Max_Accuracy.txt" #needs to be changed

## Initializes Data files
PolytopeList = InitializePolytopes(seq_list)
Slices = OrderSlices(PolytopeList)
AccuracyData = CreateAccuracyData(seq_list, PolytopeList, Slices)
NameData = CreateNameData(seq_list)
u = FindBounds(seq_list) # do we need this? Seems like not
NSequences = len(PolytopeList)

## Imports the L value for lower bound on accuracy
## The code assumes that there is an identical Lfile.txt in the same directory 
L = getLValue("L.txt")

## Creates a double array
## IndexOrder[x] is the index order for sequence x
## IndexOrder[x] lists the polytope indexes for x from highest accuracy to lowest accuracy
IndexOrder = []
for i in range(0,NSequences):
    ThisIndexOrder = []
    for j in range(len(AccuracyData[i])):
        ThisIndexOrder.append( (AccuracyData[i][j], j) )
    ThisIndexOrder.sort(reverse=True)
    ThisIndexOrder = [r[1] for r in ThisIndexOrder]
    print(ThisIndexOrder)
    IndexOrder.append(ThisIndexOrder)

## Consider each region from each sequence
## Consider all regions it intersects from each sequence
## Take the highest accuracy value from the regions it intersects from each sequence and sum them together
## This sum must be greater than L * NSequences otherwise we know that any combination including it
## Will be thrown out eventually for a bad accuracy value
## Thereby, we can reduce the regions we are considering for each sequence from ~500 down to
## about ~15, drastically increasing the runtime
## In addition, using IndexOrder created above, we can check the regions it intersects from highest accuracy
## to lowest accuracy in order to avoid having to check all regions from each intersecting sequence.
## This outputs the regions still in consideration for sequence x to "MergeData/x.txt"
for x in range(0,NSequences):
    Considering = []
    for i in range(len(Slices[x])):
        print(x,i)
        z = Slices[x][i]
        TotalBest = 0
        for j in range(0, NSequences):
            PolytopePruneList = Slices[j]
            k = 0
            while (k < len(IndexOrder[j]) and z.intersection(PolytopePruneList[IndexOrder[j][k]]).is_empty()):
                k += 1
            if (k >= len(IndexOrder[j])):
                break
            TotalBest += AccuracyData[j][IndexOrder[j][k]]
        if (TotalBest > L * NSequences):
            Considering.append(i)

    print(Considering)
    outfile = open("MergeData/" + str(x) + ".txt", "w")
    for i in range(len(Considering)):
        outfile.write("[" + str(Considering[i]) + "]" + "\n")
    outfile.close()

    ThisIndexOrder = []
    for j in Considering:
        ThisIndexOrder.append( (AccuracyData[x][j], j) )
    ## Reverse sorts greatest to smallest
    ThisIndexOrder.sort(reverse=True)
    ThisIndexOrder = [r[1] for r in ThisIndexOrder]
    IndexOrder[x] = ThisIndexOrder
    print(ThisIndexOrder)



betterAccPrune = createAccData(PolytopeList, Slices, AccuracyData,L)

## Start of main

x = "0"
for i in range(1, NSequences):
    x = x + "," + str(i)
Considering = [x]
while (len(Considering[0].split(",")) > 2):
    x = Considering.pop(0)
    x = x.split(",")
    Len = (len(x) + 1) // 2
    y = x[0]
    for i in range(1, Len):
        y = y + "," + x[i]
    z = x[Len]
    for i in range(Len+1,len(x)):
        z = z + "," + x[i]
    Considering.append(y)
    Considering.append(z)

## Joining to have necessary groups of two
for i in Considering:
    if (len(i.split(",")) == 2):
        x = i.split(",")[0]
        y = i.split(",")[1]
        z = x + "," + y
        Merge("MergeData/" + x + ".txt", "MergeData/" + y + ".txt",
              "MergeData/" + z + ".txt",
              PolytopeList, Slices, AccuracyData, L, u, betterAccPrune)
        

while (len(Considering) > 8):
    x = Considering.pop(0)
    y = Considering.pop(0)
    z = x + "," + y
    Merge("MergeData/" + x + ".txt", "MergeData/" + y + ".txt",
          "MergeData/" + z + ".txt",
          PolytopeList, Slices, AccuracyData, L, u, betterAccPrune)
    Considering.append(z)

AccSort("MergeData/" + Considering[0] + ".txt", "MergeData/S" + Considering[0] + ".txt", AccuracyData)
AccSort("MergeData/" + Considering[1] + ".txt", "MergeData/S" + Considering[1] + ".txt", AccuracyData)


## indication of which two branches are merged in the last step
## if run on 100 sequences
x="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49"
y="50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99"
RelevantIndices1 = [50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99]
RRelevantIndices2 = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49]

##if run on 50 sequences
## x="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24"
## y ="25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49"
## RelevantIndices1 = [25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49]
## RelevantIndices2 = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]

FinalMerge( "MergeData/FinalInput.txt", "MergeData/S" + x + ".txt",  "MergeData/FinalOutput.txt", PolytopeList, Slices, AccuracyData, NameData, L, "Lfile.txt", RelevantIndices1, RelevantIndices2)

           
