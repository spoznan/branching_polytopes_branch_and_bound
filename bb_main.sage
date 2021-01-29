import time
timing_data = open("time_data.txt","a")
timing_data.write("Started bb_main at " + str(time.strftime("%H:%M:%S, %D", time.localtime())))
timing_data.close()

exec(open("./Code/bb_algorithms.py").read())
exec(open("./Code/bb_merge_alg.py").read())
exec(open("./Code/rna_poly.py").read())

seq_list = "MaxAccuracyFiles/tRNA.txt" #needs to be changed

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
        

while (len(Considering) > 2):
    x = Considering.pop(0)
    y = Considering.pop(0)
    z = x + "," + y
    if (len(z.split(",")) >= 10):
        NewMerge("MergeData/" + x + ".txt", "MergeData/" + y + ".txt",
          "MergeData/" + z + ".txt",
          PolytopeList, Slices, AccuracyData, L, u, betterAccPrune)
    else:
        Merge("MergeData/" + x + ".txt", "MergeData/" + y + ".txt",
          "MergeData/" + z + ".txt",
          PolytopeList, Slices, AccuracyData, L, u, betterAccPrune)
    Considering.append(z)

AccSort("MergeData/" + Considering[0] + ".txt", "MergeData/S" + Considering[0] + ".txt", AccuracyData)
AccSort("MergeData/" + Considering[1] + ".txt", "MergeData/S" + Considering[1] + ".txt", AccuracyData)

FinalMerge("MergeData/S" + Considering[0] + ".txt", "MergeData/S" + Considering[1] + ".txt",
           "MergeData/FinalOutput.txt", PolytopeList, Slices, AccuracyData, NameData, L)

timing_data = open("time_data.txt","a")
timing_data.write("Finished bb_main at " + str(time.strftime("%H:%M:%S, %D", time.localtime())))
timing_data.close()
