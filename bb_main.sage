import time
timing_data = open("time_data.txt","a")
timing_data.write("Started bb_main at " + str(time.strftime("%H:%M:%S, %D", time.localtime())) + "\n")
timing_data.close()

exec(open("./Code/rna_poly.py").read())
exec(open("./Code/polytope_funs.py").read())
exec(open("./Code/SuboptFile.py").read())
exec(open("./Code/NopctFile.py").read())
exec(open("./Code/accuracy_analysis.py").read())
exec(open("./Code/FilenameToAccessionNumberTools.py").read())
exec(open("./Code/bb_algorithms.py").read())
exec(open("./Code/bb_merge_alg.py").read())


seq_list = "MaxAccuracyFiles/tRNA.txt" #needs to be specified by user

## Initializes Data files
PolytopeList = InitializePolytopes(seq_list)
Slices = OrderSlices(PolytopeList)
AccuracyData = CreateAccuracyData(seq_list, PolytopeList, Slices)
NameData = CreateNameData(seq_list)
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
    IndexOrder.append(ThisIndexOrder)


timing_data = open("time_data.txt","a")
timing_data.write("Started step_0 at " + str(time.strftime("%H:%M:%S, %D", time.localtime())) + "\n")
timing_data.close()

NumIntersections = 0

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
    print("Starting sequence", x, "of Step0")
    Considering = []
    for i in range(len(Slices[x])):
        print(x,i)
        z = Slices[x][i]
        TotalBest = 0
        for j in range(0, NSequences):
            PolytopePruneList = Slices[j]
            k = 0
            NumIntersections += 1
            while (k < len(IndexOrder[j]) and z.intersection(PolytopePruneList[IndexOrder[j][k]]).is_empty()):
                NumIntersections += 1
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

    ## Updates Index Order for the Sequence that was just processed
    ## Only includes surviving regions
    ThisIndexOrder = []
    for j in Considering:
        ThisIndexOrder.append( (AccuracyData[x][j], j) )
    ## Reverse sorts greatest to smallest
    ThisIndexOrder.sort(reverse=True)
    ThisIndexOrder = [r[1] for r in ThisIndexOrder]
    IndexOrder[x] = ThisIndexOrder
    print(ThisIndexOrder)

## Prints out summative data
outfile = open("MergeData/SummativeData.txt", "a")
outfile.write("Step0\n")
outfile.write("Performed Intersections: " + str(NumIntersections) + "\n\n\n")
outfile.close()

timing_data = open("time_data.txt","a")
timing_data.write("Finished step_0 at " + str(time.strftime("%H:%M:%S, %D", time.localtime())) + "\n")
timing_data.close()

print("step_0 completed")

betterAccPrune = createAccData(PolytopeList, Slices, AccuracyData,L)

timing_data = open("time_data.txt","a")
timing_data.write("Finished AccGen at " + str(time.strftime("%H:%M:%S, %D", time.localtime())) + "\n")
timing_data.close()

print("AccGen completed")

#Start of main

x = "0"
for i in range(1, NSequences):
    x = x + "," + str(i)
Considering = [x]
## Creates queue for merge ordering so it behaves like binary tree
## and finishes with sequences in order from 0 to n-1.
cont = True
while (cont):
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
    cont = False
    for i in Considering:
        if len(i.split(",")) > 2:
            cont = True
    while (not cont and Considering[0].split(",")[0] != "0"):
        Considering.append(Considering.pop(0))

print(Considering)

## Joining to have necessary groups of two
for i in Considering:
    if (len(i.split(",")) == 2):
        x = i.split(",")[0]
        y = i.split(",")[1]
        print("Starting merge of", x, "and", y)
        z = x + "," + y
        Merge("MergeData/" + x + ".txt", "MergeData/" + y + ".txt",
              "MergeData/" + z + ".txt",
              PolytopeList, Slices, AccuracyData, L, betterAccPrune)
        

while (len(Considering) > 2):
    x = Considering.pop(0)
    if (Considering[0].split(",")[0] != "0"):
        y = Considering.pop(0)
        print("Starting merge of", x, "and", y)
        z = x + "," + y
        if (len(z.split(",")) >= 10):
            NewMerge("MergeData/" + x + ".txt", "MergeData/" + y + ".txt",
              "MergeData/" + z + ".txt",
              PolytopeList, Slices, AccuracyData, L, betterAccPrune)
        else:
            Merge("MergeData/" + x + ".txt", "MergeData/" + y + ".txt",
              "MergeData/" + z + ".txt",
              PolytopeList, Slices, AccuracyData, L, betterAccPrune)
        Considering.append(z)
    else:
        Considering.append(x)

timing_data = open("time_data.txt","a")
timing_data.write("Finished normal merges at " + str(time.strftime("%H:%M:%S, %D", time.localtime())) + "\n")
timing_data.close()

print("Starting AccSort")
AccSort("MergeData/" + Considering[0] + ".txt", "MergeData/S" + Considering[0] + ".txt", AccuracyData)
AccSort("MergeData/" + Considering[1] + ".txt", "MergeData/S" + Considering[1] + ".txt", AccuracyData)

timing_data = open("time_data.txt","a")
timing_data.write("Finished AccSort at " + str(time.strftime("%H:%M:%S, %D", time.localtime())) + "\n")
timing_data.close()

print("Starting FinalMerge")

RelIndices1 = [eval(i) for i in Considering[0].split(",")]
RelIndices2 = [eval(i) for i in Considering[1].split(",")]

FinalMerge("MergeData/S" + Considering[0] + ".txt", "MergeData/S" + Considering[1] + ".txt",
           "MergeData/FinalOutput.txt", PolytopeList, Slices, AccuracyData, NameData, L,
           "Lfile.txt", RelIndices1, RelIndices2)

timing_data = open("time_data.txt","a")
timing_data.write("Finished bb_main at " + str(time.strftime("%H:%M:%S, %D", time.localtime())))
timing_data.close()
