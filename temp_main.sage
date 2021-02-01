import time

exec(open("./Code/rna_poly.py").read())
exec(open("./Code/polytope_funs.py").read())
exec(open("./Code/SuboptFile.py").read())
exec(open("./Code/NopctFile.py").read())
exec(open("./Code/accuracy_analysis.py").read())
exec(open("./Code/FilenameToAccessionNumberTools.py").read())
exec(open("./Code/bb_algorithms.py").read())
exec(open("./Code/bb_merge_alg.py").read())


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

while (len(Considering) > 2):
    print("Starting merge of", x, "and", y)
    x = Considering.pop(0)
    y = Considering.pop(0)
    z = x + "," + y
    Considering.append(z)

timing_data = open("time_data.txt","a")
timing_data.write("Started FinalMerge at " + str(time.strftime("%H:%M:%S, %D", time.localtime())) + "\n")
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
