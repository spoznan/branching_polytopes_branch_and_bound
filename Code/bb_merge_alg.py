## Takes infile1name and infile2name and merges them into outfilename
## Takes L and u as arguments
## infile1name is of the form "MergeData/x,y,...,z.txt" and has the considering combinations of
## polytopes from sequences x,y,..., and z
## infile2name is of the same form.
## outfilename is designed to be the combination of them

def Merge(infile1name, infile2name, outfilename, PolytopeList, Slices, AccuracyData, L, betterAccPrune):
    NSequences = len(PolytopeList)
    ThrownOutForAccuracy = 0
    ThrownOutForNonIntersect = 0
    NumIntersections = 0

    ## Reads in raw data
    infile1 = open(infile1name,"r")
    infile1data = infile1.readlines()
    infile1.close()
    ConsideringIntersections1 = []
    
    infile2 = open(infile2name,"r")
    infile2data = infile2.readlines()
    infile2.close()
    ConsideringIntersections2 = []

    ## Converts raw data to a list
    for i in infile1data:
        ConsideringIntersections1.append(eval(i.strip()))
    for i in infile2data:
        ConsideringIntersections2.append(eval(i.strip()))

    ## Finds relevant indices for each set
    infile1nameprocess = infile1name.split(".")[0].split("/")[1]
    RelevantIndices1 = infile1nameprocess.split(",")
    for i in range(len(RelevantIndices1)):
        RelevantIndices1[i] = eval(RelevantIndices1[i])
    
    infile2nameprocess = infile2name.split(".")[0].split("/")[1]
    RelevantIndices2 = infile2nameprocess.split(",")
    for i in range(len(RelevantIndices2)):
        RelevantIndices2[i] = eval(RelevantIndices2[i])

    ## Creates Polytope Data
    PolytopeList1 = []
    for i in range(len(ConsideringIntersections1)):
        z = Slices[RelevantIndices1[0]][ConsideringIntersections1[i][0]]
        for j in range(1,len(ConsideringIntersections1[i])):
            z = z.intersection(Slices[RelevantIndices1[j]][ConsideringIntersections1[i][j]])
            NumIntersections += 1
        PolytopeList1.append(z)

    PolytopeList2 = []
    for i in range(len(ConsideringIntersections2)):
        z = Slices[RelevantIndices2[0]][ConsideringIntersections2[i][0]]
        for j in range(1,len(ConsideringIntersections2[i])):
            z = z.intersection(Slices[RelevantIndices2[j]][ConsideringIntersections2[i][j]])
            NumIntersections += 1
        PolytopeList2.append(z)

    ## Calculates Main
    FinalIntersections = []
    for x in range(len(PolytopeList1)):
        if (x % 100 == 0):
            print("Starting region", x, "of", len(PolytopeList1))
        for y in range(len(PolytopeList2)):
            ## Checks to see if the two regions being considered intersect
            NumIntersections += 1
            if (not PolytopeList1[x].intersection(PolytopeList2[y]).is_empty()):
                ## Starts accuracy sum for the best possible accuracy for all polytopes not currently being considered
                AccSum = 0
                for i in range(NSequences):
                    if (RelevantIndices1.count(i) == 0 and RelevantIndices2.count(i) == 0):
                        Min1 = min(betterAccPrune[ SN ][ ConsideringIntersections1[x][RelevantIndices1.index(SN)]][ i ] for SN in RelevantIndices1)
                        Min2 = min(betterAccPrune[ SN ][ ConsideringIntersections2[y][RelevantIndices2.index(SN)]][ i ] for SN in RelevantIndices2)
                        AccSum += min(Min1, Min2)
                ## Adds the accuracy of the 1st region
                for i in range(len(RelevantIndices1)):
                    AccSum += AccuracyData[RelevantIndices1[i]][ConsideringIntersections1[x][i]]
                ## Adds the accuracy of the 2nd region
                for i in range(len(RelevantIndices2)):
                    AccSum += AccuracyData[RelevantIndices2[i]][ConsideringIntersections2[y][i]]
                ## Only continues to consider the region if it meets the accuracy threshold
                if (AccSum > L * NSequences):
                    FinalIntersections.append(ConsideringIntersections1[x] + ConsideringIntersections2[y])
                else:
                    ThrownOutForAccuracy += 1
            else:
                ThrownOutForNonIntersect += 1
    
    ## Prints out combined data
    outfile = open(outfilename,"w")
    for i in FinalIntersections:
        outfile.write(str(i))
        outfile.write("\n")
    outfile.close()

    ## Prints out summative data
    outfile = open("MergeData/SummativeData.txt", "a")
    outfile.write(infile1name + " and " + infile2name + "\n")
    outfile.write("Thrown out for Accuracy: " + str(ThrownOutForAccuracy) + "\n")
    outfile.write("Thrown out for Non-Intersect: " + str(ThrownOutForNonIntersect) + "\n")
    outfile.write("Performed Intersections: " + str(NumIntersections) + "\n")
    outfile.write("Length of join: " + str(len(FinalIntersections)) + "\n\n\n")
    outfile.close()
    
    return

def NewMerge(infile1name, infile2name, outfilename, PolytopeList, Slices, AccuracyData, L, betterAccPrune):
    NSequences = len(PolytopeList)
    ThrownOutForAccuracy = 0
    ThrownOutForNonIntersect = 0
    NumIntersections = 0

    ## Reads in raw data
    infile1 = open(infile1name,"r")
    infile1data = infile1.readlines()
    infile1.close()
    ConsideringIntersections1 = []
    
    infile2 = open(infile2name,"r")
    infile2data = infile2.readlines()
    infile2.close()
    ConsideringIntersections2 = []

    ## Converts raw data to a list
    for i in infile1data:
        ConsideringIntersections1.append(eval(i.strip()))
    for i in infile2data:
        ConsideringIntersections2.append(eval(i.strip()))

    ## Finds relevant indices for each set
    infile1nameprocess = infile1name.split(".")[0].split("/")[1]
    RelevantIndices1 = infile1nameprocess.split(",")
    for i in range(len(RelevantIndices1)):
        RelevantIndices1[i] = eval(RelevantIndices1[i])
    
    infile2nameprocess = infile2name.split(".")[0].split("/")[1]
    RelevantIndices2 = infile2nameprocess.split(",")
    for i in range(len(RelevantIndices2)):
        RelevantIndices2[i] = eval(RelevantIndices2[i])

    ARegions = []
    for i in range(len(ConsideringIntersections1[0])):
        ARegions.append([])
        for j in ConsideringIntersections1:
            if (not j[i] in ARegions[i]):
                ARegions[i].append(j[i])
    BRegions = []
    for i in range(len(ConsideringIntersections2[0])):
        BRegions.append([])
        for j in ConsideringIntersections2:
            if (not j[i] in BRegions[i]):
                BRegions[i].append(j[i])

    NTrue = 0 #Empty
    NFalse = 0 #Non-empty
    ## Empty reads true
    ## Non-empty reads false
    global IntersectionMatrix
    IntersectionMatrix = []
    for m in range(len(ARegions)): # Sequence Number (m) of File A
        IntersectionMatrix.append({})
        for i in ARegions[m]: # Region Number (i) of Sequence m of File A
            z = Slices[RelevantIndices1[m]][i]
            IntersectionMatrix[m][i] = []
            for j in range(len(BRegions)): # Sequence Number (j) of File B
                IntersectionMatrix[m][i].append({})
                for k in BRegions[j]: # Region Number (k) of Sequence j of File B
                    NumIntersections += 1
                    IntersectionMatrix[m][i][j][k] = z.intersection(Slices[RelevantIndices2[j]][k]).is_empty()
                    if (IntersectionMatrix[m][i][j][k]):
                        NTrue += 1
                    else:
                        NFalse += 1

    ## Calculates Main
    FinalIntersections = []
    for x in range(len(ConsideringIntersections1)):
        if (x % 100 == 0):
            print("Starting region", x, "of", len(ConsideringIntersections1))
        for y in range(len(ConsideringIntersections2)):
            ## Checks to see if the two regions being considered intersect
            if (intersects(ConsideringIntersections1[x], ConsideringIntersections2[y])):
                z = Slices[RelevantIndices1[0]][ConsideringIntersections1[x][0]]
                for k in range(1, len(ConsideringIntersections1[x])):
                    NumIntersections += 1
                    z = z.intersection(Slices[RelevantIndices1[k]][ConsideringIntersections1[x][k]])
                for k in range(0, len(ConsideringIntersections2[y])):
                    NumIntersections += 1
                    z = z.intersection(Slices[RelevantIndices2[k]][ConsideringIntersections2[y][k]])
                if (not z.is_empty()):
                    ## Starts accuracy sum for the best possible accuracy for all polytopes not currently being considered
                    AccSum = 0
                    for i in range(NSequences):
                        if (RelevantIndices1.count(i) == 0 and RelevantIndices2.count(i) == 0):
                            Min1 = min(betterAccPrune[ SN ][ ConsideringIntersections1[x][RelevantIndices1.index(SN)]][ i ] for SN in RelevantIndices1)
                            Min2 = min(betterAccPrune[ SN ][ ConsideringIntersections2[y][RelevantIndices2.index(SN)]][ i ] for SN in RelevantIndices2)
                            AccSum += min(Min1, Min2)
                    ## Adds the accuracy of the 1st region
                    for i in range(len(RelevantIndices1)):
                        AccSum += AccuracyData[RelevantIndices1[i]][ConsideringIntersections1[x][i]]
                    ## Adds the accuracy of the 2nd region
                    for i in range(len(RelevantIndices2)):
                        AccSum += AccuracyData[RelevantIndices2[i]][ConsideringIntersections2[y][i]]
                    ## Only continues to consider the region if it meets the accuracy threshold
                    if (AccSum > L * NSequences):
                        FinalIntersections.append(ConsideringIntersections1[x] + ConsideringIntersections2[y])
                    else:
                        ThrownOutForAccuracy += 1
                else:
                    ThrownOutForNonIntersect += 1
            else:
                ThrownOutForNonIntersect += 1
    
    ## Prints out combined data
    outfile = open(outfilename,"w")
    for i in FinalIntersections:
        outfile.write(str(i))
        outfile.write("\n")
    outfile.close()

    ## Prints out summative data
    outfile = open("MergeData/SummativeData.txt", "a")
    outfile.write(infile1name + " and " + infile2name + "\n")
    outfile.write("Thrown out for Accuracy: " + str(ThrownOutForAccuracy) + "\n")
    outfile.write("Thrown out for Non-Intersect: " + str(ThrownOutForNonIntersect) + "\n")
    outfile.write("Performed Intersections: " + str(NumIntersections) + "\n")
    outfile.write("Length of join: " + str(len(FinalIntersections)) + "\n\n\n")
    outfile.close()
    
    return

def AccSort(infilename, outfilename, AccuracyData):
    infile = open(infilename, "r")
    infiledata = infile.readlines()
    infile.close()
    ConsideringIntersections = []
    for i in infiledata:
        ConsideringIntersections.append(eval(i.strip()))
    infilenameprocess = infilename.split(".")[0].split("/")[1]
    RelevantIndices = infilenameprocess.split(",")
    for i in range(len(RelevantIndices)):
        RelevantIndices[i] = eval(RelevantIndices[i])
    AccValue = []
    for i in range(len(ConsideringIntersections)):
        AccValue.append(0)
        for j in range(len(ConsideringIntersections[i])):
            AccValue[i] += AccuracyData[RelevantIndices[j]][ConsideringIntersections[i][j]]
        ConsideringIntersections[i] = (AccValue[i], ConsideringIntersections[i])
    ConsideringIntersections.sort(reverse=True)
    outfile = open(outfilename, "w")
    for i in ConsideringIntersections:
        outfile.write(str(i[0]) + " = " + str(i[1]) + "\n")
    outfile.close()
    return

def intersects(x,y):
    global IntersectionMatrix
    for i in range(len(x)):
        for j in range(len(y)):
            if (IntersectionMatrix[i][x[i]][j][y[j]]):
                return False
    return True

    
def FinalMerge(infile1name, infile2name, outfilename, PolytopeList, Slices, AccuracyData, NameData, L, Lfile, RelevantIndices1, RelevantIndices2):
    NSequences = len(PolytopeList)
    NumIntersections = 0

    ## Reads in raw data
    infile1 = open(infile1name,"r")
    infile1data = infile1.readlines()
    infile1.close()
    ConsideringIntersections1 = []
    AccuracyIntersections1 = []
    
    infile2 = open(infile2name,"r")
    infile2data = infile2.readlines()
    infile2.close()
    ConsideringIntersections2 = []
    AccuracyIntersections2 = []

    ## Converts raw data to a list
    for i in infile1data:
        ConsideringIntersections1.append(eval(i.strip().split(" = ")[1]))
    print("ConsideringIntersections 1 generated")
    for i in infile2data:
        ConsideringIntersections2.append(eval(i.strip().split(" = ")[1]))
    print("ConsideringIntersections 2 generated")
    
    for i in infile1data:
        AccuracyIntersections1.append(eval(i.strip().split(" = ")[0]))
    print("AccuracyIntersections 1 generated")
    for i in infile2data:
        AccuracyIntersections2.append(eval(i.strip().split(" = ")[0]))
    print("AccuracyIntersections 2 generated")

    ARegions = []
    for i in range(len(ConsideringIntersections1[0])):
        ARegions.append([])
        for j in ConsideringIntersections1:
            if (not j[i] in ARegions[i]):
                ARegions[i].append(j[i])
    
    BRegions = []
    for i in range(len(ConsideringIntersections2[0])):
        BRegions.append([])
        for j in ConsideringIntersections2:
            if (not j[i] in BRegions[i]):
                BRegions[i].append(j[i])
    
    NTrue = 0 #Empty
    NFalse = 0 #Non-empty
    ## Empty reads true
    ## Non-empty reads false
    global IntersectionMatrix
    IntersectionMatrix = []
    for m in range(len(ARegions)): # Sequence Number (m) of File A
        IntersectionMatrix.append({})
        for i in ARegions[m]: # Region Number (i) of Sequence m of File A
            z = Slices[RelevantIndices1[m]][i]
            IntersectionMatrix[m][i] = []
            for j in range(len(BRegions)): # Sequence Number (j) of File B
                IntersectionMatrix[m][i].append({})
                for k in BRegions[j]: # Region Number (k) of Sequence j of File B
                    NumIntersections += 1
                    IntersectionMatrix[m][i][j][k] = z.intersection(Slices[RelevantIndices2[j]][k]).is_empty()
                    if (IntersectionMatrix[m][i][j][k]):
                        NTrue += 1
                    else:
                        NFalse += 1
    print("Number Pairwise Empty:", NTrue)
    print("Number Pariwise Non-Empty:", NFalse)

    i = 0
    needsUpdate = True
    while (i < len(infile1data)):
        #x = infile1data[i].strip().split(" = ")[1]
        #y = infile1data[i].strip().split(" = ")[0]
        ## Runtime creating of PolytopeList1 rather than pre-initialization
        #z = Slices[RelevantIndices1[0]][ConsideringIntersections1[i][0]]
        #for j in range(1,len(ConsideringIntersections1[i])):
        #    z = z.intersection(Slices[RelevantIndices1[j]][ConsideringIntersections1[i][j]])

        ## L update function

        if (i % 100 == 0):
            needsUpdate = True
            print("Starting region", i, "of", len(infile1data))
            
        if (needsUpdate):
            infile = open(Lfile, "r")
            raw = infile.readlines()
            infile.close()
            try:
                for k in range(len(raw)):
                    raw[k] = eval(raw[k].strip())
                currentL = max(raw)
            except ValueError:
                print("Updated L lookup failed")
                currentL = L
            if (L < currentL):
                print("Received new L of", currentL)
                L = currentL
            elif (L > currentL):
                print("Uploaded new L of", L)
                outfile = open(Lfile, "a")
                outfile.write(str(L) + "\n")
                outfile.close()

        needsUpdate = False

        j = 0
        while (j < len(ConsideringIntersections2)):
            if (AccuracyIntersections1[i] + AccuracyIntersections2[j] < L * NSequences):
                break
            if (intersects(ConsideringIntersections1[i], ConsideringIntersections2[j])):
                ## Check if intersection is non-empty
                z = Slices[RelevantIndices1[0]][ConsideringIntersections1[i][0]]
                for k in range(1, len(ConsideringIntersections1[i])):
                    NumIntersections += 1
                    z = z.intersection(Slices[RelevantIndices1[k]][ConsideringIntersections1[i][k]])
                for k in range(0, len(ConsideringIntersections2[j])):
                    NumIntersections += 1
                    z = z.intersection(Slices[RelevantIndices2[k]][ConsideringIntersections2[j][k]])
                if (not z.is_empty()):
                    L = (AccuracyIntersections1[i] + AccuracyIntersections2[j]) / NSequences
                    print(L)
                    CurrentBest = (ConsideringIntersections1[i], ConsideringIntersections2[j], i, j)
                    print(CurrentBest)
                    needsUpdate = True
                    break
            j += 1
        i += 1

    ## Prints out summative data
    outfile = open("MergeData/SummativeData.txt", "a")
    outfile.write("Final Merge\n")
    outfile.write("Performed Intersections: " + str(NumIntersections) + "\n\n\n")
    outfile.close()
    
    outfile = open(outfilename, "w")
    outfile.write("Accuracy: "  + str(L) + "\n")
    outfile.write(str(CurrentBest) + "\n")
    R = CurrentBest[0] + CurrentBest[1]
    S = RelevantIndices1 + RelevantIndices2
    Region = Polyhedron(lines = [(1,0,0),(0,1,0),(0,0,1)])
    for i in range(0,NSequences):
        d1Slice = Slices[S[i]][R[i]]
        Region = Region.intersection(d1Slice)
        polyname = NameData[S[i]]
        print("Polytope: ", polyname)
        outfile.write("Polytope: " + polyname +"\n")
        sl_center = slice_center(d1Slice)
        print("slice center: ", sl_center)
        outfile.write("slice center: " + str(sl_center) + "\n")
        vert = [vector(i) for i in d1Slice.vertices()]
        rays = [vector(i) for i in d1Slice.rays()]
        acc = AccuracyOf(d1Slice, polyname)
        print("Accuracy: ", acc)
        outfile.write("Accuracy: " + str(acc) + "\n")
        print("Slice vertices: ", vert)
        outfile.write("Polyhedron( vertices = " + str(vert) + ", rays = " + str(rays) +")\n")
        print("Slice rays: ", rays)
        rnasuboptfile = polyname
        try:
            rnasuboptfile += str(slice_center(d1Slice)[0]).split("/")[0] + ":" + str(slice_center(d1Slice)[0]).split("/")[1] + "_"
        except:
            rnasuboptfile += str(slice_center(d1Slice)[0]) + "_"
        try:
            rnasuboptfile += str(slice_center(d1Slice)[1]).split("/")[0] + ":" + str(slice_center(d1Slice)[1]).split("/")[1] + "_"
        except:
            rnasuboptfile += str(slice_center(d1Slice)[1]) + "_"
        try:
            rnasuboptfile += str(slice_center(d1Slice)[2]).split("/")[0] + ":" + str(slice_center(d1Slice)[2]).split("/")[1]
        except:
            rnasuboptfile += str(slice_center(d1Slice)[2])
        rnasuboptfile += ".rnasubopt"
        print(rnasuboptfile)
        outfile.write(rnasuboptfile + "\n")
    print("Center of region: ", Region.center())
    print("Rays: ", Region.rays())
    print("Vertices: ", Region.vertices())
    outfile.write("Center of region:" + str(Region.center()) + "\n")
    outfile.write("Vertices: " + str(Region.vertices()) + "\n")
    outfile.write("Rays: " + str(Region.rays()) + "\n")
    outfile.close()
    return

def createAccData(PolytopeList, Slices, AccuracyData, L):
    NSequences = len(PolytopeList)
    NumIntersections = 0
    
    IndexOrder = []
    for i in range(0,NSequences):
        ThisIndexOrder = []
        for j in range(len(AccuracyData[i])):
            ThisIndexOrder.append( (AccuracyData[i][j], j) )
        ThisIndexOrder.sort(reverse=True)
        ThisIndexOrder = [r[1] for r in ThisIndexOrder]
        IndexOrder.append(ThisIndexOrder)
        
    BetterAccPrune = []
    for x in range(0,NSequences):
        print("Starting sequence", x, "of AccGen")
        BetterAccPrune.append({})
        infile = open("MergeData/" + str(x) + ".txt", "r")
        Considering = infile.readlines()
        infile.close()
        Considering = [eval(r.strip())[0] for r in Considering]
        for i in Considering:
            NewList = []
            z = Slices[x][i]
            for j in range(0, NSequences):
                PolytopePruneList = Slices[j]
                k = 0
                NumIntersections += 1
                while (z.intersection(PolytopePruneList[IndexOrder[j][k]]).is_empty()):
                    k += 1
                    NumIntersections += 1
                NewList.append(AccuracyData[j][IndexOrder[j][k]])
            BetterAccPrune[x][i] = NewList
            outfile = open("AccData/" + str(x) + "," + str(i) + ".txt", "w")
            for j in BetterAccPrune[x][i]:
                outfile.write(str(j) + "\n")
            outfile.close()

    ## Prints out summative data
    outfile = open("MergeData/SummativeData.txt", "a")
    outfile.write("AccGen\n")
    outfile.write("Performed Intersections: " + str(NumIntersections) + "\n\n\n")
    outfile.close()
    
    return(BetterAccPrune)

def readAccData(PolytopeList, Slices, AccuracyData, L):
    NSequences = len(PolytopeList)

    IndexOrder = []
    for i in range(0,NSequences):
        ThisIndexOrder = []
        for j in range(len(AccuracyData[i])):
            ThisIndexOrder.append( (AccuracyData[i][j], j) )
        ThisIndexOrder.sort(reverse=True)
        ThisIndexOrder = [r[1] for r in ThisIndexOrder]
        IndexOrder.append(ThisIndexOrder)

    BetterAccPruneData = []
    for x in range(0,NSequences):
        BetterAccPruneData.append({})
        infile = open("MergeData/" + str(x) + ".txt", "r")
        Considering = infile.readlines()
        infile.close()
        Considering = [eval(r.strip())[0] for r in Considering]
        for i in Considering:
            accfile = open("AccData/" + str(x) + "," + str(i) + ".txt", "r")
            NewList =  accfile.readlines()
            NewList = [eval(r.strip()) for r in NewList]
            accfile.close()
            BetterAccPruneData[x][i] = NewList
    return(BetterAccPruneData)
