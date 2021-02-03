infile = open("SplitTest.txt", "r")
raw = infile.readlines()
infile.close()
filenames = [i.split(" = ")[0] for i in raw]

infile = open("5S_50_results_paper.csv", "r")
data = infile.readlines()[1:]
infile.close()
infile = open("tRNA_50_results_paper.csv", "r")
data = data + infile.readlines()[1:]
infile.close()
i = 0
while (i < len(data)):
    if (data[i] == "\n"):
        data.pop(i)
    else:
        i += 1
x = [0] * len(data[0].split(","))
for i in data:
    if filenames.count(i.split(",")[0]) > 0:
        for j in range(1, len(i.split(","))):
            x[j] += eval(i.split(",")[j])

print(x)
print(max(x))
