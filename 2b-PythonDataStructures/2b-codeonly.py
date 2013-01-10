voltageList = [-2.0, -1.0, 0.0, 1.0, 2.0]
# <codecell>
currentList = [-1.0, -0.5, 0.0, 0.5, 1.0]
# <codecell>
type(voltageList)
# <codecell>
voltageList[0]
# <codecell>
voltageList[2]
# <codecell>
currentList[-1]
# <codecell>
currentList[-2]
# <codecell>
voltageList[1:4]
# <codecell>
voltageList[2:]
# <codecell>
dir(list)
# <codecell>
voltageList.append(3.)
# <codecell>
voltageList.append(4.)
# <codecell>
voltageList
# <codecell>
currentList.extend([1.5, 2.0])
# <codecell>
currentList
# <codecell>
len(voltageList)
# <codecell>
dataList = ["experiment: current vs. voltage", \
# <codecell>
a = [1,2]
# <codecell>
b = a
# <codecell>
a.append(10)
# <codecell>
b
# <codecell>
f = open("data.dat")
# <codecell>
ivdata = f.readlines()
# <codecell>
f.close()
# <codecell>
ivdata
# <codecell>
tup = ("red", "white", "blue")
# <codecell>
type(tup)
# <codecell>
fruit = set(["apple", "banana", "pear", "banana"]) # You have to use a list to create a set.
# <codecell>
firstBowl = set(["apple", "banana", "pear", "peach"])
# <codecell>
secondBowl = set(["peach", "watermelon", "orange", "apple"])
# <codecell>
set.intersection(firstBowl, secondBowl)
# <codecell>
dataDict = {"experiment": "current vs. voltage", \
# <codecell>
dataDict["run"]
# <codecell>
dataDict["voltage"]
# <codecell>
dataDict["current"][-1]
# <codecell>
dataDict["temperature"] = 3275.39
# <codecell>
dataDict["user"] = "Johann G. von Ulm"
# <codecell>
dataDict.keys()
# <codecell>
dataDict.values()
