"""
Protein Sequencing Project
Name:
Roll Number:
"""

import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    f=open(filename)
    x=f.read().splitlines()
    string="".join(x)
    # print(string)
    # print(len(string))
    # print(len(x))
    return string


'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    lst=[]
    for i in range(startIndex,len(dna),3):
        lst.append(dna[i:i+3])
        if dna[i:i+3]== 'TAG' or dna[i:i+3] == 'TAA' or dna[i:i+3] == 'TGA' :
            break
    replace=[string.replace("T","U")for string in lst]
    return replace


'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    import json
    d={}
    f=open(filename)
    x=json.load(f)
    for k,v in x.items():
        for i in v:
            d[i.replace("T","U")]=k
    return d


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    # print("this is", codons)
    # print("this is", codonD)
    empty=[]
    if codons[0] == "AUG":
        empty.append("Start")
    for i in range(1,len(codons)):
        if codons[i] in codonD.keys():
            empty.append(codonD[codons[i]])
    # print("this is" ,empty)
    return empty


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    x=readFile(dnaFilename)
    y=makeCodonDictionary(codonFilename)
    count=0
    lst=[]
    i=0
    while i < len(x):
        if x[i:i+3]=="ATG":
            z=dnaToRna(x, i)
            # print(z)      
            protein=generateProtein(z,y)
            # print(z)
            lst.append(protein)
            i=i+3*len(z)
        else:
            i=i+1
            count+=1
    # print(lst)
    return lst


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    lst=[]
    for row in proteinList1:
        for col in proteinList2:
            if row==col and row not in lst:
                lst.append(row)
    return lst


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    lst=[]
    for i in proteinList:
        for j in i:
            lst.append(j)
    return lst


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    dict={}
    for i in aaList:
        if i not in dict:
            dict[i]=1
        else:
            dict[i]+=1
    return dict


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    lst1=combineProteins(proteinList1)
    lst2=combineProteins(proteinList2)
    dict1=aminoAcidDictionary(lst1)
    dict2=aminoAcidDictionary(lst2)
    freq1={}
    freq2={}
    templst=[]
    finddiff=[]
    for i in dict1:
        freq1[i]=dict1[i]/len(lst1)
        # print(freq1[i])
        if i not in templst and i!='Start' and i!='Stop':
            templst.append(i)
    for j in dict2:
        freq2[j]=dict2[j]/len(lst2)
        if j not in templst and j!='Start' and j!='Stop':
            templst.append(j)
    # print(temp)
    for k in templst:
        frequency1=0
        frequency2=0
        if k in freq1:
            frequency1=freq1[k]
        if k in freq2:
            frequency2=freq2[k]
        difference=frequency2-frequency1
        if difference < -cutoff or difference > cutoff:
            #difflst=[k,frequency1,frequency2]
            finddiff.append([k,frequency1,frequency2])
    return finddiff



'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    print("The following proteins occurred in both DNA Sequences:")
    for i in commonalities:
        count=0
        proteins=""
        lst=i[1:len(i)-1]
        for j in lst:
            proteins+=j
            count+=1
            if count!=len(lst):
                proteins+="-"
        if len(proteins)!=0:
            print(proteins)
    print("The following amino acids occurred at very different rates in the two DNA sequences:")
    for i in differences:
        x=i[0]
        frequency1=round(i[1]*100,2)
        frequency2=round(i[2]*100,2)
        print(str(x)+" "+str(frequency1)+" % in Seq1"+","+str(frequency2)+"% in Seq2")
    return



def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    x=combineProteins(proteinList1)
    y=combineProteins(proteinList2)
    lst1=aminoAcidDictionary(x)
    lst2=aminoAcidDictionary(y)
    lst=[]
    for i in lst1:
        if i not in lst:
            lst.append(i)
    for k in lst2:
        if k not in lst:
            lst.append(k)
    lst.sort()
    return lst


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    x=combineProteins(proteinList)
    y=aminoAcidDictionary(x)
    lst=[]
    for i in labels:
        if i in y:
            z=y[i]/len(x)
            lst.append(z)
        else:
            lst.append(0)
    return lst


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    import numpy as np
    w=0.4
    xvalues=np.arange(len(xLabels))
    plt.bar(xvalues,freqList1,width=-w,align='edge',label=label1,edgecolor=edgeList)
    plt.bar(xvalues,freqList2,width=w,align='edge',label=label2,edgecolor=edgeList)
    plt.xticks(ticks=list(range(len(xLabels))),labels=xLabels,rotation="horizontal")
    plt.legend()
    plt.title("Comparision of Frequencies")
    plt.show()
    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    # print(labels)
    # print(biggestDiffs)
    lst=[]
    words=[]
    for i in range(len(biggestDiffs)):
        words.append(biggestDiffs[i][0])
    for i in range(len(labels)):
        if labels[i] in words:
            lst.append("black")
        else:
            lst.append("white")
    # print(lst)
    return lst


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    humanproteins=synthesizeProteins("data/human_p53.txt","data/codon_table.json")
    elephantproteins=synthesizeProteins("data/elephant_p53.txt","data/codon_table.json")
    x=commonProteins(humanproteins,elephantproteins)
    y=findAminoAcidDifferences(humanproteins,elephantproteins,0.005)
    displayTextResults(x,y)
    labels=makeAminoAcidLabels(humanproteins,elephantproteins)
    f1=setupChartData(labels,humanproteins)
    f2=setupChartData(labels,elephantproteins)
    edges=makeEdgeList(labels,y)
    createChart(labels, f1, "Human", f2, "Elephant", edgeList=edges)
    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    # print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    # test.week1Tests()
    # print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    # runWeek1()
    # test.testMakeEdgeList()

    ## Uncomment these for Week 2 ##
    
    # print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    # test.week2Tests()
    # print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    # runWeek2()
    

    ## Uncomment these for Week 3 ##
    
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    
