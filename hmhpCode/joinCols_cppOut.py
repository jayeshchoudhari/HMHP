import numpy as np
from collections import defaultdict

origFile = open("../wuvPerGroup_sorted_wuvVals_log_1.5.txt")

origWuv = defaultdict(float)

for line in origFile:
        flds = line.strip().split()
        gid = flds[0].strip()
        wuvVal = float(flds[1].strip())
        origWuv[gid] = wuvVal

origFile.close()


# predFile = open("B5_GroupedUserUserInf_OurModel.txt")
# predFile = open("wuvEstFile_fromModeParent.txt")
predFile = open("wuvEstFile_fromLastIterationParent.txt")
# predFile = open("B5_GroupedUserUserInf_OurModel_levels_restricted.txt")
# predFile = open("B5_GroupedUserUserInf_OurModel_1L_Events.txt")
# predFile = open("B5_GroupedUserUserInf_OurModel_All_Events.txt")
# predFile = open("B5_GroupedUserUserInf_OurModel_All_Events_unb.txt")
# predFile = open("B5_GroupedUserUserInf_OurModel_10M_Events.txt")
# predFile = open("B5_GroupedUserUserInf_OurModel_1Million_Events.txt")
# predFile = open("B5_GroupedUserUserInf_OurModel_1Lakh_Events.txt")
# predFile = open("B5_GroupedUserUserInf_OurModel_10Thou_Events.txt")

predWuv = defaultdict(float)
predNuv = defaultdict(int)
predN = defaultdict(float)
predNEdge = defaultdict(float)
predNuvByNu = defaultdict(float)
triggeredEdges = defaultdict(int)
untriggeredEdges = defaultdict(int)
edgesInGroup = defaultdict(int)
numTimesNuAdded = defaultdict(int)
zerosInGroup = defaultdict(int)
avgPoisson = defaultdict(float)

noOfNonZeroNus = defaultdict(float)
meanNonZeroNus = defaultdict(float) 
medianNonZeroNus = defaultdict(float)
perc75Val = defaultdict(float)
maxNonZeroNus = defaultdict(float)
aboveMedWuv = defaultdict(float)

for line in predFile:
        flds = line.strip().split()
        gid = flds[0].strip()
        wuvVal = float(flds[1].strip())
        predWuv[gid] = wuvVal
        predNuv[gid] = int(flds[2].strip())
        predN[gid] = float(flds[3].strip())
        '''
        predNEdge[gid] = float(flds[4].strip())
        predNuvByNu[gid] = float(flds[5].strip())
        triggeredEdges[gid] = int(flds[6].strip())
        untriggeredEdges[gid] = int(flds[7].strip())
        edgesInGroup[gid] = int(flds[8].strip())
        numTimesNuAdded[gid] = int(flds[9].strip())
        zerosInGroup[gid] = int(flds[10].strip())
        avgPoisson[gid] = float(flds[11].strip())
        

        noOfNonZeroNus[gid] = float(flds[4].strip())
        meanNonZeroNus[gid] = float(flds[5].strip())
        medianNonZeroNus[gid] = float(flds[6].strip())
        perc75Val[gid] = float(flds[7].strip())
        maxNonZeroNus[gid] = float(flds[8].strip())
        aboveMedWuv[gid] = float(flds[9].strip())
        '''

predFile.close()

# joinedFile = open("origPredWuvEstCpp.txt", "w")
joinedFile = open("origPredWuvEstLastIterationParent.txt", "w")
# joinedFile = open("origPredWuvPerGroup_10Thou_Events.txt", "w")
# joinedFile = open("origPredWuvPerGroup_1Lakh_Events.txt", "w")
# joinedFile = open("origPredWuvPerGroup_1Million_Events.txt", "w")
# joinedFile = open("origPredWuvPerGroup_10M_Events.txt", "w")
# joinedFile = open("origPredWuvPerGroup_10M_Events_SumPoisson.txt", "w")
# joinedFile = open("origPredWuvPerGroup_AllEvents_1M_unb.txt", "w")
# joinedFile = open("origPredWuvPerGroup_AllEvents_1M.txt", "w")
# joinedFile = open("origPredWuvPerGroup_1L_Events.txt", "w")
# joinedFile = open("origPredWuvPerGroup_level_restricted.txt", "w")

joinedFile.write("gid, origWuv, PredWuv, SumNuv, SumNuOverAllEdges, ape\n")

apeList = []
apeListCount = 0

above100SumNuTotalApe = []
above100SumNuCount = 0

for gid in predWuv:
        # untriggered = edgesInGroup[gid] - triggeredEdges[gid]
        if origWuv[gid] > 0:
            ape = abs((predWuv[gid]) - origWuv[gid])*1.0/origWuv[gid]
            joinedFile.write(gid + " " + str(origWuv[gid]) + " " + str(predWuv[gid]) +  " " + str(predNuv[gid]) + " " + str(predN[gid]) + " " + str(ape) + "\n")

            apeList.append(ape)
            apeListCount += 1


            if predN[gid] > 100:
                above100SumNuTotalApe.append(ape)
                above100SumNuCount += 1
            
        # if untriggered > 0: 
        '''
        if triggeredEdges[gid] > 0: 
             ratio = predNuv[gid]*1.0/triggeredEdges[gid]
        else:
             ratio = 0
        '''

        # joinedFile.write(gid + " " + str(origWuv[gid]) + " " + str(predWuv[gid]) + " " + str(predNuv[gid]) + " " + str(predN[gid]) + " " + str(predNEdge[gid]) + " " + str(predNuvByNu[gid]) + " " + str(triggeredEdges[gid]) + " " + str(untriggeredEdges[gid]) + " " + str(edgesInGroup[gid]) + " " + str(ape) + " " + str(numTimesNuAdded[gid]) + " " + str(zerosInGroup[gid]) + " " + str(numTimesNuAdded[gid] - zerosInGroup[gid]) + " " + str(avgPoisson[gid]) + "\n")
        
        # joinedFile.write(gid + " OrigWuv = " + str(origWuv[gid]) + " PredWuv = " + str(predWuv[gid]) + " NuvSumNume = " + str(predNuv[gid]) + " NusumDenom = " + str(predN[gid]) + " ape = " + str(ape) + " length = " + str(noOfNonZeroNus[gid]) + " mean = " + str(meanNonZeroNus[gid]) + " median = " + str(medianNonZeroNus[gid]) + " perc75val = " + str(perc75Val[gid])  + " max = " + str(maxNonZeroNus[gid]) + " aboveMedWuv = " + str(aboveMedWuv[gid]) +  "\n")
        # joinedFile.write(gid + " " + str(origWuv[gid]) + " " + str(predWuv[gid]) +  " " + str(predNuv[gid]) + " " + str(predN[gid]) + " " + str(ape) + "\n")

        # jodiff = origWuv[gid] - predWuv[gid]
        # if jodiff < 0:
            # print(gid, jodiff)
	
		
print "Avg Ape = ", np.sum(apeList)*1.0/apeListCount, "\nMedian Ape = ", np.median(apeList), "\ngroups = ", apeListCount, "\nAbove 100 -- Avg Ape = ", np.sum(above100SumNuTotalApe)*1.0/above100SumNuCount, "\nabove 100 median ape = ", np.median(above100SumNuTotalApe), "\ngroups = ", above100SumNuCount,

joinedFile.close()
    



