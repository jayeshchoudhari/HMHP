#! /bin/bash

# Evaluation for HMHP results...

echo -e 'Get the mode parent and topic assignments -- '

g++ -std=c++11 -g -o modeAssign getModeValueFromAssignments.cpp

./modeAssign ./B5_topicAssignments_OurModel.txt ./B5_topicAssignments_ModeAssignment.txt
./modeAssign ./B5_parentAssignments_OurModel.txt ./B5_parentAssignments_ModeAssignment.txt

echo -e 'Got the mode assignments'


echo -e 'Evaluating Topic Assignments...'

g++ -std=c++11 -g -o evalTopics docDocTopicAssignmentEvaluation.cpp
./evalTopics ./B5_topicAssignments_ModeAssignment_estAll.txt


echo -e 'Evaluting Parent Assignments...'

g++ -std=c++11 -g -o accpred accPredParents.cpp
./accpred ./B5_parentAssignments_ModeAssingment_estAll.txt


echo -e 'Getting Recall values for Parent Assignments...'

g++ -std=c++11 -g -o precAtK precAtKParentAssignment.cpp
./precAtK ./B5_avgParProb_OurModel_estAll.txt ./B5_parentAssignments_RecallAtK_estAll.txt


echo -e 'Recall Values @ 1, 3, 5, 7, 10, 50, 100'

cut -d' ' -f2 ./B5_parentAssignments_RecallAtK_estAll.txt | awk '{t += $1} END {print t/NR}' &&  cut -d' ' -f3 ./B5_parentAssignments_RecallAtK_estAll.txt | awk '{t += $1} END {print t/NR}' && cut -d' ' -f4 ./B5_parentAssignments_RecallAtK_estAll.txt | awk '{t += $1} END {print t/NR}' &&  cut -d' ' -f5 ./B5_parentAssignments_RecallAtK_estAll.txt | awk '{t += $1} END {print t/NR}' && cut -d' ' -f6 ./B5_parentAssignments_RecallAtK_estAll.txt | awk '{t += $1} END {print t/NR}' && cut -d' ' -f7 ./B5_parentAssignments_RecallAtK_estAll.txt | awk '{t += $1} END {print t/NR}'


echo -e 'Wuv evaluation....'
python joinCols_cppOut.py