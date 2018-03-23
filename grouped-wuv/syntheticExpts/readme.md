For synthetic experiments:

The code for generating the (semi-)synthetic data is in file named "generateSemiSyntheticData_usingDists.py". 
The input files required for the code are:
	-- User-User Influence file
	-- user-topic distributions
	-- topic-topic distributions
	-- topic-word distributions
	-- user base rates

The code outputs events file and documents file. Each line of events file describes an event with 5 space separated values as follows:

event_time user_node parent_event_id topic_id level_info

The number of lines in the document file is equal to the number of lines in the event file. 
Each line describes a document as space separated word-ids, and this document line corresponds to the event on the same line in the events file.

The data generation code can be executed as follows:

python generateSemiSyntheticData_usingDists.py 302400 1 12397232

The various arguments are argv[1] = Tmax (time horizon), argv[2] = 1 (constant alpha=1 for time-kernel), argv[3] = random_seed