
# Readme for synthetic experiments:

The code for generating the (semi-)synthetic data is in file named "generateSemiSyntheticData_usingDists.py". 

The input files required for the code are:
    - User-User Influence file
    - user-topic distributions
	- topic-topic distributions
	- topic-word distributions
	- user base rates

The code outputs events file and documents file. Each line of events file describes an event with 5 space separated values as follows:

> event_time user_node parent_event_id topic_id level_info

The number of lines in the document file is equal to the number of lines in the event file. 
Each line describes a document as space separated word-ids, and this document line corresponds to the event on the same line in the events file.

The data generation code can be executed as follows:

> python generateSemiSyntheticData_usingDists.py 1000 1 1231241 1000 5 500  

The various arguments are **argv[1]** = *Tmax (time horizon)*, **argv[2]** = *1 (some constant which is not used for now)*, **argv[3]** = *random_seed*, **argv[4]** = *max number of events to be generated*, **argv[5]** = *no. of topics*, **argv[6]** = *vocabsize*


## File format for input files:

1. *User-User Influence file*: Each line specifies the followers (**v**) of each user (**u**) along with the influence of user (**u**) on each of its followers (**v**)
	
	> no_Of_followers user_u_id user_v1_id W(u,v1) user_v2_id W(u,v2) ... user_vn_id W(u,vn)

2. *User-topic Distribution*: Each line give the user-id followed by No_Of_Topics many probability values, where i-th value indicates the users preference for i-th topic
	
	> user_id pr(t1) pr(t2) ... pr(tk)

3. *topic-topic distribution*: Each line gives the topic-id followed by No_Of_Topics many probability values, where i-th value indicates the correlation(probability) of topic-id with i-th topic
	
	> topic_id pr(topic_id,t1) pr(topic_id,t2) ... pr(topic_id,tk)  

4. *topic-word distribution*: Each line gives the topic-id followed by *Vocabsize* many probability values, where i-th value indicates the probability of i-th word under topic-id
	
	> topic_id pr(w1) pr(w2) ... pr(wm)

5. *User Base Rate*: Each line is 2 space separated values. First value is user-id and second value is base rate of user-id
	> user_id base_rate