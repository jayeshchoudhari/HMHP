# HMHP
Hidden Markov Hawkes Process is probabilistic model proposed in *Discovering Topical Interactions in Text-based Cascades using Hidden Markov Hawkes Processes*. 

The folder named 

1. *syntDataGeneration* contains the code for generating the synthetic data as per the proposed model
2. *hmhpCode* contains the inference code for the HMHP model 
3. *nhmCode* contains the inference code for the model that performs *network reconstruction* and *parent identification* task simultaneously. This is similar to Network Hawkes Model
4. *hawkesLDAMM* is the model that does the *topic identification* using LDAMM and *parent identification* depending on the time kernel. The *network reconstruction* task is done once the *parent identification* is in place
5. *hawkesDIAG* is the model that is similar to that of *HMHP*. The only difference is that *hawkesDIAG* considers identity *topic-topic interactions*
6. *inputFiles* contains the files that are required to run each of the models