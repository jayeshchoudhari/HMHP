#ifndef NAMESPACE_H
#define NAMESPACE_H

#include <iostream>				//for basic C++ functions 
#include <cstdio>				//for std C functions
#include <fstream>				//for i/o stream
#include <sstream>				//for parse data using stringstreams...
#include <string>				//for C++ string functions
#include <cstring>				//we need this for memset... and string functions from C
#include <unordered_map>		//to maintain the map of user-tweet_count, and user_user_tweet_count
#include <map>		            //to maintain the map of user-tweet_count, and user_user_tweet_count
#include <ctime>				//for time realted functions... to convert user-formatted time to timestamp
#include <vector>
#include <assert.h>				//for assertions
#include <limits>		
#include <cmath>
#include <algorithm>
#include <numeric>	
#include <random>
#include <chrono>
#include <cfenv>
#include <iomanip>


using ll = long long int;
using li = long int;
using ui = unsigned int;

#define INF std::numeric_limits<int>::max()					//use limits for infinity
#define PI 3.14159265

#endif