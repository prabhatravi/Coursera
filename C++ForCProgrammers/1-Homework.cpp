//  Convert this program to C++
//  change to C++ io
//  change to one line comments
//  change defines of constants to const
//  change array to vector<>
//  inline any short function

#include <iostream>
#include <vector>
using namespace std;

const int N=40;

// Initialize vector data
inline void init(int& accum, vector<int>& data)
{
  for(int i = 0; i < data.size(); ++i)
    data.at(i) = i;	// Initialize each position with its index value
}

// Sum all elements of vector data
inline void sum(int& accum, vector<int>& data)
{
  for(int i = 0; i < data.size(); ++i)
    accum += data.at(i);	// Accumulate the sum of each position
}

int main()
{
  int accum = 0;	// Create accumulator
  vector<int> data(N);	// Create vector data
  
  init(accum, data);	// Initialize vector data
  sum(accum, data);	// Sum all elements of vector data
  cout << "sum is " << accum << endl;	// Show the sum of all elements of vector data
  
  return 0;
}
