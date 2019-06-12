#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <string>
#include <cstddef>
#include <vector>
#include <map>
#include <algorithm>
#include <set>
#include <iostream>
#include <assert.h>
#include "../include/sirgl32.h"
#include "../include/sirg.h"
#define BUFFERSIZE 1024
using namespace std;

struct CompareByDegree{
  CompareByDegree (int* degree) {
    this->degree = degree;
  }
  bool operator() (int i, int j){
    if (degree[i] > degree[j])
      return false;
    else if (degree[i] < degree[j])
      return true;
    
    return i < j;
  }
  int* degree;
};

void read_input_degree(const char* fileName, int* numRows, int* numCols, int** nonZeroIndices, int** numColIndices, int* numElems, int** nextVertexIndices, int* maxDeg, int** degree)
{
	FILE *fin = fopen(fileName, "r");
	if (fin == NULL)
	{
		printf("File %s not found\n", fileName);
		return ;
	}
	char* buffer = NULL;
	size_t len = BUFFERSIZE;
	int m,n,numEdges;
	while (getline(&buffer, &len, fin) != 0)
	{
		if (buffer[0] == '%')
			continue;
		sscanf(buffer, "%d %d %d\n", &m, &n, &numEdges);
		break;
	}
	*numRows = m;
	*numCols = n;
	assert(*numRows == *numCols);

	vector<set<int> > hashmap(*numRows);
	for (int i = 0; i < numEdges; ++i)
	{
		int u,v;
		getline(&buffer, &len, fin);
		sscanf(buffer, "%d %d\n", &u, &v);
    if (u == v)
      continue;
		u = u - 1;
		v = v - 1;
    hashmap[u].insert(v);
    hashmap[v].insert(u);
	}

  *numColIndices = (int *)malloc(sizeof(int)*(*numRows+1));
	memset(*numColIndices, 0, *numRows + 1);
  *degree = (int *)malloc(sizeof(int)*(*numRows));
  memset(*degree, 0, *numRows);
	*numElems = 0;
	for (int it = 0; it < *numRows; it++)
	{
		(*numColIndices)[it + 1] = hashmap[it].size();
    (*degree)[it] = hashmap[it].size();
		*numElems += hashmap[it].size();
	}

	for (int i = 1; i < *numRows + 1; ++i)
	{
		(*numColIndices)[i] += (*numColIndices)[i-1];
	}
  
  *nextVertexIndices = (int *)malloc(sizeof(int)*(*numRows));
	*nonZeroIndices = (int *)malloc(sizeof(int)*(*numElems));
	int index = 0;
  *maxDeg = 0;
	for (int i = 0; i < *numRows; ++i)
	{
      if (hashmap[i].size() > *maxDeg)
        *maxDeg = hashmap[i].size();

      vector<int> v(hashmap[i].begin(), hashmap[i].end());
      sort(v.begin(), v.end(), CompareByDegree(*degree));

      int higherIndex = -1;
			for (vector<int>::iterator it = v.begin(); it != v.end(); it++)
			{
				(*nonZeroIndices)[index] = *it;
        if (higherIndex == -1 && (*degree)[*it] > (*degree)[i])
          higherIndex = index;
        else if (higherIndex == -1 && (*degree)[*it] == (*degree)[i]){
          if (*it > i)
            higherIndex = index;
        }
				index++;
			}
      if (higherIndex == -1)
        higherIndex = index;
      (*nextVertexIndices)[i] = higherIndex;
	}
}

void read_input_vertex(const char* fileName, int* numRows, int* numCols, int** nonZeroIndices, int** numColIndices, int* numElems, int** nextVertexIndices, int* maxDeg)
{
	FILE *fin = fopen(fileName, "r");
	if (fin == NULL)
	{
		printf("File %s not found\n", fileName);
		return ;
	}
	char* buffer = NULL;
	size_t len = BUFFERSIZE;
	int m,n,numEdges;
	while (getline(&buffer, &len, fin) != 0)
	{
		if (buffer[0] == '%')
			continue;
		sscanf(buffer, "%d %d %d\n", &m, &n, &numEdges);
		break;
	}
	*numRows = m;
	*numCols = n;
	assert(*numRows == *numCols);

	vector<set<int> > hashmap(*numRows);
	for (int i = 0; i < numEdges; ++i)
	{
		int u,v;
		getline(&buffer, &len, fin);
		sscanf(buffer, "%d %d\n", &u, &v);
    if (u == v)
      continue;
		u = u - 1;
		v = v - 1;
    hashmap[u].insert(v);
    hashmap[v].insert(u);
	}
  *numColIndices = (int *)malloc(sizeof(int)*(*numRows+1));
	memset(*numColIndices, 0, *numRows + 1);
	*numElems = 0;
	for (int it = 0; it < *numRows; it++)
	{
		(*numColIndices)[it + 1] = hashmap[it].size();
		*numElems += hashmap[it].size();
	}

	for (int i = 1; i < *numRows + 1; ++i)
	{
		(*numColIndices)[i] += (*numColIndices)[i-1];
	}
  
  *nextVertexIndices = (int *)malloc(sizeof(int)*(*numRows));
	*nonZeroIndices = (int *)malloc(sizeof(int)*(*numElems));
	int index = 0;
  *maxDeg = 0;
	for (int i = 0; i < *numRows; ++i)
	{
      if (hashmap[i].size() > *maxDeg)
        *maxDeg = hashmap[i].size();
      int higherIndex = -1;
			for (set<int>::iterator it = hashmap[i].begin(); it != hashmap[i].end(); it++)
			{
				(*nonZeroIndices)[index] = *it;
        if (higherIndex == -1 && *it > i)
          higherIndex = index;
				index++;
			}
      if (higherIndex == -1)
        higherIndex = index;
      (*nextVertexIndices)[i] = higherIndex;
	}
}


bool verifyGraphColoring(int* colors, int numRows, int* nonZeroIndices, int* numColIndices){
  for (int i = 0; i < numRows; i++){
    for(int j = numColIndices[i]; j < numColIndices[i+1]; j++){
      if (colors[i] == colors[nonZeroIndices[j]]){
        return false;
      }
    }
  }
  return true;
}

int countColors(int* colors, int numRows){
  set<int> hash;
  for (int i = 0; i < numRows; i++){
    hash.insert(colors[i]);
  }
  return hash.size();
}

void checkColorHoles(int* colors, int numRows, int* nonZeroIndices, int* numColIndices){
  int adjColors[300];
  for (int j = 0; j < 300; j++)
      adjColors[j] = 0;

  for (int i = 0; i < numRows; i++){
    adjColors[colors[i]] = 1;
  }

  for (int i = 0; i < 300; i++){
    if (adjColors[i] == 0)
        cout << i << endl;
  }

}

int main(int argc, char const *argv[])
{
	if (argc != 3)
  {
		printf("Use: ./a.out <fileName> <0/1> (0 for degree, 1 for vertex)\n");
    return 0;
  }

	int numRows, numCols, numElems, maxDeg;
	int *nonZeroIndices, *numColIndices, *nextVertexIndices, *degree;

	if (strstr(argv[1], ".mtx") != NULL)
  {
    if (atoi(argv[2]) == 0)
  		read_input_degree(argv[1], &numRows, &numCols, &nonZeroIndices, &numColIndices, &numElems, &nextVertexIndices, &maxDeg, &degree);
    else
      read_input_vertex(argv[1], &numRows, &numCols, &nonZeroIndices, &numColIndices, &numElems, &nextVertexIndices, &maxDeg);
  }
  else
    return 0;

  string filename(argv[1]);
  size_t found = filename.find_last_of("/\\");
  cout << filename.substr(found+1, filename.length()-4) << " ";
/*  
  cout << "Reading Graph " << argv[1] << endl;
  cout << "Number of Vertices: " << numRows << endl;
  cout << "Number of Edges: " << numElems << endl;
  cout << "Average Degree: " << (numElems * 1.0)/(numRows * 1.0) << endl;
  cout << "Maximum Degree: " << (maxDeg) << endl;
*/
  int *colors;
  if (maxDeg < 32){
    colors = colorGraphSIRGL32(numRows, nonZeroIndices, numColIndices, numElems, nextVertexIndices, maxDeg);
  }else{
    colors = colorGraphSIRG(numRows, nonZeroIndices, numColIndices, numElems, nextVertexIndices, maxDeg);
  }
/*
  if (verifyGraphColoring(colors, numRows, nonZeroIndices, numColIndices)){
    cout << "Graph Colored Correctly" << endl;
  }
  else{
    cout << "Incorrect Graph Coloring" << endl;
  }
*/  
	return 0;
}
