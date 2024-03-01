#ifndef IR_PROJECTS_SORTER_H
#define IR_PROJECTS_SORTER_H

#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <vector>
#include <algorithm>
#define nodeId_t unsigned int

using namespace std;

class Sorter {
private:
    int blocksCount;
    string runtimeFolder;
    vector<ifstream> blockStreams;
    vector<ifstream> blockTStreams;
    vector<pair<nodeId_t, nodeId_t> > front;
    vector<pair<nodeId_t, nodeId_t> > frontT;
    vector<bool> frontAvailable;
    vector<bool> frontTAvailable;

    void readLineOfStream(int index);
    void readLineOfStreamT(int index);

public:
    Sorter(const string &inputFilename, const string &runtimeFolder, int arcPerBlock);

    bool next(pair<nodeId_t, nodeId_t> &value);
    bool nextTransposed(pair<nodeId_t, nodeId_t> &value);

};


#endif //IR_PROJECTS_SORTER_H
