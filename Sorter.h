#ifndef IR_PROJECTS_SORTER_H
#define IR_PROJECTS_SORTER_H

#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <vector>
#include <algorithm>
#include <filesystem>
#define nodeId_t unsigned int

using namespace std;

class Sorter {
private:
    int blocksCount;
    string runtimeFolder;
    vector<ifstream> blockStreams;
    vector<pair<nodeId_t, nodeId_t>> front;
    vector<bool> frontAvailable;

    void readLineOfStream(int index);

public:
    Sorter(const string &inputFilename, const string &runtimeFolder, int arcPerBlock);

    bool next(pair<nodeId_t, nodeId_t> &value);
};


#endif //IR_PROJECTS_SORTER_H
