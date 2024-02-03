#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <format>
#include <filesystem>
#include "Sorter.h"

#define blocksFolder "C:\\Users\\Admin\\source\\repos\\IR_Projects\\blocks"
#define arcPerBlock 200000

using namespace std;

int main() {

    Sorter s("C:\\Users\\Admin\\source\\repos\\IR_Projects\\data\\web-NotreDame.txt", blocksFolder, arcPerBlock);

    pair<nodeId_t, nodeId_t> tmp;
    vector<pair<nodeId_t, nodeId_t>> aa;
    while(s.next(tmp))
    {
        aa.push_back(tmp);
    }


    return 0;
}
