#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <filesystem>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include "Sorter.h"

#define ull unsigned long long
#define runtimeFolder "/mnt/c/Users/Admin/source/repos/IR_Projects/runtime"
#define arcPerBlock 200000

using namespace std;

ull writeAtFile(Sorter &s, string const& colFilename, string const& rowFilename)
{
    std::fstream atFileCol, atFileRow;
    atFileCol.open(colFilename, std::ios::app | std::ios::binary);
    atFileRow.open(rowFilename, std::ios::app | std::ios::binary);

    ull colInd = 0;
    nodeId_t currentRow = 0;
    pair<nodeId_t, nodeId_t> tmp;
    while(s.next(tmp))
    {
        if(tmp.second != currentRow)
        {
            for(nodeId_t i = currentRow; i < tmp.second; i++)
                atFileRow.write(reinterpret_cast<char*>(&colInd), sizeof(colInd));
            currentRow = tmp.second;
        }
        colInd++;
        atFileCol.write(reinterpret_cast<char*>(&tmp.first), sizeof(tmp.first));
    }
    atFileRow.close();
    return currentRow;
}

int main() {
    filesystem::remove_all(runtimeFolder);
    filesystem::create_directory(runtimeFolder);

    Sorter s("/mnt/c/Users/Admin/source/repos/IR_Projects/data/super_small.txt", runtimeFolder, arcPerBlock);

    string atColFilename = runtimeFolder + string("/atCol.txt");
    string atRowFilename = runtimeFolder + string("/atRow.txt");
    ull nColumns = writeAtFile(s, atColFilename, atRowFilename);

    int atRowFile = open(atRowFilename.c_str(), O_RDONLY);
    ull *atRowMap = (ull *) mmap(nullptr, nColumns*sizeof(ull), PROT_READ, MAP_SHARED, atRowFile, 0);

    if (atRowMap == MAP_FAILED)
    {
        perror("Error mmapping the file");
        exit(EXIT_FAILURE);
    }

    for(int i = 0; i < nColumns; i++)
    {
        cout << atRowMap[i] << endl;
    }

    return 0;
}