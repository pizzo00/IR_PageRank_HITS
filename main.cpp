#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <filesystem>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <cmath>
#include "Sorter.h"

#define ull unsigned long long
#define runtimeFolder "/home/pizzolato/IR_Project/runtime"
#define arcPerBlock 200000

using namespace std;

nodeId_t writeMatrixFiles(Sorter &s, bool transposed, string const& colFilename, string const& rowFilename)
{
    std::fstream colFile, rowFile;
    colFile.open(colFilename, std::ios::app | std::ios::binary);
    rowFile.open(rowFilename, std::ios::app | std::ios::binary);

    ull colInd = 0;
    nodeId_t maxId = 0;
    nodeId_t currentRow = 0;
    pair<nodeId_t, nodeId_t> tmp;
    while(transposed ? s.nextTransposed(tmp) : s.next(tmp))
    {
        nodeId_t row = transposed? tmp.second : tmp.first;
        nodeId_t col = transposed? tmp.first : tmp.second;
        if(row != currentRow)
        {
            for(nodeId_t i = currentRow; i < row; i++)
                rowFile.write(reinterpret_cast<char*>(&colInd), sizeof(colInd));
            currentRow = row;
        }
        maxId = max(maxId, col);
        colInd++;
        colFile.write(reinterpret_cast<char*>(&col), sizeof(col));
    }

    maxId = max(maxId, currentRow);
    for(nodeId_t i = currentRow; i < maxId+1; i++)
        rowFile.write(reinterpret_cast<char*>(&colInd), sizeof(colInd));

    rowFile.close();
    return maxId+1;
}

template <typename T>
void writeVector(string const& filename, ull length, T value)
{
    std::fstream f;
    f.open(filename, std::ios::app | std::ios::binary);

    for(ull i = 0; i < length; i++)
        f.write(reinterpret_cast<char*>(&value), sizeof(value));

    f.close();
}

double multiply(
        string const& outFilename,
        nodeId_t nNodes,
        ull const *aRow, nodeId_t const *aCol, double const *aValues,
        double const *b, double const& bDivider = 1.0)
{
    std::fstream f;
    f.open(outFilename, std::ios::app | std::ios::binary);
    double outSum = 0;
    for(nodeId_t r = 0; r < nNodes; r++)
    {
        double res = 0;
        for(ull cIdx = aRow[r]; cIdx < aRow[r+1]; cIdx++)
        {
            nodeId_t col = aCol[cIdx];
            double colValue = aValues == nullptr ? 1 : aValues[cIdx];
            res += colValue * (b[col] / bDivider);
        }
        f.write(reinterpret_cast<char*>(&res), sizeof(res));
        outSum += res;
    }

    f.close();
    return outSum;
}

double error(nodeId_t nNodes, double const *a, double const *newA, double const& aDivider = 1.0, double const& newADivider = 1.0)
{
    double out = 0;
    for(nodeId_t r = 0; r < nNodes; r++)
    {
        out += pow((newA[r] / newADivider) - (a[r] / aDivider), 2);
    }
    return sqrt(out);
}

int main() {
    filesystem::remove_all(runtimeFolder);
    filesystem::create_directory(runtimeFolder);

    Sorter s("/home/pizzolato/IR_Project/data/web-NotreDame.txt", runtimeFolder, arcPerBlock);

    cout << "Begin writing" << endl;
    string aFilename = runtimeFolder + string("/a_0.txt");
    string hFilename = runtimeFolder + string("/h_0.txt");
    string mColFilename = runtimeFolder + string("/mCol.txt");
    string mRowFilename = runtimeFolder + string("/mRow.txt");
    string mtColFilename = runtimeFolder + string("/mtCol.txt");
    string mtRowFilename = runtimeFolder + string("/mtRow.txt");
    cout << "    Matrix" << endl;
    nodeId_t nNodes = writeMatrixFiles(s, false, mColFilename, mRowFilename);
    cout << "    Transposed Matrix" << endl;
    writeMatrixFiles(s, true, mtColFilename, mtRowFilename);
    cout << "    A" << endl;
    writeVector(aFilename, nNodes, (double)1.0);
    cout << "    H" << endl;
    writeVector(hFilename, nNodes, (double)1.0);

    int mRowFile = open(mRowFilename.c_str(), O_RDONLY);
    int mColFile = open(mColFilename.c_str(), O_RDONLY);
    int mtRowFile = open(mtRowFilename.c_str(), O_RDONLY);
    int mtColFile = open(mtColFilename.c_str(), O_RDONLY);
    ull *mRowMap = (ull *) mmap(nullptr, nNodes*sizeof(ull), PROT_READ, MAP_SHARED, mRowFile, 0);
    nodeId_t *mColMap = (nodeId_t *) mmap(nullptr, mRowMap[nNodes-1]*sizeof(nodeId_t), PROT_READ, MAP_SHARED, mColFile, 0);
    ull *mtRowMap = (ull *) mmap(nullptr, nNodes*sizeof(ull), PROT_READ, MAP_SHARED, mtRowFile, 0);
    nodeId_t *mtColMap = (nodeId_t *) mmap(nullptr, mtRowMap[nNodes-1]*sizeof(nodeId_t), PROT_READ, MAP_SHARED, mtColFile, 0);


    double errA = 1000;
    double errH = 1000;
    double aSum = nNodes * 1.0;
    double hSum = nNodes * 1.0;
    double newASum;
    double newHSum;
    ull step = 0;
    int aFile = open(aFilename.c_str(), O_RDONLY);
    int hFile = open(hFilename.c_str(), O_RDONLY);
    double *aMap = (double *) mmap(nullptr, nNodes * sizeof(double), PROT_READ, MAP_SHARED, aFile, 0);
    double *hMap = (double *) mmap(nullptr, nNodes * sizeof(double), PROT_READ, MAP_SHARED, hFile, 0);
    while (step < 50 && (errA > 1e-10 || errH > 1e-10))
    {
        step++;
        string newAFilename = runtimeFolder + string("/a_") + to_string(step) + string(".txt");
        string newHFilename = runtimeFolder + string("/h_") + to_string(step) + string(".txt");

        cout << "Begin computation" << endl;
        newASum = multiply(newAFilename, nNodes, mtRowMap, mtColMap, nullptr, hMap, hSum);
        newHSum = multiply(newHFilename, nNodes, mRowMap, mColMap, nullptr, aMap, aSum);

        int newAFile = open(newAFilename.c_str(), O_RDONLY);
        double *newAMap = (double *) mmap(nullptr, nNodes * sizeof(double), PROT_READ, MAP_SHARED, newAFile, 0);
        int newHFile = open(newHFilename.c_str(), O_RDONLY);
        double *newHMap = (double *) mmap(nullptr, nNodes * sizeof(double), PROT_READ, MAP_SHARED, newHFile, 0);
        errA = error(nNodes, aMap, newAMap, aSum, newASum);
        errH = error(nNodes, hMap, newHMap, hSum, newHSum);

        //Close old file
        close(aFile);
        close(hFile);

        //Swap files
        aFile = newAFile;
        hFile = newHFile;
        aMap = newAMap;
        hMap = newHMap;
        aSum = newASum;
        hSum = newHSum;

        cout << "error A: " << errA << endl;
        cout << "error H: " << errH << endl;

        close(newAFile);
    }

    close(aFile);
    close(hFile);
    close(mRowFile);
    close(mColFile);
    close(mtRowFile);
    close(mtColFile);



//    if (atRowMap == MAP_FAILED || atColMap == MAP_FAILED)
//    {
//        perror("Error mmapping the files");
//        exit(EXIT_FAILURE);
//    }

    for(int i = 0; i < 10; i++)
    {
//        cout << atColMap[i] << endl;
    }

    return 0;
}