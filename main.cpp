#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <queue>
#include <chrono>
#include <set>
#include <filesystem>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <cmath>
#include <thread>
#include <future>
#include "Sorter.h"

#define ull unsigned long long
#define runtimeFolder "./runtime"
#define arcPerBlock 200000
#define WORKERS 10

using namespace std;

nodeId_t writeMatrixFiles(
        Sorter &s,
        bool transposed,
        string const& colFilename,
        string const& rowFilename,
        string const& valueFilename,
        string const& danglingFilename
        )
{
    double one = 1;
    double zero = 0;

    std::fstream colFile, rowFile, valueFile, dangFile;
    colFile.open(colFilename, std::ios::app | std::ios::binary);
    rowFile.open(rowFilename, std::ios::app | std::ios::binary);
    valueFile.open(valueFilename, std::ios::app | std::ios::binary);
    dangFile.open(danglingFilename, std::ios::app | std::ios::binary);

    ull startColInd = 0;
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
            for(int i = 0; i < (colInd - startColInd); i++)
            {
                double val = 1.0 / (double)(colInd - startColInd);
                valueFile.write(reinterpret_cast<char *>(&val), sizeof(val));
            }
            startColInd = colInd;

            if(!transposed)
                dangFile.write(reinterpret_cast<char*>(&one), sizeof(one));
            rowFile.write(reinterpret_cast<char*>(&colInd), sizeof(colInd));
            for(nodeId_t i = currentRow + 1; i < row; i++)
            {
                if(!transposed)
                    dangFile.write(reinterpret_cast<char*>(&zero), sizeof(zero));
                rowFile.write(reinterpret_cast<char *>(&colInd), sizeof(colInd));
            }
            currentRow = row;
        }

        maxId = max(maxId, col);
        colInd++;
        colFile.write(reinterpret_cast<char*>(&col), sizeof(col));
    }

    for(int i = 0; i < (colInd - startColInd); i++)
    {
        double val = 1.0 / (double)(colInd - startColInd);
        valueFile.write(reinterpret_cast<char *>(&val), sizeof(val));
    }

    maxId = max(maxId, currentRow);
    if(!transposed)
        dangFile.write(reinterpret_cast<char*>(&one), sizeof(one));
    rowFile.write(reinterpret_cast<char*>(&colInd), sizeof(colInd));
    for(nodeId_t i = currentRow + 1; i < maxId+1; i++)
    {
        if(!transposed)
            dangFile.write(reinterpret_cast<char*>(&zero), sizeof(zero));//dangling
        rowFile.write(reinterpret_cast<char*>(&colInd), sizeof(colInd));
    }

    colFile.close();
    rowFile.close();
    valueFile.close();
    dangFile.close();
    return maxId+1; // Row count
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
        double* outFile,
        nodeId_t start,
        nodeId_t end,
        ull const *aRow, nodeId_t const *aCol, double const *aValues,
        double const *b, double const& bMultiplier = 1.0, double const& outputSum = 0)
{
    double outSum = 0;
    for(nodeId_t r = start; r < end; r++)
    {
        double res = 0;
        for(ull cIdx = aRow[r]; cIdx < aRow[r+1]; cIdx++)
        {
            nodeId_t col = aCol[cIdx];
            double colValue = aValues == nullptr ? 1 : aValues[cIdx];
            res += colValue * b[col] * bMultiplier;
        }
        res += outputSum;
        outFile[r] = res;
        outSum += res;
    }

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

double sumIf(nodeId_t nNodes, double const *valArr, bool const *ifArr, double const& divider = 1.0)
{
    double out = 0;
    for(nodeId_t r = 0; r < nNodes; r++)
    {
        if(ifArr[r])
            out += valArr[r] / divider;
    }
    return sqrt(out);
}

template <typename T, typename T2>
set<T> getSetOfFirstK(int k, vector<pair<T, T2>> const& in)
{
    set<T> out;
    for(int i = 0; i < k; i++)
        out.insert(in[i].first);
    return out;
}

template <typename T>
double jaccardIndex(set<T> const& s1, set<T> const& s2)
{
    double s1Size = s1.size();
    double s2Size = s2.size();

    // Get the intersection set
    set<T> intersection;
    set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), inserter(intersection, intersection.begin()));

    double intersectionSize = intersection.size();

    return intersectionSize / (s1Size + s2Size - intersectionSize);
}

vector<pair<nodeId_t, double>> getVectorTopK(nodeId_t nNodes, int k, double const *a, double const& aDivider = 1.0)
{
    priority_queue<pair<double, nodeId_t>, vector<pair<double, nodeId_t>>, greater<>> pq; //Min-heap
    for (nodeId_t r = 0; r < nNodes; r++) {
        pq.emplace(a[r], r);  //add to min-heap

        if (pq.size() > k)
            pq.pop();   //remove the top element (smallest) once the queue reaches the size K
    }

    vector<pair<nodeId_t, double>> out;
    for (int i = 0; i < k; i++) {
        out.emplace_back(pq.top().second, pq.top().first / aDivider);
        pq.pop();
    }
    sort(out.begin(), out.end(), [](pair<nodeId_t, double> const& l, pair<nodeId_t, double> const& r) {
        return make_pair(l.second, l.first) > make_pair(r.second, r.first);
    });

    return out;
}

vector<pair<nodeId_t, nodeId_t>> inDegree(int k, nodeId_t nNodes, ull *mtRowMap)
{
//    cout << "In Degree" << endl;
    //<qty, id>
    priority_queue<pair<nodeId_t, nodeId_t>, vector<pair<nodeId_t, nodeId_t>>, greater<>> pq; //Min-heap
    pq.emplace(mtRowMap[0], 0);  //add to min-heap
    for (nodeId_t r = 1; r < nNodes; r++)
    {
        pq.emplace(mtRowMap[r] - mtRowMap[r-1], r);  //add to min-heap

        if (pq.size() > k)
            pq.pop();   //remove the top element (smallest) once the queue reaches the size K
    }

    //<id, qty>
    vector<pair<nodeId_t, nodeId_t>> out;
    for (int i = 0; i < k; i++) {
        out.emplace_back(pq.top().second, pq.top().first);
        pq.pop();
    }
    sort(out.begin(), out.end(), [](pair<nodeId_t, double> const& l, pair<nodeId_t, double> const& r) {
        return make_pair(l.second, l.first) > make_pair(r.second, r.first);
    });

    return out;
}

pair<vector<pair<nodeId_t, double>>, vector<pair<nodeId_t, double>>> hits(int k, nodeId_t nNodes, ull *mRowMap, nodeId_t *mColMap, ull *mtRowMap, nodeId_t *mtColMap)
{
    string aFilename = runtimeFolder + string("/a_0.txt");
    string hFilename = runtimeFolder + string("/h_0.txt");

    writeVector(aFilename, nNodes, (double)1.0);
    writeVector(hFilename, nNodes, (double)1.0);

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

        int newAFile = open(newAFilename.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0x0777);
        int newHFile = open(newHFilename.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0x0777);
        ftruncate(newAFile, nNodes * sizeof(double));
        ftruncate(newHFile, nNodes * sizeof(double));
        double *newAMap = (double *)mmap (nullptr, nNodes * sizeof(double), PROT_READ | PROT_WRITE, MAP_SHARED, newAFile, 0);
        double *newHMap = (double *)mmap (nullptr, nNodes * sizeof(double), PROT_READ | PROT_WRITE, MAP_SHARED, newHFile, 0);


//        cout << "Begin computation " << step << endl;
        std::future<double> futures[WORKERS];
        nodeId_t nodePerWorkers = nNodes / WORKERS;

        for(int w = 0; w < WORKERS; w++) {
            nodeId_t start = w * nodePerWorkers;
            nodeId_t end = w == WORKERS-1 ? nNodes : (w+1) * nodePerWorkers;
            futures[w] = std::async(launch::async, multiply, newAMap, start, end, mtRowMap, mtColMap, nullptr, hMap,
                                    1.0 / hSum, 0);
        }
        newASum = 0;
        for(int w = 0; w < WORKERS; w++)
            newASum += futures[w].get();

        for(int w = 0; w < WORKERS; w++) {
            nodeId_t start = w * nodePerWorkers;
            nodeId_t end = w == WORKERS-1 ? nNodes : (w+1) * nodePerWorkers;
            futures[w] = std::async(launch::async, multiply, newHMap, start, end, mRowMap, mColMap, nullptr, aMap,
                                    1.0 / aSum, 0);
        }
        newHSum = 0;
        for(int w = 0; w < WORKERS; w++)
            newHSum += futures[w].get();

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

//        cout << "error A: " << errA << endl;
//        cout << "error H: " << errH << endl;

        close(newAFile);
    }

    auto topA = getVectorTopK(nNodes, k, aMap, aSum);
    auto topH = getVectorTopK(nNodes, k, hMap, hSum);

    close(aFile);
    close(hFile);

    return {topA, topH};
}

vector<pair<nodeId_t, double>> pagerank(int k, double d, nodeId_t nNodes, ull *mtRowMap, double *mtValueMap, nodeId_t *mtColMap, bool *mDanglingMap)
{
    double telep = (1-d)/nNodes;
    string pFilename = runtimeFolder + string("/p_0.txt");

    writeVector(pFilename, nNodes, (double)1.0/nNodes);

    double err = 1000;
    double pSum = 1.0;
    double newPSum;
    ull step = 0;
    int pFile = open(pFilename.c_str(), O_RDONLY);
    double *pMap = (double *) mmap(nullptr, nNodes * sizeof(double), PROT_READ, MAP_SHARED, pFile, 0);
    while (step < 50 && err > 1e-10)
    {
        step++;
        string newPFilename = runtimeFolder + string("/p_") + to_string(step) + string(".txt");

        int newPFile = open(newPFilename.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0x0777);
        ftruncate(newPFile, nNodes * sizeof(double));
        double *newPMap = (double *)mmap (nullptr, nNodes * sizeof(double), PROT_READ | PROT_WRITE, MAP_SHARED, newPFile, 0);

        double dangSum = sumIf(nNodes, pMap, mDanglingMap, nNodes);
//        cout << "Begin computation " << step << endl;

        std::future<double> futures[WORKERS];
        nodeId_t nodePerWorkers = nNodes / WORKERS;
        for(int w = 0; w < WORKERS; w++) {
            nodeId_t start = w * nodePerWorkers;
            nodeId_t end = w == WORKERS-1 ? nNodes : (w+1) * nodePerWorkers;
            futures[w] = std::async(launch::async, multiply, newPMap, start, end, mtRowMap, mtColMap, nullptr, pMap,
                                    d / pSum, telep + dangSum);
        }
        newPSum = 0;
        for(int w = 0; w < WORKERS; w++)
            newPSum += futures[w].get();

        err = error(nNodes, pMap, newPMap);

        //Close old file
        close(pFile);

        //Swap files
        pFile = newPFile;
        pMap = newPMap;
        pSum = newPSum;

//        cout << "error: " << err << endl;

        close(newPFile);
    }

    auto top = getVectorTopK(nNodes, k, pMap);

    close(pFile);

    return top;
}

int main(int argc, char *argv[])
{
    if(argc != 3)
    {
        cout << "Require 2 params, the graph path and the k value" << endl;
        return -1;
    }

    int maxK = atoi(argv[2]);
    string graphFile = argv[1];

    filesystem::remove_all(runtimeFolder);
    filesystem::create_directory(runtimeFolder);

    cout << "Begin sorting" << endl;
    Sorter s(graphFile, runtimeFolder, arcPerBlock);

    cout << "Begin writing" << endl;
    string mColFilename = runtimeFolder + string("/mCol.txt");
    string mRowFilename = runtimeFolder + string("/mRow.txt");
    string mValueFilename = runtimeFolder + string("/mValue.txt");
    string mtColFilename = runtimeFolder + string("/mtCol.txt");
    string mtRowFilename = runtimeFolder + string("/mtRow.txt");
    string mtValueFilename = runtimeFolder + string("/mtValue.txt");
    string mDanglingFilename = runtimeFolder + string("/mDanglings.txt");
    cout << " - Matrix" << endl;
    nodeId_t nNodes = writeMatrixFiles(s, false, mColFilename, mRowFilename, mValueFilename, mDanglingFilename);
    cout << " - Transposed Matrix" << endl;
    writeMatrixFiles(s, true, mtColFilename, mtRowFilename, mtValueFilename, mDanglingFilename);


    int mColFile = open(mColFilename.c_str(), O_RDONLY);
    int mRowFile = open(mRowFilename.c_str(), O_RDONLY);
    int mtColFile = open(mtColFilename.c_str(), O_RDONLY);
    int mtRowFile = open(mtRowFilename.c_str(), O_RDONLY);
    int mtValueFile = open(mtValueFilename.c_str(), O_RDONLY);
    int mDanglingFile = open(mDanglingFilename.c_str(), O_RDONLY);
    ull *mRowMap = (ull *) mmap(nullptr, nNodes*sizeof(ull), PROT_READ, MAP_SHARED, mRowFile, 0);
    nodeId_t *mColMap = (nodeId_t *) mmap(nullptr, mRowMap[nNodes-1]*sizeof(nodeId_t), PROT_READ, MAP_SHARED, mColFile, 0);
    ull *mtRowMap = (ull *) mmap(nullptr, nNodes*sizeof(ull), PROT_READ, MAP_SHARED, mtRowFile, 0);
    double *mtValueMap = (double *) mmap(nullptr, mtRowMap[nNodes-1]*sizeof(double), PROT_READ, MAP_SHARED, mtValueFile, 0);
    nodeId_t *mtColMap = (nodeId_t *) mmap(nullptr, mtRowMap[nNodes-1]*sizeof(nodeId_t), PROT_READ, MAP_SHARED, mtColFile, 0);
    bool *mDanglingMap = (bool *) mmap(nullptr, nNodes*sizeof(bool), PROT_READ, MAP_SHARED, mDanglingFile, 0);

    cout << "Pagerank ";
    auto startPR = std::chrono::high_resolution_clock::now();
    auto pagerankRes = pagerank(maxK, 0.85, nNodes, mtRowMap, mtValueMap, mtColMap, mDanglingMap);
    auto endPR = std::chrono::high_resolution_clock::now();
    auto millisPR = std::chrono::duration_cast<std::chrono::milliseconds>(endPR-startPR).count();
    cout << millisPR << endl;


    cout << "HITS ";
    auto startH = std::chrono::high_resolution_clock::now();
    auto hitsRes = hits(maxK, nNodes, mRowMap, mColMap, mtRowMap, mtColMap);
    auto endH = std::chrono::high_resolution_clock::now();
    auto millisH = std::chrono::duration_cast<std::chrono::milliseconds>(endH-startH).count();
    cout << millisH << endl;

    auto topA = hitsRes.first;
    auto topH = hitsRes.second;
//    cout<< "Top A" << endl;
//    for(auto i : topA)
//        cout << i.first << " " << i.second << endl;
//    cout<< "Top H" << endl;
//    for(auto i : topH)
//        cout << i.first << " " << i.second << endl;

    cout << "In degree ";
    auto startIn = std::chrono::high_resolution_clock::now();
    auto inDegreeRes = inDegree(maxK, nNodes, mtRowMap);
    auto endIn = std::chrono::high_resolution_clock::now();
    auto millisIn = std::chrono::duration_cast<std::chrono::milliseconds>(endIn-startIn).count();
    cout << millisIn << endl;
//    cout<< "Top" << endl;
//    for(auto i : inDegreeRes)
//        cout << i.first << " " << i.second << endl;

    for(int k = 10; k <= maxK; k+=10)
    {
        double jaccPrHt = jaccardIndex(getSetOfFirstK(k, pagerankRes), getSetOfFirstK(k, topA));
        double jaccPrIn = jaccardIndex(getSetOfFirstK(k, pagerankRes), getSetOfFirstK(k, inDegreeRes));
        double jaccHtIn = jaccardIndex(getSetOfFirstK(k, topA), getSetOfFirstK(k, inDegreeRes));
        cout << "------- Jaccard K = " << k << " -------" << endl;
        cout << "PageRank - Hits:      " << jaccPrHt << endl;
        cout << "PageRank - In Degree: " << jaccPrIn << endl;
        cout << "Hits     - In Degree: " << jaccHtIn << endl;
    }
    cout << "-------------------------------" << endl;

    /*
    for(int k = 10; k <= maxK; k+=10)
    {
        double jaccPrHt = jaccardIndex(getSetOfFirstK(k, pagerankRes), getSetOfFirstK(k, topA));
        cout << "(" << k << "," << jaccPrHt << ")";
    }
    cout << endl;
    for(int k = 10; k <= maxK; k+=10)
    {
        double jaccPrIn = jaccardIndex(getSetOfFirstK(k, pagerankRes), getSetOfFirstK(k, inDegreeRes));
        cout << "(" << k << "," << jaccPrIn << ")";
    }
    cout << endl;
    for(int k = 10; k <= maxK; k+=10)
    {
        double jaccHtIn = jaccardIndex(getSetOfFirstK(k, topA), getSetOfFirstK(k, inDegreeRes));
        cout << "(" << k << "," << jaccHtIn << ")";
    }
    cout << endl;
    */

    close(mRowFile);
    close(mColFile);
    close(mtRowFile);
    close(mtColFile);

    return 0;
}