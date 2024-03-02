#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <queue>
#include <set>
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

template <typename T, typename T2>
set<T> getSetOfFirst(vector<pair<T, T2>> const& in)
{
    set<T> out;
    for(auto const& i : in)
        out.insert(i.first);
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

//    cout << " - A" << endl;
    writeVector(aFilename, nNodes, (double)1.0);
//    cout << " - H" << endl;
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

//        cout << "Begin computation " << step << endl;
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

int main() {
    filesystem::remove_all(runtimeFolder);
    filesystem::create_directory(runtimeFolder);

    Sorter s("/home/pizzolato/IR_Project/data/web-NotreDame.txt", runtimeFolder, arcPerBlock);

    cout << "Begin writing" << endl;
    string mColFilename = runtimeFolder + string("/mCol.txt");
    string mRowFilename = runtimeFolder + string("/mRow.txt");
    string mtColFilename = runtimeFolder + string("/mtCol.txt");
    string mtRowFilename = runtimeFolder + string("/mtRow.txt");
    cout << " - Matrix" << endl;
    nodeId_t nNodes = writeMatrixFiles(s, false, mColFilename, mRowFilename);
    cout << " - Transposed Matrix" << endl;
    writeMatrixFiles(s, true, mtColFilename, mtRowFilename);


    int mRowFile = open(mRowFilename.c_str(), O_RDONLY);
    int mColFile = open(mColFilename.c_str(), O_RDONLY);
    int mtRowFile = open(mtRowFilename.c_str(), O_RDONLY);
    int mtColFile = open(mtColFilename.c_str(), O_RDONLY);
    ull *mRowMap = (ull *) mmap(nullptr, nNodes*sizeof(ull), PROT_READ, MAP_SHARED, mRowFile, 0);
    nodeId_t *mColMap = (nodeId_t *) mmap(nullptr, mRowMap[nNodes-1]*sizeof(nodeId_t), PROT_READ, MAP_SHARED, mColFile, 0);
    ull *mtRowMap = (ull *) mmap(nullptr, nNodes*sizeof(ull), PROT_READ, MAP_SHARED, mtRowFile, 0);
    nodeId_t *mtColMap = (nodeId_t *) mmap(nullptr, mtRowMap[nNodes-1]*sizeof(nodeId_t), PROT_READ, MAP_SHARED, mtColFile, 0);

    for(int k = 10; k <= 100; k+=10) {
        auto hitsRes = hits(k, nNodes, mRowMap, mColMap, mtRowMap, mtColMap);

        auto topA = hitsRes.first;
        auto topH = hitsRes.second;
//    cout<< "Top A" << endl;
//    for(auto i : topA)
//        cout << i.first << " " << i.second << endl;
//    cout<< "Top H" << endl;
//    for(auto i : topH)
//        cout << i.first << " " << i.second << endl;

        auto inDegreeRes = inDegree(k, nNodes, mtRowMap);
//    cout<< "Top" << endl;
//    for(auto i : inDegreeRes)
//        cout << i.first << " " << i.second << endl;

        double jacc = jaccardIndex(getSetOfFirst(topA), getSetOfFirst(inDegreeRes));
        cout << "Jaccard " << k << ": " << jacc << endl;
    }

    close(mRowFile);
    close(mColFile);
    close(mtRowFile);
    close(mtColFile);

    return 0;
}