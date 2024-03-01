#include "Sorter.h"

Sorter::Sorter(string const& inputFilename, string const& runtimeFolder, int arcPerBlock)
: runtimeFolder(runtimeFolder)
{
    auto saveBlock = [&runtimeFolder](vector<pair<nodeId_t, nodeId_t> > & block, int const blockNum)
    {
        sort(block.begin(), block.end());
        ofstream blockFile;
        blockFile.open(runtimeFolder + string("/block_") + to_string(blockNum) + string(".txt"));
        for(auto const& i : block)
        {
            blockFile << i.first << "\t" << i.second << endl;
        }
        blockFile.close();

        //Transposed
        sort(block.begin(), block.end(),
             [](pair<nodeId_t, nodeId_t> const& l, pair<nodeId_t, nodeId_t> const& r)
             {
                 return make_pair(l.second, l.first) < make_pair(r.second, r.first);
             }
        );
        ofstream blockTFile;
        blockTFile.open(runtimeFolder + string("/blockT_") + to_string(blockNum) + string(".txt"));
        for(auto const& i : block)
        {
            blockTFile << i.first << "\t" << i.second << endl;
        }
        blockTFile.close();

        block.clear();
    };

    ifstream input(inputFilename);

    vector<pair<nodeId_t, nodeId_t> > block;
    int arcInBlock = 0;
    int blockNum = 0;
    string line;
    while(getline(input, line))
    {
        if(line[0] != '#')
        {
            auto sep = line.find('\t');
            string strSrc = line.substr(0, sep);
            string strDst = line.substr(sep+1);
            nodeId_t src = stoul(strSrc);
            nodeId_t dst = stoul(strDst);
            block.emplace_back(src, dst);
            if(++arcInBlock >= arcPerBlock)
            {
                saveBlock(block, blockNum);

                arcInBlock = 0;
                blockNum++;
            }
        }
    }

    saveBlock(block, blockNum);

    arcInBlock = 0;
    blockNum++;

    blocksCount = blockNum;
    blockStreams.resize(blocksCount);
    blockTStreams.resize(blocksCount);
    front.resize(blocksCount);
    frontT.resize(blocksCount);
    frontAvailable.resize(blocksCount);
    frontTAvailable.resize(blocksCount);
    for(int i = 0; i < blocksCount; i++)
    {
        blockStreams[i] = ifstream(runtimeFolder + string("/block_") + to_string(i) + string(".txt"));
        blockTStreams[i] = ifstream(runtimeFolder + string("/blockT_") + to_string(i) + string(".txt"));
        readLineOfStream(i);
        readLineOfStreamT(i);
    }
}

bool Sorter::next(pair<unsigned int, unsigned int> &value)
{
    int blockNum = -1;
    for(int i = 0; i < blocksCount; i++)
    {
        if(frontAvailable[i] && (blockNum == -1 || front[i] < front[blockNum]))
            blockNum = i;
    }

    if(blockNum != -1)
    {
        value = front[blockNum];
        readLineOfStream(blockNum);
    }

    return blockNum != -1;
}

bool Sorter::nextTransposed(pair<unsigned int, unsigned int> &value)
{
    int blockNum = -1;
    for(int i = 0; i < blocksCount; i++)
    {
        if(frontTAvailable[i] && (blockNum == -1 ||
        make_pair(frontT[i].second, frontT[i].first) < make_pair(frontT[blockNum].second, frontT[blockNum].first)))
            blockNum = i;
    }

    if(blockNum != -1)
    {
        value = frontT[blockNum];
        readLineOfStreamT(blockNum);
    }

    return blockNum != -1;
}

void Sorter::readLineOfStream(int index)
{
    string line;
    if(getline(blockStreams[index], line))
    {
        frontAvailable[index] = true;

        auto sep = line.find('\t');
        string strSrc = line.substr(0, sep);
        string strDst = line.substr(sep + 1);
        nodeId_t src = stoul(strSrc);
        nodeId_t dst = stoul(strDst);
        front[index] = make_pair(src, dst);
    }
    else
    {
        frontAvailable[index] = false;
    }
}

void Sorter::readLineOfStreamT(int index)
{
    string line;
    if(getline(blockTStreams[index], line))
    {
        frontTAvailable[index] = true;

        auto sep = line.find('\t');
        string strSrc = line.substr(0, sep);
        string strDst = line.substr(sep + 1);
        nodeId_t src = stoul(strSrc);
        nodeId_t dst = stoul(strDst);
        frontT[index] = make_pair(src, dst);
    }
    else
    {
        frontTAvailable[index] = false;
    }
}
