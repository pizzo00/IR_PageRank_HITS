#include "Sorter.h"

Sorter::Sorter(string const& inputFilename, string const& blocksFolder, int arcPerBlock)
: blocksFolder(blocksFolder)
{
    filesystem::remove_all(blocksFolder);
    filesystem::create_directory(blocksFolder);
    ifstream input(inputFilename);

    vector<pair<nodeId_t, nodeId_t>> block;
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
            block.emplace_back(dst, src); //Reversed for correct sorting
            if(++arcInBlock >= arcPerBlock)
            {
                sort(block.begin(), block.end());
                ofstream blockFile;
                blockFile.open(format("{}\\{}.txt", blocksFolder, blockNum));
                for(auto const& i : block)
                {
                    blockFile << i.second << "\t" << i.first << endl;
                }
                blockFile.close();

                block.clear();
                arcInBlock = 0;
                blockNum++;
            }
        }
    }

    blocksCount = blockNum;
    blockStreams.resize(blocksCount);
    front.resize(blocksCount);
    frontAvailable.resize(blocksCount);
    for(int i = 0; i < blocksCount; i++)
    {
        blockStreams[i] = ifstream(format("{}\\{}.txt", blocksFolder, i));
        readLineOfStream(i);
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
        value = make_pair(front[blockNum].second, front[blockNum].first);
        readLineOfStream(blockNum);
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
        front[index] = make_pair(dst, src);
    }
    else
    {
        frontAvailable[index] = false;
    }
}
