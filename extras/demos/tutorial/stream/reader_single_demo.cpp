#include <iostream>
#include <fstream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;

int main(int argc, char const ** argv)
{
    if (argc != 2)
        return 1;
    std::fstream stream(argv[1], std::ios::binary | std::ios::in);
    if (!stream.good())
        return 1;
    
    RecordReader<std::fstream, SinglePass<> > reader(stream);
    StringSet<CharString> result;
    
    while (!atEnd(reader))
    {
        resize(result, length(result) + 1);
        int res = readLine(back(result), reader);
        if (res != 0)
            return 1;
    }

    return 0;
}
