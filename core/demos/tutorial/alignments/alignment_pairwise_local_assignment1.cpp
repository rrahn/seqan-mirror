//FRAGMENT(main)
#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
//FRAGMENT(init)
	Align< String<AminoAcid> > ali;
	resize(rows(ali), 1);
	assignSource(row(ali, 0), "PNCFDAKQRTASRPL");
	assignSource(row(ali, 1), "CFDKQKNNRTATRDTA");

//FRAGMENT(ali)
	LocalAlignmentFinder<> finder(ali);
	Score<int> sc(3,-2,-1,-5);
	unsigned count = 0;
	LocalAlignmentEnumerator<Score<int>, Unbanded> enumerator(sc);
	while (nextLocalAlignment(ali, enumerator) && count < 3)
	{
		::std::cout << "Score = " << getScore(finder) << ::std::endl;
		::std::cout << ali;
		++count;
	}
	return 0;
}
