#include <cstring>
#include <string>
#include <iostream>
#include <algorithm>
#include <functional>
#include <stdlib.h>
#include <stdio.h>
#include "defs.h"

using std::string;

void trim(char *dna_seq, int32 seq_len, string pat);

int main(int argc, char **argv)
{
    char *dna_seq = argv[1];
    int32 seq_len = strlen(dna_seq);

    string pat = "TATA";

    // char *it = std::search(dna_seq, dna_seq+seq_len, std::boyer_moore_searcher(pat, pat+pat_len));
    // std::cout << "'" << it << "'" << "\n";

    trim(dna_seq, seq_len, pat);

    exit(EXIT_SUCCESS);
}

void trim(char *dna_seq, int32 seq_len, string pat)
{
    // Do not check 
    if (seq_len < pat.size())
    {
        printf("Sequence length (%d) shorter than pattern (%lu)\n", seq_len, pat.size());
        return;
    }

    // Check if beginning of dna sequence has the pattern
    if (pat.compare(0, pat.size(), dna_seq, 0, pat.size()) == 0)
        printf("5' end contains adapter sequence\n");
    else
        printf("No 5' match with adapter sequence\n");
    
    // Check if end of dna sequence has the pattern
    if (pat.compare(0, pat.size(), dna_seq, seq_len - pat.size(), pat.size()) == 0)
        printf("3' end contains adapter sequence\n");
    else
        printf("No 3' match with adapter sequence\n");
}