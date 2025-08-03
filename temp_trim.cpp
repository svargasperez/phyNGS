#include <cstring>
#include <string>
#include <iostream>
#include <algorithm>
#include <functional>
#include <stdlib.h>
#include <stdio.h>
#include "defs.h"

using std::string;

void trim(char *trimmed, const char *dna_seq, int32 seq_len, const string pat)
{
    // TODO: What should the final logic be regarding small sequences?
    // Do not trim if sequence is smaller than pattern
    if (seq_len < pat.size())
    {
        printf("Sequence length (%d) shorter than pattern (%lu)\n", seq_len, pat.size());
        trimmed = strcpy(trimmed, dna_seq);
        return;
    }

    // TODO: Possible optimization with tracking a new_len variable instead of end_pos
    int32 new_start_pos = 0;     // Inclusive start index
    int32 new_end_pos = seq_len; // Exclusive end index

    // Check if beginning of dna sequence has the pattern
    if (pat.compare(0, pat.size(), dna_seq, 0, pat.size()) == 0)
    {
        printf("5' end contains adapter sequence\n");
        new_start_pos = pat.size();
    }
    else
        printf("No 5' match with adapter sequence\n");

    // Check if end of dna sequence has the pattern
    if (pat.compare(0, pat.size(), dna_seq, seq_len - pat.size(), pat.size()) == 0)
    {
        printf("3' end contains adapter sequence\n");
        new_end_pos = seq_len - pat.size();
    }
    else
        printf("No 3' match with adapter sequence\n");

    // trimmed = dna_seq[new_start_pos, new_length]
    int32 new_length = new_end_pos - new_start_pos;
    strncpy(trimmed, dna_seq + new_start_pos, new_length);
    trimmed[new_length] = '\0';
}

int main(int argc, char **argv)
{
    char *dna_seq = argv[1];
    int32 seq_len = strlen(dna_seq);

    string pat = "TATA";

    // char *it = std::search(dna_seq, dna_seq+seq_len, std::boyer_moore_searcher(pat, pat+pat_len));
    // std::cout << "'" << it << "'" << "\n";

    // Possibly overallocate on stack b/c trimmed seq will never exceed full seq
    char trimmed_dna_seq[seq_len];

    trim(trimmed_dna_seq, dna_seq, seq_len, pat);
    printf("Trimmed sequence: %s\n", trimmed_dna_seq);

    exit(EXIT_SUCCESS);
}