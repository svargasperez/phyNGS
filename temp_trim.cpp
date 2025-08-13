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
    const char *seq_end = dna_seq + seq_len;

    // TODO: What should the final logic be regarding small sequences?
    // Do not trim if sequence is smaller than pattern
    // TODO: Fix when seq_len == pat.size() (or other edge cases, check start/end indicies at end)
    if (seq_len < pat.size())
    {
        printf("Sequence length (%d) shorter than pattern (%lu)\n", seq_len, pat.size());
        trimmed = strcpy(trimmed, dna_seq);
        return;
    }

    // TODO: Possible optimization with tracking pointers instead of indicies
    int32 new_start_pos = 0;     // Inclusive start index
    int32 new_end_pos = seq_len; // Exclusive end index

    // Check if start of dna sequence has the pattern
    // TODO: Improve naive algorithm (KMP?)
    const char *search_start = strstr(dna_seq, pat.c_str());
    int32 index = search_start ? search_start - dna_seq : -1;
    printf("First index of pattern in string: %d\n", index);

    // Trim if end of dna sequence has the pattern
    if (search_start)
    {
        printf("5' end contains adapter sequence\n");
        new_start_pos = std::distance(dna_seq, search_start) + pat.size();
    }
    else
        printf("No 5' match with adapter sequence\n");

    // Check if end of dna sequence has the pattern
    int32 search_offset = MAX(seq_len / 2, pat.size());
    auto bm_searcher = std::boyer_moore_searcher(pat.begin(), pat.end());
    const char *search_end = std::search(seq_end - search_offset, seq_end, bm_searcher);

    if (search_end != seq_end)
        index = std::distance(dna_seq, search_end);
    else
        index = -1;
    printf("Last index of pattern in string: %d\n", index);

    // Trim if end of dna sequence has the pattern
    if (search_end != seq_end)
    {
        printf("3' end contains adapter sequence\n");
        new_end_pos = std::distance(dna_seq, search_end);
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

    string pat = argv[2];
    // std::cout << "'" << it << "'" << "\n";

    // Possibly overallocate on stack b/c trimmed seq will never exceed full seq
    char trimmed_dna_seq[seq_len];

    trim(trimmed_dna_seq, dna_seq, seq_len, pat);
    printf("Trimmed sequence: %s\n", trimmed_dna_seq);

    exit(EXIT_SUCCESS);
}