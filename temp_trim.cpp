#include <cstring>
#include <string>
#include <iostream>
#include <algorithm>
#include <functional>
#include <stdlib.h>
#include <stdio.h>
#include "defs.h"

using std::string;

/**
 * @brief Attempts to find part of a pattern on an end of text,
 *        as if the pattern continued off the edge of the text.
 *        At least half of the pattern must be present.
 *        
 *        Example: text="o bar baz qux" with pat="foo bar" (left side);
 *        would return 5 because "(fo)o bar". has 5/7 characters.
 * 
 * @param text C string of the text to be searched.
 * @param text_end Pointer to the end of the C string to be searched.
 * @param pat Pattern to search for in the C string text.
 * @param left_side Flag, whether to search the left or right side.
 * @return The length of the partial pattern within the text, or -1.
 */
int32 partial_search(const char *text, const char *text_end, const string pat, bool left_side)
{
    if (left_side)
        // Try to match pat[i:] to beginning of text
        for (int32 i = 1; i <= pat.size() / 2; i++)
        {
            int32 len = pat.size() - i;
            if (pat.compare(i, len, text, len) == 0)
                return len;
        }
    else
        // Try to match pat[0:len] to end of text
        for (int32 i = 1; i <= pat.size() / 2; i++)
        {
            int32 len = pat.size() - i;
            if (pat.compare(0, len, text_end - len, len) == 0)
                return len;
        }

    return -1; // No match
}

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

    int32 index;
    int32 partial_len;
    int32 search_offset = MAX(seq_len / 2, pat.size());

    // Check if start of dna sequence has the pattern
    auto bm_searcher = std::boyer_moore_searcher(pat.begin(), pat.end());
    const char *search_start = std::search(dna_seq, dna_seq + search_offset, bm_searcher);
    if (search_start != dna_seq + search_offset)
        index = std::distance(dna_seq, search_start);
    else
        index = -1;
    printf("  First index of pattern in string: %d\n", index);

    // Trim if start of dna sequence has the pattern
    if (search_start != dna_seq + search_offset)
    {
        printf("5' end contains adapter sequence (index %d)\n", index);
        new_start_pos = std::distance(dna_seq, search_start) + pat.size();
    }
    // Trim if dna sequence has part of the pattern
    else if ((partial_len = partial_search(dna_seq, seq_end, pat, true)) != -1)
    {
        printf("5' end contains partial adapter sequence (%d/%lu)\n",
               partial_len, pat.size());
        new_start_pos = partial_len;
    }
    else
        printf("No 5' match with adapter sequence\n");

    // Check if end of dna sequence has the pattern
    // TODO: Improve naive algorithm (possibly custom right-most Boyer Moore)
    const char *search_end = std::find_end(seq_end - search_offset, seq_end, pat.begin(), pat.end());

    if (search_end != seq_end)
        index = std::distance(dna_seq, search_end);
    else
        index = -1;
    printf("  Last index of pattern in string: %d\n", index);

    // Trim if end of dna sequence has the pattern
    if (search_end != seq_end)
    {
        printf("3' end contains adapter sequence (index %d)\n", index);
        new_end_pos = std::distance(dna_seq, search_end);
    }
    else if ((partial_len = partial_search(dna_seq, seq_end, pat, false)) != -1)
    {
        printf("3' end contains partial adapter sequence (%d/%lu)\n",
               partial_len, pat.size());
        new_end_pos = seq_len - partial_len;
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

    // Possibly overallocate on stack b/c trimmed seq will never exceed full seq
    char trimmed_dna_seq[seq_len];

    trim(trimmed_dna_seq, dna_seq, seq_len, pat);
    printf("Trimmed sequence: %s\n", trimmed_dna_seq);

    exit(EXIT_SUCCESS);
}