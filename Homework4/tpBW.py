# >>===================================================<<
# ||+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+||
# ||--------------------|B|O|X|-----------------------|||
# |||B|u|r|r|o|w|s|-|W|h|e|e|l|e|r|-|T|r|a|n|s|f|o|r|m|||
# ||+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+||
# >>===================================================<<

# Template by Leo Ackermann (2025)
# Code by Samuel Fosse
import copy
from ks import simple_kark_sort
import numpy as np

# RLE compression
def rle(seq):
    l = 0
    r = 0
    ans = ""
    while(l<len(seq)):
        while(r<len(seq) and seq[l]==seq[r]):
            r+=1
        ans = "".join([ans, seq[l] + str(r-l)])
        l = r
    return ans
# >>==========================================================================<<
# >>                              FM-index class                              <<
# >>==========================================================================<<


# Implements a BW transform, along with functions to reverse the transform, and use it to do word search

class FMindex:

    # >>===[ Class constructor ]================================================
    def __init__(self, seq, verbose=False):
        self.sa = simple_kark_sort(seq)
        self.bwt = self.__get_bwt(seq)
        self.fm_rank = self.__get_fm_rank()
        self.count, self.next_smallest_letter = self.__get_count_smallest_letter()
        self.fm_ranks = self.__get_fm_ranks()
    
    # >>===[ Attribute initialisation functions ]===============================
    def __get_bwt(self, seq):
        return [seq[-1]] + [seq[self.sa[i]-1] if (self.sa[i]>0) else "$" for i in range(len(seq))]

    def __get_fm_rank(self):
        last_seen = [0]*(ord("z")+1)
        answer = [0]*len(self.bwt)
        for i in range(len(self.bwt)):
            last_seen[ord(self.bwt[i])] += 1
            answer[i] = last_seen[ord(self.bwt[i])]
        return answer
    
    def __get_count_smallest_letter(self):
        count = {}
        firstletters = sorted(self.bwt)
        current_pos = -1
        j = 0
        next_smallest_letter = {}
        while(j<len(firstletters)):
            current_letter = firstletters[j]
            current_pos = j
            while(j<len(firstletters) and firstletters[j]==current_letter):
                count[current_letter] = current_pos
                j += 1
            if(j<len(firstletters)):
                next_smallest_letter[current_letter] = firstletters[j]
            else:
                next_smallest_letter[current_letter] = "~"

        count["~"] = j

        return count, next_smallest_letter

    def __get_fm_ranks(self):
        fm_ranks = np.zeros((ord("z") + 1, len(self.bwt)), dtype = int)

        for character in range((ord("z") + 1)):
            times_seen = 0
            for i in range(len(self.bwt)):
                if(ord(self.bwt[i])==character):
                    times_seen += 1
                fm_ranks[character][i] = times_seen
        
        return fm_ranks
    
    # >>===[ Pattern matching functions ]=======================================

    def membership(self, pattern):
        i_min = self.count[pattern[-1]]
        i_max = self.count[self.next_smallest_letter[pattern[-1]]]-1
        for k in range(len(pattern)-2, -1, -1):
            i_min = self.count[pattern[k]] + self.fm_ranks[ord(pattern[k])][i_min - 1]
            i_max = self.count[pattern[k]] + self.fm_ranks[ord(pattern[k])][i_max] - 1
            if(i_min > i_max):
                return False

        return True

    def counts(self, pattern):
        i_min = self.count[pattern[-1]]
        i_max = self.count[self.next_smallest_letter[pattern[-1]]]-1
        for k in range(len(pattern)-2, -1, -1):
            i_min = self.count[pattern[k]] + self.fm_ranks[ord(pattern[k])][i_min - 1]
            i_max = self.count[pattern[k]] + self.fm_ranks[ord(pattern[k])][i_max] - 1
            if(i_min > i_max):
                return 0

        return i_max - i_min + 1
    
    def locate(self, pattern):
        i_min = self.count[pattern[-1]]
        i_max = self.count[self.next_smallest_letter[pattern[-1]]]-1
        for k in range(len(pattern)-2, -1, -1):
            i_min = self.count[pattern[k]] + self.fm_ranks[ord(pattern[k])][i_min - 1]
            i_max = self.count[pattern[k]] + self.fm_ranks[ord(pattern[k])][i_max] - 1
            if(i_min > i_max):
                return []

        seen_pos = []
        for i in range(i_max-i_min + 1):
            seen_pos.append(self.sa[i_min - 1 + i])
        return sorted(seen_pos)


    def get_string__naive(self):
        bw_matrix = copy.deepcopy(self.bwt)
        bw_matrix.sort()
        for i in range(len(self.bwt)-1):
            bw_matrix = ["".join([self.bwt[i], bw_matrix[i]]) for i in range(len(self.sa) + 1)]
            bw_matrix.sort()
        return bw_matrix[0][1:]

    def get_string(self):
        S = "$"
        i_row = 0
        for i in range(len(self.bwt)-1):
            prec = self.bwt[i_row]
            S = prec + S
            i_row = self.lfmapping(i_row)
        return S[:-1]

    def lfmapping(self, i):
        return self.count[self.bwt[i]] + self.fm_rank[i] - 1
    

ex = FMindex("BANANA")
print(ex.sa)
print(ex.fm_rank)
print(ex.count)
print(ex.get_string__naive())
print(rle("BANANA"))
print(ex.get_string())
print(ex.locate("AB"))
print(ex.next_smallest_letter)
print(ex.fm_ranks[ord("B")])
print(rle("".join(ex.bwt)))

