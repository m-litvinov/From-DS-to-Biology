
def PatternCount(Text:str, Pattern:str):
    """Counts occurencies of Pattern in Text (with overlap)"""
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count += 1
    return count

def ApproximatePatternCount(Text, Pattern, d):
    """Counts occurencies of Pattern (with the most d mismatches) in Text (with overlap)"""
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
            count += 1
    return count

def FrequencyMap(Text:str, k:int):
    """Return dictionary with frequency of all k-mers in Text"""
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = 0
    for pattern in freq.keys():
        for i in range(n-k+1):
            if Text[i:i+k] == pattern:
                freq[pattern] += 1
    return freq

def ApproximateFrequencyMap(Text:str, k:int, d:int, complementary=True):
    """Return dictionary with frequency of all k-mers in Text, with the most d mismatches and optional complementary"""
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = 0
    for pattern in freq.keys():
        for i in range(n-k+1):
            if complementary:
                if HammingDistance(Text[i:i+k], pattern) <= d or HammingDistance(Text[i:i+k], ReverseComplement(pattern)) <= d:
                    freq[pattern] += 1
            else:
                if HammingDistance(Text[i:i+k], pattern) <= d:
                    freq[pattern] += 1

    return freq

def FrequentWords(Text:str, k:int):
    """Return the most frequent k-mers in Text"""
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for word in freq:
        if freq[word] == m:
        # if freq[word] >= 3:
            words.append(word)
    return words

def ApproximateFrequentWords(Text:str, k:int, d:int, **kwargs):
    """Returns the most frquent k-mers with the most d mismatches in Text"""
    words = []
    freq = ApproximateFrequencyMap(Text, k, d, **kwargs)
    m = max(freq.values())
    for word in freq:
        if freq[word] == m:
        # if freq[word] >= 3:
            words.append(word)
    return words

def Reverse(Pattern):
    Pattern = Pattern[::-1]
    return Pattern

def Complement(Pattern):
    complement_map = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
        "a": "t",
        "t": "a",
        "c": "g",
        "g": "c",

    }
    complement_str = ''
    for letter in Pattern:
        if letter not in complement_map:
            raise Exception("DNA Error: Unknown nucleotide")
        complement_str += complement_map[letter]
    return complement_str

def ReverseComplement(Pattern:str):
    """Create complementary and reverse string based on Pattern"""
    Pattern = Reverse(Pattern)
    Pattern = Complement(Pattern)
    return Pattern

def PatternMatching(Text:str, Pattern:str):
    """Return list with starting indexes where Pattern occures in Text"""
    positions = []
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            positions.append(i)
    return positions

def ApproximatePatternMatching(Text:str, Pattern:str, d:int):
    """Return list with starting indexes where Pattern (with the most d mismatches) occures in Text"""
    positions = []
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
            positions.append(i)
    return positions

def SymbolArray(Genome:str, symbol:str):
    """Return an array of number with occurrences of symbol in half-genome window"""
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    array[0] = PatternCount(symbol, Genome[0:n//2])

    for i in range(1, n):
        array[i] = array[i-1]
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array

def SkewArray(Genome:str):
    """Create list (skew array) which shows difference between total amount of guanine and cytosine"""
    Skew = []
    Skew.append(0)
    n = len(Genome)
    for i in range(n):
        if Genome[i].upper() == 'A' or Genome[i].upper() == 'T':
            Skew.append(Skew[i])
        elif Genome[i].upper() == 'G':
            Skew.append(Skew[i]+1)
        elif Genome[i].upper() == 'C':
            Skew.append(Skew[i]-1)
        else:
            raise Exception("DNA Error: Unknown nucleotide")
    
    return Skew

def MinimumSkew(Genome:str):
    """Return list of positions in Genome where difference between total amount of guanine and cytosine is minimum"""
    skew = SkewArray(Genome)
    min_value = min(skew)

    positions = []
    for i in range(len(skew)):
        if skew[i] == min_value:
            positions.append(i)
    
    return positions

def HammingDistance(p:str, q:str):
    """Calculate Hamming distance (count of letter mismatches) between two string"""
    h_distance = 0
    if len(p) != len(q):
        h_distance += abs(len(p) - len(q))
    
    for i in range(len(p)):
        try:
            if p[i] != q[i]:
                h_distance += 1
        except IndexError:
            pass
    return h_distance

