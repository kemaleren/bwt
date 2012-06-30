"""
bwt.py
Author: Kemal Eren

Functions for transforming and searching strings using the
Burrows-Wheeler transform.

"""

EOS = "\0"

def get_bwt(s):
    """
    Returns the Burrows-Wheeler transform of 's'.

    Based on the Python implementation given on Wikipedia:
    http://en.wikipedia.org/wiki/Burrows%E2%80%93Wheeler_transform

    Examples:
    ---------

    >>> get_bwt('banana\0')
    'annb\x00aa'

    """
    table = [s[i:] + s[:i] for i in range(len(s))]  # rotations of string
    table = sorted(table)
    last_column = [row[-1] for row in table]  # Last characters of each row
    return "".join(last_column)  # Convert list of characters into string


def get_occ(bwt):
    """
    Returns occurrence information for letters in the string 'bwt'.
    occ[letter][i] = the number of occurrences of 'letter' in
    bwt[0, i + 1].

    Examples:
    ---------

    >>> get_occ('annb\x00aa')
    {'\x00': [0, 0, 0, 0, 1, 1, 1],
    'a': [1, 1, 1, 1, 1, 2, 3],
    'b': [0, 0, 0, 1, 1, 1, 1],
    'n': [0, 1, 2, 2, 2, 2, 2]}

    """
    letters = set(bwt)
    occ = {}
    for letter in letters:
        occ[letter] = []
        for i in range(len(bwt)):
            occ[letter].append(len([j for j in bwt[:i + 1] if j == letter]))
    return occ


def get_count(s):
    """
    Returns count information for the letters in the string 's'.

    count[letter] contains the number of symbols in 's' that are
    lexographically smaller than 'letter'.

    Examples:
    ---------

    >>> get_count('sassy')
    {'a': 0, 's': 1, 'y': 4}
    
    """
    letters = set(s)
    count = {}
    for letter in letters:
        count[letter] = len([i for i in s if i < letter])
    return count


def get_sa(s):
    """
    Returns the suffix array of 's'

    Examples:
    ---------
    
>>> get_sa('banana\0')
    [6, 5, 3, 1, 0, 4, 2]

    """
    suffixes = [s[i:] for i in range(len(s))]
    return [suffixes.index(i) for i in sorted([s[j:] for j in range(len(s))])]
    

def use_occ(occ, letter, i, length):
    """
    Handles retrieving occurrence information; in particular, deals
    with overflows.

    """
    if i == -1:
        return 0
    if i == length:
        return occ[letter][-1]
    return occ[letter][i]


def bwt_interval(query, occ, count, length):
    """
    Returns the interval [start, stop) in the suffix array that
    corresponds to exact matches of 'query'.

    occ    : occurrence information calculated with get_occ.
    count  : count information calculated with get_count.
    length : the number of characters in the reference string.

    """
    begin = 0
    end = length - 1
    query = query[::-1] #reverse the query
    
    for letter in query:
        begin = count[letter] + use_occ(occ, letter, begin - 1, length) + 1 
        end = count[letter] + use_occ(occ, letter, end, length)
        if begin > end:
            return None, None
    return begin, end + 1


def mutations(s, dist, alphabet, used=None):
    """
    Yields all strings that can be derived from string 's' with
    Hamming distance at most 'dist', using a given alphabet.

    Examples:
    ---------
    
    >>> mutations('bad', 1, set(['a', 'b', 'd'])
    ['bad', 'aad', 'baa', 'bbd', 'bab', 'dad', 'bdd']

    """
    if used is None:
        used = set([])
    if dist == 0:
        if not s in used:
            used.add(s)
            yield s
    for d in range(dist): #from least to most distance
        for letter in alphabet:
            for pos in range(len(s)):
                for result in mutations(s[:pos] + letter + s[pos + 1:], dist - 1, alphabet, used):
                    yield result


def get_bwt_data(reference):
    """Returns the data structures needed to perform BWT searches"""
    alphabet = set(reference)
    assert EOS not in reference
    reference += EOS
    bwt = get_bwt(reference)
    occ = get_occ(bwt)
    count = get_count(reference[:-1])
    sa = get_sa(reference)
    return alphabet, bwt, occ, count, sa


def bwt_inexact_match(query, reference, mismatches=None):
    """
    Find all matches of the string 'query' in the string 'reference', with at most
    'mismatch' mismatches

    Examples:
    ---------
    
    >>> bwt_inexact_match('abc', 'abcabd', 1)
    [0, 3]

    """

    alphabet, bwt, occ, count, sa = get_bwt_data(reference)
    results = []

    for s in mutations(query, mismatches, alphabet):
        begin, end = bwt_interval(s, occ, count, len(bwt))
        if begin is not None:
            results += sa[begin:end]
    return sorted(set(results))
        

def bwt_exact_match(query, reference):
    """
    Find all exact matches of the string 'query' in the string 'reference'.

    Examples:
    ---------
    >>> bwt_exact_match('abc', 'abcabcabc')
    [0, 3, 6]
    
    >>> bwt_exact_match('gef', 'abcabcabc')
    []

    """
    return bwt_inexact_match(query, reference, mismatches=0)
    

