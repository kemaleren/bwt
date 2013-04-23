"""
bwt.py
Author: Kemal Eren

Functions for transforming and searching strings using the
Burrows-Wheeler transform.

"""

EOS = "\0"

def get_bwt(s):
    r"""
    Returns the Burrows-Wheeler transform of 's'.

    Based on the Python implementation given on Wikipedia:
    http://en.wikipedia.org/wiki/Burrows%E2%80%93Wheeler_transform

    Examples:
    ---------

    >>> get_bwt('banana\0')
    'annb\x00aa'

    """
    table = sorted([s[i:] + s[:i] for i in range(len(s))])
    last_column = [row[-1] for row in table]
    return "".join(last_column)


def get_occ(bwt):
    r"""
    Returns occurrence information for letters in the string 'bwt'.
    occ[letter][i] = the number of occurrences of 'letter' in
    bwt[0, i + 1].

    Examples:
    ---------

    >>> get_occ('annb\x00aa')['\x00']
    [0, 0, 0, 0, 1, 1, 1]

    >>> get_occ('annb\x00aa')['a']
    [1, 1, 1, 1, 1, 2, 3]

    >>> get_occ('annb\x00aa')['b']
    [0, 0, 0, 1, 1, 1, 1]

    >>> get_occ('annb\x00aa')['n']
    [0, 1, 2, 2, 2, 2, 2]

    """
    return {letter : list(sum(1 for j in bwt[:i + 1] if j == letter)
                          for i in range(len(bwt)))
            for letter in set(bwt)}


def get_count(s):
    """
    Returns count information for the letters in the string 's'.

    count[letter] contains the number of symbols in 's' that are
    lexographically smaller than 'letter'.

    Examples:
    ---------

    >>> get_count('sassy') == {'a': 0, 'y': 4, 's': 1}
    True

    """
    return {letter : sum(1 for i in s if i < letter)
            for letter in set(s)}


def get_sa(s):
    r"""
    Returns the suffix array of 's'

    Examples:
    ---------

    >>> get_sa('banana\0')
    [6, 5, 3, 1, 0, 4, 2]

    """
    suffixes = list(s[i:] for i in range(len(s)))
    return list(suffixes.index(suffix) for suffix in sorted(suffixes))


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


def update_range(begin, end, letter, occ, count, length):
    """update (start, end) given a new letter"""
    newbegin = count[letter] + use_occ(occ, letter, begin - 1, length) + 1
    newend = count[letter] + use_occ(occ, letter, end, length)
    return newbegin, newend


def get_bwt_data(reference, eos=EOS):
    """Returns the data structures needed to perform BWT searches"""
    alphabet = set(reference)
    assert eos not in reference
    reference += eos
    bwt = get_bwt(reference)
    occ = get_occ(bwt)
    count = get_count(reference[:-1])
    sa = get_sa(reference)
    return alphabet, bwt, occ, count, sa


def bwt_match(query, reference, mismatches=0, bwt_data=None):
    """
    Find all matches of the string 'query' in the string 'reference', with at most
    'mismatch' mismatches

    Examples:
    ---------

    >>> bwt_match('abc', 'abcabcabc')
    [0, 3, 6]

    >>> bwt_match('gef', 'abcabcabc')
    []

    >>> bwt_match('abc', 'abcabd', mismatches=1)
    [0, 3]

    >>> bwt_match('abdd', 'abcabd', mismatches=1)
    []

    """
    assert len(query) > 0

    if bwt_data is None:
        bwt_data = get_bwt_data(reference)
    alphabet, bwt, occ, count, sa = bwt_data
    assert len(alphabet) > 0

    if not set(query) <= alphabet:
        return []

    length = len(bwt)
    results = []

    # a stack of partial matches
    class Partial(object):
        def __init__(self, **kwargs):
            self.__dict__.update(kwargs)

    partials = [Partial(query=query, begin=0, end=len(bwt) - 1,
                        mismatches=mismatches)]

    while len(partials) > 0:
        p = partials.pop()
        query = p.query[:-1]
        curletter = p.query[-1:]
        letters = [curletter] if p.mismatches == 0 else alphabet
        for letter in letters:
            mm = p.mismatches if letter == curletter else max(0, p.mismatches - 1)
            begin, end = update_range(p.begin, p.end, letter, occ, count, length)
            if begin <= end:
                if len(query) == 0:
                    results.extend(sa[begin : end + 1])
                else:
                    partials.append(Partial(query=query, begin=begin,
                                            end=end, mismatches=mm))
    return sorted(set(results))
