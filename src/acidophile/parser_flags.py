"""The detector for the parser's one real failure mode, shared by the scripts that need it.

The rule-based parser behind 92% of the phenotypes is exact for 98.8% of minima against the hand-curated rows, and its one
reproducible failure is a sentence that reports only an optimum - "grew optimally at 35 C, at pH 7.0-8.0" - whose optimum is
recorded as the growth range, biasing the minimum upward. `flag` detects exactly that, on four conditions together:
an "optim" stem within LOOKBACK characters before the parser's own range match; a matched interval no wider than MAXWIDTH pH
units, because an "optimum" spanning four units is a growth range whatever the sentence calls it; no optimum within TRAILING
characters after the match, which would signal the legitimate "range (optimum ...)" form; and no clause stating a growth range
independently of the optimum.
It was validated in scripts/61 against fifteen statements read by hand and recorded there verbatim with their verdicts, where it
scores 8 true positives, 0 false positives and 0 misses. It lives here rather than in that script so that scripts/64 and anything
later test the same detector rather than a copy of it; the hand-read ground truth stays in scripts/61, which is where it is run.
"""
import re

N = r"(\d+(?:\.\d+)?)"
SEP = r"\s*(?:–|-|—|−|to|and)\s*"
RANGE = re.compile(r"pH(?:\s*values?)?(?:\s*range)?(?:\s*(?:of|from|between|for growth|ranging from))?\s*\(?\s*" + N + SEP + N)
OPTSTEM = re.compile(r"optim", re.I)
# a clause that states a growth range independently of any optimum; if one of these is present the sentence is not optimum-only
RANGEWORD = re.compile(r"(?:growth (?:was )?(?:observed|occurred|detected|possible|range)|grew (?:at|between|over|from)|able to grow|occurred at|observed (?:at|between|over)|pH range|range of pH|ranging from|between pH)", re.I)
LOOKBACK = 60     # characters before the range match searched for an "optim" stem
MAXWIDTH = 2.0    # an "optimum" wider than this is a growth range, whatever the sentence calls it
TRAILING = 40     # characters after the match searched for a following optimum, which signals "range (optimum ...)"


def flag(statement: str) -> bool:
    """True if the parser's range match is an optimum clause misread as the growth range. See the module docstring for the four conditions."""
    if not isinstance(statement, str):
        return False
    m = RANGE.search(statement)
    if not m:
        return False
    if float(m.group(2)) - float(m.group(1)) > MAXWIDTH:
        return False
    if not OPTSTEM.search(statement[max(0, m.start() - LOOKBACK):m.start()]):
        return False
    if OPTSTEM.search(statement[m.end():m.end() + TRAILING]):
        return False
    return not RANGEWORD.search(statement)


def flagged_pmids(candidates, used_pmids):
    """The PMIDs whose statement trips the detector, among the candidate rows behind the parser tier."""
    inuse = candidates[candidates.pmid.isin(used_pmids) & candidates.usable_noopt]
    return set(inuse.pmid[inuse.statement.map(flag)]), inuse
