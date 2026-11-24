from main import indexed_match, join_units

def test_indexed_match_finds_two_runs():
    assert indexed_match("(?:GAA)+", "GAAGAATTGAA") == [(0, 6), (8, 11)]

def test_tract_kept_when_units_equal_min_units():
    tracts = join_units([(0, 9)], 9, 3, "GAAGAAGAA", False)
    assert len(tracts) == 1
    assert tracts[0]["units"] == 3
    