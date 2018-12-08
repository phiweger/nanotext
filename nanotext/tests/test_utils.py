import pkg_resources

from pybedtools import BedTool
import pytest

from nanotext.utils import get_interval_uid
from nanotext.utils import remove_overlap, deduplicate, create_domain_sequence


@pytest.fixture
def mock_domains():
    '''
    https://docs.pytest.org/en/latest/fixture.html
    stackoverflow.com/questions/779495
    '''
    fp = pkg_resources.resource_filename('nanotext', 'data/') + 'domains.bed'
    intervals = {}
    for i in BedTool(fp):
        intervals[get_interval_uid(i)] = i
    return intervals


def test_remove_overlap(mock_domains):
    text = []
    for i in remove_overlap(mock_domains):
        text.append(i.fields[3])
    assert text == [
        'zebra', 'zebra', 'monkey', 'monkey', 'monkey', 'giraffe', 'giraffe']


def test_deduplicate(mock_domains):
    text = []
    for i in deduplicate(remove_overlap(mock_domains)):
        text.append(i.fields[3])
    assert text == [
        'zebra', 'monkey', 'monkey', 'giraffe']


def test_create_domain_sequence_unknown(mock_domains):
    seq = create_domain_sequence(
        deduplicate(remove_overlap(mock_domains)), keep_unknown=True)
    assert seq['A'] == ['unknown', 'unknown', 'zebra']
    assert seq['B'] == ['monkey', 'unknown', 'unknown', 'monkey', 'giraffe']


def test_create_domain_sequence(mock_domains):
    seq = create_domain_sequence(
        deduplicate(remove_overlap(mock_domains)), keep_unknown=False)
    assert seq['A'] == ['zebra']
    assert seq['B'] == ['monkey', 'monkey', 'giraffe']
