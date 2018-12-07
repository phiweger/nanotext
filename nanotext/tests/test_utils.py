import pkg_resources

from pybedtools import BedTool
import pytest

from nanotext.utils import get_interval_uid, remove_overlap, deduplicate


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
    assert text == ['zebra', 'zebra', 'monkey', 'elefant']


def test_deduplicate(mock_domains):
    text = []
    for i in deduplicate(remove_overlap(mock_domains)):
        text.append(i.fields[3])
    assert text == ['zebra', 'monkey', 'elefant']


