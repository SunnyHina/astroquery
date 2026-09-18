# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Regression cases from historical LAMOST table responses."""

from pathlib import Path

from astropy import units as u
import pytest

from astroquery.exceptions import TableParseError
from ..core import LamostClass
from .helpers import create_mock_response


DATA = Path(__file__).parent / 'data'


@pytest.mark.parametrize('body', [b'', b' \t\r\n', b'\xef\xbb\xbf \n'])
@pytest.mark.parametrize('content_type', ['text/csv', 'text/plain', 'application/json', 'application/x-votable+xml'])
def test_blank_response_is_not_a_zero_row_table(body, content_type):
    client = LamostClass(token='')
    with pytest.raises(TableParseError, match='empty response body'):
        client._parse_result(create_mock_response(content=body, content_type=content_type))
    assert client.response.content == body


@pytest.mark.parametrize('body,names,values', [
    ('name,value\nA|B,2\n', ['name', 'value'], ['A|B', 2]),
    ('name|value\nA,B|2\n', ['name', 'value'], ['A,B', 2]),
    ('"name|unit",value\n"A|B",2\n', ['name|unit', 'value'], ['A|B', 2]),
    ('"name,unit"|value\n"A,B"|2\n', ['name,unit', 'value'], ['A,B', 2]),
    ('"name|unit"\n"A|B"\n', ['name|unit'], ['A|B']),
    ('name\nA|B\n', ['name'], ['A|B']),
])
def test_delimiter_evidence_comes_from_header(body, names, values):
    result = LamostClass(token='')._parse_csv_result(create_mock_response(content=body.encode()))
    assert result.colnames == names
    assert list(result[0]) == values


@pytest.mark.parametrize('body, diagnostic', [
    ('id,ra|dec\n1,2|3\n', 'Ambiguous table header'),
    ('id|ra|dec\n1|2\n', 'has 2 fields; expected 3'),
    ('id,ra,dec\n1,2\n', 'has 2 fields; expected 3'),
    ('id,ra\n1,2,3\n', 'has 3 fields; expected 2'),
    ('id|id\n1|2\n', 'empty or duplicate column names'),
    ('id,,dec\n1,2,3\n', 'empty or duplicate column names'),
    ('id|ra\n1|"unterminated\n', 'unexpected end of data'),
    ('id,ra\n1,"quoted"junk\n', 'expected after'),
])
def test_damaged_or_ambiguous_csv_is_rejected(body, diagnostic):
    with pytest.raises(TableParseError, match=diagnostic):
        LamostClass(token='')._parse_csv_result(create_mock_response(content=body.encode()))


@pytest.mark.parametrize('name,rows,columns', [('med_catalogue', 10, 6), ('med_stellar', 1, 80)])
def test_preserved_historical_sql_response(name, rows, columns):
    response = create_mock_response(content=(DATA / f'{name}_dr8.txt').read_bytes(), content_type='text/csv')
    schema = {'obsid': {'datatype': 'long'}, 'rv_br0': {'datatype': 'double', 'unit': 'km/s'}}
    if name == 'med_stellar':
        schema['gaia_source_id'] = {'datatype': 'char'}
    result = LamostClass(token='')._parse_result(response, column_schema=schema)
    assert len(result) == rows
    assert len(result.colnames) == columns
    assert set(result['obsid']) == {635003103}
    assert result['rv_br0'].unit == u.km / u.s
    assert result['rv_br0'][0] == (-3.97 if name == 'med_catalogue' else -1.08)
    if name == 'med_stellar':
        assert result['gaia_source_id'][0] == '3700975728440669184'
        assert result['mobsid'][0] == '635003103R'
