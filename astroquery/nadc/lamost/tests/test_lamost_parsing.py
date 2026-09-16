# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Regression cases from historical LAMOST SQL responses and MRS products."""

from pathlib import Path

from astropy import units as u
from astropy.io import fits
import numpy as np
import pytest

from astroquery.exceptions import TableParseError
from ..core import LamostClass, parse_lrs_spectrum, parse_mrs_spectrum
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


@pytest.mark.parametrize('filename', ['mrs_dr8_excerpt.fits.gz', 'mrs_dr10_excerpt.fits.gz'])
def test_real_mrs_layouts_preserve_all_extensions(filename):
    path = DATA / filename
    result = parse_mrs_spectrum(path)
    with fits.open(path) as hdus:
        assert list(result) == [hdu.name for hdu in hdus[1:]]
        for hdu in hdus[1:]:
            if 'LOGLAM' in hdu.columns.names:
                expected_wave = np.power(10., np.asarray(hdu.data['LOGLAM'], dtype=float))
                expected_flux = hdu.data['FLUX']
            else:
                expected_wave = hdu.data['WAVELENGTH'][0]
                expected_flux = hdu.data['FLUX'][0]
            np.testing.assert_array_equal(result[hdu.name]['wavelength'], expected_wave)
            np.testing.assert_array_equal(result[hdu.name]['flux'], expected_flux)


@pytest.mark.parametrize('layout, diagnostic', [
    ('no_extensions', 'at least one spectrum extension'),
    ('image', 'nonempty binary table'), ('empty', 'nonempty binary table'),
    ('missing_flux', 'requires FLUX and WAVELENGTH'), ('missing_wave', 'requires FLUX and WAVELENGTH'),
    ('vector_loglam', 'nonempty numeric arrays of equal length'),
    ('scalar_wavelength', 'one table row of vector arrays'),
    ('mismatched', 'nonempty numeric arrays of equal length'),
    ('nonnumeric', 'nonempty numeric arrays of equal length'),
    ('ambiguous', 'ambiguous WAVELENGTH and LOGLAM'), ('duplicate', 'duplicate spectrum extension name'),
])
def test_mrs_rejects_unsupported_layouts(tmp_path, layout, diagnostic):
    columns = [fits.Column(name='FLUX', format='E', array=[1., 2.]),
               fits.Column(name='LOGLAM', format='D', array=[3.7, 3.8])]
    if layout == 'missing_flux':
        columns = columns[1:]
    elif layout == 'missing_wave':
        columns = columns[:1]
    elif layout == 'ambiguous':
        columns.append(fits.Column(name='WAVELENGTH', format='D', array=[5000., 6000.]))
    elif layout == 'scalar_wavelength':
        columns[1] = fits.Column(name='WAVELENGTH', format='D', array=[5000., 6000.])
    elif layout in ('vector_loglam', 'mismatched'):
        columns[0] = fits.Column(name='FLUX', format='2E', array=[[1., 2.]])
        columns[1] = fits.Column(name='LOGLAM', format='2D', array=[[3.7, 3.8]])
        if layout == 'mismatched':
            columns[1] = fits.Column(name='WAVELENGTH', format='3D', array=[[5000., 6000., 7000.]])
    elif layout == 'nonnumeric':
        columns[0] = fits.Column(name='FLUX', format='A', array=['a', 'b'])
    table = fits.BinTableHDU.from_columns(columns, name='B-123')
    if layout == 'empty':
        table = fits.BinTableHDU(data=table.data[:0])
    hdus = [fits.PrimaryHDU(), table]
    if layout == 'no_extensions':
        hdus = hdus[:1]
    elif layout == 'image':
        hdus[1] = fits.ImageHDU(data=np.ones(2))
    elif layout == 'duplicate':
        hdus.append(table.copy())
    path = tmp_path / 'invalid.fits'
    with fits.HDUList(hdus) as hdul:
        hdul.writeto(path)
    with pytest.raises(ValueError, match=diagnostic):
        parse_mrs_spectrum(path)


def test_historical_layout_is_mrs_only(tmp_path):
    path = tmp_path / 'historical.fits'
    with fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU.from_columns([
        fits.Column(name='FLUX', format='E', array=[2., 1.]),
        fits.Column(name='LOGLAM', format='E', array=[3.8, 3.7]),
    ])]) as hdus:
        hdus.writeto(path)
    result = parse_mrs_spectrum(path)['Extension_1']
    np.testing.assert_array_equal(result['flux'], [2., 1.])
    assert result['wavelength'][0] > result['wavelength'][1]
    assert result['wavelength'].dtype == np.float64
    with pytest.raises(ValueError):
        parse_lrs_spectrum(path)
