# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

from .. import conf
from ..core import LamostClass


pytestmark = pytest.mark.remote_data


@pytest.fixture
def lamost(tmp_path):
    with conf.set_temp('server', 'https://www.lamost.org/openapi'):
        client = LamostClass(token='', data_release='dr10', sub_version='v2.0')
    client.cache_location = tmp_path
    return client


def test_query_region_remote(lamost):
    center = SkyCoord(10.0004738, 40.9952444, unit='deg', frame='icrs')
    result = lamost.query_region(center, 0.2 * u.deg, cache=False)
    independent = lamost.query_sql(
        'SELECT obsid, ra, dec FROM catalogue '
        "WHERE spos(ra,dec) @ scircle '<(10.0004738d,40.9952444d),0.2d>' LIMIT 1000",
        output_format='csv', cache=False,
    )

    assert isinstance(result, Table)
    assert len(result) == 195
    assert 176604010 in result['obsid']
    assert len(independent) < 1000
    assert set(result['obsid']) == set(independent['obsid'])
    positions = SkyCoord(result['ra'], result['dec'], unit='deg')
    assert (center.separation(positions) <= 0.2 * u.deg).all()


@pytest.mark.parametrize('obsid, expected_rows', [(176604010, 1), (0, 0)])
def test_query_sql_remote(lamost, obsid, expected_rows):
    result = lamost.query_sql(
        f'SELECT obsid, ra, dec, teff, feh FROM combined WHERE obsid = {obsid} LIMIT 2',
        column_schema={
            'obsid': {'datatype': 'long'},
            'ra': {'datatype': 'double', 'unit': 'deg'},
            'dec': {'datatype': 'double', 'unit': 'deg'},
            'teff': {'datatype': 'float', 'unit': 'K'},
            'feh': {'datatype': 'float'},
        },
        cache=False,
    )

    assert isinstance(result, Table)
    assert len(result) == expected_rows
    assert result['obsid'].dtype.kind in 'iu'
    assert result['teff'].dtype.kind == 'f'
    assert result['teff'].unit == u.K
    if expected_rows:
        assert result['obsid'][0] == obsid
        assert result['ra'][0] == pytest.approx(10.008848, abs=1e-6)
        assert result['dec'][0] == pytest.approx(40.969976, abs=1e-6)


def test_query_ssap_remote(lamost):
    center = SkyCoord(10.008848, 40.969976, unit='deg', frame='icrs')
    result = lamost.query_ssap(center, radius='5 arcsec', cache=False)

    assert isinstance(result, Table)
    assert len(result) >= 1
    assert 176604010 in result['obsid']
    positions = SkyCoord(result['ra'], result['dec'], unit='deg')
    assert (center.separation(positions) <= 5 * u.arcsec).all()


@pytest.mark.parametrize('obsid, expected_rows', [(176604010, 1), (0, 0)])
def test_query_catalog_remote(lamost, obsid, expected_rows):
    result = lamost.query_catalog(
        'combined', columns=['obsid', 'ra', 'dec'],
        column_constraints=[{'column_name': 'obsid', 'operation': 'equal', 'constraint': str(obsid)}],
        max_rows=2, cache=False,
    )

    assert isinstance(result, Table)
    assert len(result) == expected_rows
    assert result.meta['catalog'] == 'combined'
    assert result['obsid'].dtype.kind in 'iu'
    if expected_rows:
        assert result['obsid'][0] == obsid
        assert result['ra'][0] == pytest.approx(10.008848, abs=1e-6)
        assert result['dec'][0] == pytest.approx(40.969976, abs=1e-6)


def test_query_spectra_nearest_remote(lamost):
    center = SkyCoord(10.0004738, 40.9952444, unit='deg', frame='icrs')
    result = lamost.query_spectra(
        center, 0.2 * u.deg, nearest_only=True,
        columns=['obsid', 'ra', 'dec'], cache=False,
    )

    assert len(result) == 1
    assert result['obsid'][0] == 176604010
    assert result['ra'][0] == pytest.approx(10.008848, abs=1e-6)
    assert result['dec'][0] == pytest.approx(40.969976, abs=1e-6)


def test_get_dr_versions_remote(lamost):
    result = lamost.get_dr_versions()

    assert isinstance(result, list)
    assert len(result) > 0
    assert {"dr_version", "sub_version", "public_status"}.issubset(result[0])
    assert any(version['dr_version'] == 'dr10' and version['sub_version'] == 'v2.0' for version in result)


def test_get_metadata_and_stellar_parameters_remote(lamost):
    metadata = lamost.get_metadata(176604010, cache=False)
    assert len(metadata) == 1
    assert int(metadata['obsid'][0]) == 176604010
    result = lamost.query_stellar_parameters(
        SkyCoord(10.008848, 40.969976, unit='deg'), '5 arcsec', nearest_only=True, cache=False,
    )
    assert list(result['obsid']) == [176604010]
    # Numeric conversion is required; units are attached only if the service
    # declares them (the reference combined schema currently omits this unit).
    assert result['teff'].dtype.kind == 'f'
    assert float(result['teff'][0]) == pytest.approx(float(metadata['teff'][0]), abs=.01)


def test_related_observations_remote(lamost):
    result = lamost.query_repeat_observations(obsid=176604010, cache=False)
    assert result['unique_id']
    assert 176604010 in list(map(int, result['related_obsids_low']))
    assert len(set(result['related_obsids_low'])) == len(result['related_obsids_low'])
