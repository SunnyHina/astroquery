# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Invalid MRS wavelengths and per-file batch outcomes."""

import gzip
from pathlib import Path

import numpy as np
import pytest

from .. import parse_mrs_spectrum, parse_mrs_spectra
from .. import core
from .helpers import DATA, write_spectrum


@pytest.mark.parametrize('loglam', [np.nan, np.inf, -np.inf, -400., 400.],
                         ids=['nan', 'inf', 'minus_inf', 'underflow', 'overflow'])
def test_mrs_invalid_loglam_has_file_extension_and_pixel_diagnostics(tmp_path, loglam):
    path = write_spectrum(tmp_path / 'bad-loglam.fits', [3.7, loglam, 3.8], [1., 2., 3.], loglam=True)
    # Conversion failures must remain ValueError even with strict NumPy flags.
    with np.errstate(all='raise'), pytest.raises(ValueError) as caught:
        parse_mrs_spectrum(path)
    message = str(caught.value)
    assert str(path) in message
    assert "extension 'COADD_B' (HDU 1)" in message
    assert '1 invalid pixels' in message and 'zero-based index 1' in message


@pytest.mark.parametrize('wavelength', [np.nan, np.inf, -np.inf, 0., -1.])
def test_mrs_direct_wavelengths_must_be_finite_and_positive(tmp_path, wavelength):
    path = write_spectrum(tmp_path / 'bad-wave.fits', [5000., wavelength, 6000.], [1., 2., 3.])
    with pytest.raises(ValueError, match='finite and positive') as caught:
        parse_mrs_spectrum(path)
    assert 'zero-based index 1' in str(caught.value)


@pytest.mark.parametrize('loglam', [False, True])
def test_mrs_valid_wavelength_order_and_raw_flux_are_preserved(tmp_path, loglam):
    wavelength = np.array([.1, 10., 1., 5000.])
    flux = np.array([np.nan, np.inf, -1., 0.])
    path = write_spectrum(tmp_path / 'raw-flux.fits', np.log10(wavelength) if loglam else wavelength,
                          flux, loglam=loglam)
    result = parse_mrs_spectrum(path)['COADD_B']
    np.testing.assert_allclose(result['wavelength'], wavelength, rtol=1e-15)
    np.testing.assert_array_equal(result['flux'], flux)


def test_mrs_batch_records_failures_and_continues_in_input_order(tmp_path):
    bad = write_spectrum(tmp_path / 'bad.fits', [np.nan, 5000.], [1., 2.])
    good = DATA / 'mrs_dr8_excerpt.fits.gz'
    missing = tmp_path / 'missing.fits'
    filenames = [bad, good, missing, good]
    original = {path: path.read_bytes() for path in [bad, good]}
    spectra, manifest = parse_mrs_spectra(iter(filenames))
    assert manifest.colnames == ['Local Path', 'Status', 'Message']
    assert list(manifest['Local Path']) == [str(path) for path in filenames]
    assert list(manifest['Status']) == ['ERROR', 'COMPLETE', 'ERROR', 'COMPLETE']
    assert spectra[0] is None and spectra[2] is None
    assert 'ValueError:' in manifest['Message'][0] and 'COADD_B' in manifest['Message'][0]
    assert manifest['Message'][1] == manifest['Message'][3] == ''
    assert 'FileNotFoundError:' in manifest['Message'][2]
    expected = parse_mrs_spectrum(good)
    for result in [spectra[1], spectra[3]]:
        assert list(result) == list(expected)
        for extension in expected:
            np.testing.assert_array_equal(result[extension]['wavelength'], expected[extension]['wavelength'])
            np.testing.assert_array_equal(result[extension]['flux'], expected[extension]['flux'])
    for path, content in original.items():
        assert path.read_bytes() == content


@pytest.mark.parametrize('count', [0, 2])
def test_mrs_batch_empty_or_all_failed(tmp_path, count):
    spectra, manifest = parse_mrs_spectra([tmp_path / 'missing.fits'] * count)
    assert spectra == [None] * count
    assert len(manifest) == count
    assert manifest.colnames == ['Local Path', 'Status', 'Message']
    assert all(manifest['Status'] == 'ERROR')


def test_mrs_batch_records_truncated_compressed_file(tmp_path):
    good = write_spectrum(tmp_path / 'good.fits', [5000., 6000.], [1., 2.])
    truncated = tmp_path / 'truncated.fits.gz'
    compressed = gzip.compress(good.read_bytes())
    # Remove compressed FITS content, not only the unused gzip trailer.
    truncated.write_bytes(compressed[:len(compressed) // 2])
    spectra, manifest = parse_mrs_spectra([truncated, good])
    assert list(manifest['Status']) == ['ERROR', 'COMPLETE']
    assert spectra[0] is None and spectra[1] is not None


@pytest.mark.parametrize('filename', ['spectrum.fits', Path('spectrum.fits')])
def test_mrs_batch_rejects_single_path(filename):
    with pytest.raises(TypeError, match='iterable of paths'):
        parse_mrs_spectra(filename)


def test_mrs_batch_does_not_hide_unexpected_errors(monkeypatch):
    def unexpected_error(filename):
        raise RuntimeError('Unexpected parser failure')

    monkeypatch.setattr(core, 'parse_mrs_spectrum', unexpected_error)
    with pytest.raises(RuntimeError, match='Unexpected parser failure'):
        parse_mrs_spectra(['spectrum.fits'])
