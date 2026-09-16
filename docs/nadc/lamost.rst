NADC LAMOST Queries
===================

``astroquery.nadc.lamost`` provides access to the LAMOST archive for catalog
queries, metadata lookups, and spectrum-related utilities.

Configuration
-------------

The base URL, timeout in seconds, default data release, sub-version, and token
are read when a `~astroquery.nadc.lamost.LamostClass` instance is created.
The imported ``Lamost`` object is an instance created at import time. Changing
``conf`` does not update existing instances; create a new instance after
changing configuration:

.. doctest::

  >>> from astroquery.nadc.lamost import LamostClass, conf
  >>> with conf.set_temp('timeout', 120):
  ...     lamost = LamostClass(token='', data_release='dr10', sub_version='v2.0')
  >>> lamost.TIMEOUT
  120

``LamostClass`` resolves tokens in this order: an explicit ``token`` argument,
``conf.token``, environment variables such as
``ASTROQUERY_NADC_LAMOST_TOKEN``, and an explicitly requested pylamost-style
configuration file. Pass ``token=''`` to force anonymous access, including
when a token is configured elsewhere. An empty ``conf.token`` permits the
environment/config-file fallback; it does not force anonymous access.

Configure a token using ``astroquery.cfg``, an environment variable, or
``conf`` before constructing the client. For example:

.. doctest::

  >>> conf.token = 'your-token'  # doctest: +SKIP
  >>> authenticated = LamostClass()  # doctest: +SKIP
  >>> configured = LamostClass(pylamost_config='~/pylamost.ini')  # doctest: +SKIP

The last example reads ``token=your-token`` from the named file only if no
higher-priority source supplied a token. The client does not search for this
file automatically. Authenticated requests disable response caching.
``conf.server`` must include the OpenAPI base path, for example
``https://www.lamost.org/openapi``. The timeout applies to queries and data
downloads.

Basic Usage
-----------

``query_region`` and ``query_ssap`` request CSV by default. Some archive VOTable
responses declare string columns too short for their data; the client raises
`~astroquery.exceptions.TableParseError` instead of returning truncated
identifiers when VOTable is explicitly requested.

.. doctest::

  >>> import astropy.units as u
  >>> from astropy.coordinates import SkyCoord
  >>> coord = SkyCoord(10.0004738, 40.9952444, unit='deg', frame='icrs')
  >>> payload = lamost.query_region(
  ...     coord, radius=0.2*u.deg, output_format='csv', get_query_payload=True)
  >>> payload['ra'], payload['dec'], payload['output.fmt']
  (10.0004738, 40.9952444, 'csv')
  >>> matches = lamost.query_region(  # doctest: +SKIP
  ...     coord, radius=0.2*u.deg, output_format='csv')

Catalog query methods return `~astropy.table.Table` objects. Coordinates are
transformed to ICRS. Radii accept angular quantities such as ``5*u.arcsec``
and angle strings such as ``'5 arcsec'``. Bare numbers mean degrees for
``query_region``, ``query_ssap``, and ``query_repeat_observations``; they mean
arcseconds for the structured ``query_spectra`` and
``query_stellar_parameters`` methods. Use explicit units to avoid ambiguity.
Invalid or non-angular radii raise `~astroquery.exceptions.InvalidQueryError`.

CSV avoids the known VOTable format problem; it does not establish that the
service returned every match for any release or search size. A single query
does not automatically retrieve additional pages.

SQL-style queries and structured catalog requests are also available:

.. doctest::

  >>> sql_payload = lamost.query_sql('SELECT * FROM combined LIMIT 5', get_query_payload=True)
  >>> sql_payload['output.fmt']
  'json'
  >>> catalog_payload = lamost.query_catalog(
  ...     'combined',
  ...     columns=['obsid', 'ra', 'dec'],
  ...     max_rows=5,
  ...     get_query_payload=True,
  ... )
  >>> catalog_payload['rows']
  5

``get_query_payload=True`` returns a request mapping with token values
redacted, without submitting the LAMOST query or fetching schema metadata.
For authenticated structured catalog queries, the mapping contains ``json``
and ``params`` entries for the body and query parameters. Pass explicit
coordinates, as above, to avoid online name resolution of an object name.

For ``query_catalog`` and its wrappers, ``max_rows`` limits a single page
(default: 100), and ``page`` selects the one-based page number. These methods
do not aggregate pages.

Data Release and Metadata
-------------------------

Use ``get_dr_versions`` to inspect available data-release and sub-version
combinations. The instance's ``data_release`` and ``sub_version`` select the
archive endpoint used by query and data-product methods.
``get_tables_metadata`` returns schema metadata, ``get_tap_url`` returns a
dictionary of TAP connection information, and ``get_footprint`` returns
image bytes for the selected resolution.

Release and Format Boundaries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The reference configuration is the public DR10/v2.0 service. A release listed
by ``get_dr_versions`` does not establish that every endpoint or output format
works for that release. Select both version components explicitly when
reproducing an observation; do not replace a historical release with a newer
one merely because its endpoint responds.

The following boundaries were checked on 2026-09-14:

.. list-table:: Representative coverage
   :header-rows: 1
   :widths: 22 35 43

   * - Release
     - Operation
     - Scope
   * - DR10/v2.0
     - SQL, cone search, SSAP, structured queries
     - JSON/CSV are the reference table transports. Explicit SQL TXT can
       return a service error; VOTable can contain invalid string declarations.
   * - DR7/v2.0
     - Metadata and LRS spectrum
     - Metadata is available for observation 54901214.
   * - DR8/v1.0, DR9
     - Historical SQL and local MRS parsing
     - Legacy pipe-delimited SQL text and scalar FLUX/LOGLAM FITS tables
       have offline regression samples. DR9/v0 SQL can return an empty body.
   * - DR6
     - Table discovery
     - Historical table endpoints can return non-table or service-error
       responses. This does not establish that DR6 spectrum downloads fail.

Format support in the parser is distinct from service availability. Supported
table inputs are JSON, CSV, tab-separated TXT, and valid VOTable. Historical
pipe-delimited text labelled as CSV is recognized from unquoted delimiters in
its header. Quoted commas, pipes, whitespace, and embedded newlines remain
field content, including blank lines and the original CR/LF line endings.
Mixed unquoted header delimiters, missing fields, duplicate or empty column
names, and broken quoting raise `~astroquery.exceptions.TableParseError`.
When a column name itself contains a delimiter, the service must quote it.

Executed ``query_catalog`` requests validate output, constraint, position, and
sort columns against ``get_tables_metadata`` before submitting the query.
Use ``cache=False`` to refresh cached metadata and query responses. The
returned table must also contain every requested output column.
Structured results record their source as ``table.meta['catalog']``. Catalog
names are specific to the LAMOST release; discover them from this metadata.

.. doctest::

  >>> versions = lamost.get_dr_versions()  # doctest: +SKIP
  >>> versions[0]["dr_version"]  # doctest: +SKIP
  'dr10'
  >>> metadata = lamost.get_tables_metadata()  # doctest: +SKIP
  >>> "tables" in metadata  # doctest: +SKIP
  True

Column Types and Units
----------------------

``query_region``, ``query_catalog``, and its spectral-query wrappers use catalog metadata to
convert numeric columns and attach recognized units. For example, ``obsid``
is an integer when declared ``long``; character identifiers such as
``gaia_source_id`` remain strings, preserving leading zeros. Missing numeric
values are masked. Available schema information also determines the column
types of empty results.

``query_sql`` uses column definitions included in the response. Raw JSON
without these definitions preserves the service's types: a value such as
``teff='5770'`` remains a string. The client does not infer types from SQL
expressions or aliases. Supply ``column_schema`` using the actual result
column names when known types are required:

.. doctest::

  >>> column_schema = {
  ...     'obsid': {'datatype': 'long'},
  ...     'temperature': {'datatype': 'double', 'unit': 'K'},
  ...     'feh': {'datatype': 'double'},
  ...     'gaia_source_id': {'datatype': 'char'},
  ... }
  >>> stars = lamost.query_sql(  # doctest: +SKIP
  ...     'SELECT obsid, teff AS temperature, feh, gaia_source_id '
  ...     'FROM combined LIMIT 5', column_schema=column_schema)
  >>> hot_stars = stars[stars['temperature'] > 5500]  # doctest: +SKIP

The comparison above uses the column's numeric values in kelvin. The schema
must describe the returned expression, including any SQL unit conversion.
With ``column_schema``, CSV/TXT fields are read as strings before conversion;
columns without a declared datatype remain strings. This preserves character
identifiers such as ``'00123'``.
CSV character fields retain leading/trailing whitespace and embedded line
endings. Under a character schema, whitespace-only strings remain values;
empty fields are masked. Whitespace-only numeric fields remain masked.
Values that cannot be converted to their declared datatype raise
`~astroquery.exceptions.TableParseError`. Floating-point NaN and infinity
are preserved; conversion does not establish scientific validity. Check
column masks and ``numpy.isfinite`` when selecting numeric data for analysis.
Temperature cuts use kelvin, ``logg`` cuts use the base-10 logarithm of surface
gravity in cm/s\ :sup:`2`, and ``feh`` cuts use [Fe/H] in dex. S/N thresholds
are dimensionless. Units are attached only when declared and recognized;
check ``table[column].unit`` before combining results from different sources.

Spectral Sample Queries
-----------------------

Use ``query_spectra`` when the task is to build
a LAMOST spectral sample with common quality cuts.  It translates query
parameters such as SNR and stellar-parameter ranges into the structured
``query_catalog`` payload.

The current structured cone-query service may return only the nearest match
even when ``nearest_only=False``. Executing such a query raises
`~astroquery.exceptions.RemoteServiceError`; a one-row response cannot be
treated as all matches. This applies to ``query_spectra``,
``query_stellar_parameters``, and direct ``query_catalog`` cone constraints
without ``cone_nearestonly=True``. Payload inspection still works. To search a
field, use ``query_region(..., output_format='csv')`` as shown above. Use
``nearest_only=True`` only when the nearest match is the intended result.
Structured queries without a spatial constraint remain available.

.. doctest::

  >>> spectral_payload = lamost.query_spectra(
  ...     coord,
  ...     5*u.arcsec,
  ...     snr_min=20,
  ...     teff_range=(4500, 6500),
  ...     logg_range=(3.5, 5.0),
  ...     feh_range=(-1.0, 0.5),
  ...     nearest_only=True,
  ...     columns=["obsid", "ra", "dec", "snrg", "teff", "logg", "feh"],
  ...     get_query_payload=True,
  ... )
  >>> spectral_payload["column_constraints"][0]["column_name"]
  'snrg'

Use ``query_stellar_parameters`` when the
desired output is focused on stellar atmospheric parameters.  The default
columns are ``obsid``, ``ra``, ``dec``, ``teff``, ``logg``, ``feh``, and the
selected SNR column for LRS. The default LRS S/N field is ``snrg``; callers can
select ``snru``, ``snrr``, ``snri``, or ``snrz`` explicitly with
``snr_column``. MRS uses ``snr``, ``teff_lasp``, ``logg_lasp``, and
``feh_lasp``.

.. doctest::

  >>> stellar_payload = lamost.query_stellar_parameters(
  ...     teff_range=(4500, 6500),
  ...     snr_min=30,
  ...     get_query_payload=True,
  ... )
  >>> stellar_payload["showcol"]
  ['obsid', 'ra', 'dec', 'teff', 'logg', 'feh', 'snrg']

Use ``query_repeat_observations`` to resolve one observation ID, or one
position and radius, to related observations. These two input forms are
mutually exclusive. Its return value is a dictionary containing
``unique_id``, ``related_obsids``, ``related_obsids_low``, and
``related_obsids_medium``; the latter lists distinguish LRS and MRS IDs.

.. doctest::

  >>> repeat_payload = lamost.query_repeat_observations(
  ...     coordinates=coord,
  ...     radius=3*u.arcsec,
  ...     get_query_payload=True,
  ... )
  >>> repeat_payload["ra"], repeat_payload["dec"]
  (10.0004738, 40.9952444)

Paged Query Results
-------------------

The four ``sqlid`` methods below are provisional compatibility interfaces
requiring an externally supplied server-side query ID. The
`LAMOST OpenAPI specification <https://www.lamost.org/openapi/openapi.yaml>`_
reviewed on 2026-09-14 documents result retrieval but no submission endpoint
producing such an ID. ``query_sql`` and ``query_catalog`` return parsed tables,
not IDs. These interfaces have offline regression tests; a complete live
count/page/export workflow has not been validated.

Use ``get_query_result_count`` and ``get_query_result_by_page`` for explicit
pagination of an externally supplied server-side ``sqlid``, or
``get_query_result`` to fetch all pages with a fixed ``page_size``.
Page numbers are one-based. Keep ``rows`` fixed across pages, including the
last page, because it determines the service's page offset. A JSON page or
aggregate returns a list of records; the other formats return an Astropy
table. ``download_query_result`` writes the result set to a local file and
returns its path.

``download_query_result`` supports ``output_format='csv'``, ``'json'``,
``'votable'``, and ``'txt'``. TXT exports are UTF-8 tab-separated files with a
column-name header, readable with
``Table.read(path, format='ascii.tab', encoding='utf-8')``.
Each export fetches the complete result once in the required transport
format. CSV, JSON, and TXT exports fetch JSON records to preserve character
values such as identifiers with leading zeros; VOTable exports retain
VOTable column metadata. Text readers can infer numeric types when reopening
a TXT file; pass ``converters={'gaia_source_id': str}`` to ``Table.read``
when that column must remain a string.
For a zero-row result, JSON exports ``[]``. CSV and TXT cannot reconstruct
column names from an empty JSON record list, so they contain no table header;
check the result count before reopening these files as tables.
The exporter writes to a temporary file on the destination filesystem and
replaces the destination only after writing succeeds. A query or write
failure preserves an existing destination and removes temporary output.

Page sizes and page numbers must be positive integers; counts must be
nonnegative integers. An invalid service count or a page length inconsistent
with that count raises `~astroquery.exceptions.RemoteServiceError` instead
of returning partial results. These checks cannot detect changed records
when the service preserves the same row counts.

.. doctest::

  >>> count = lamost.get_query_result_count(12345)  # doctest: +SKIP
  >>> page = lamost.get_query_result_by_page(  # doctest: +SKIP
  ...     12345, count, rows=1000, page=1, output_format="json"
  ... )
  >>> results = lamost.get_query_result(  # doctest: +SKIP
  ...     12345, output_format="json", page_size=1000
  ... )
  >>> path = lamost.download_query_result(  # doctest: +SKIP
  ...     12345, "lamost-results.csv", output_format="csv"
  ... )

Data Products
-------------

Data-product methods use an observation ID and ``resolution='low'`` for LRS
or ``resolution='medium'`` for MRS. ``get_metadata`` returns a table,
``get_spectra`` returns a list of `~astropy.io.fits.HDUList` objects,
``get_images`` returns a list of PNG byte strings, and ``get_fits_csv``
returns CSV text. ``get_spectrum_list`` and ``get_image_list`` return URL
lists without downloading files. Those URLs include the token when
authenticated access is configured; treat them as credentials and avoid
sharing or logging them. Their ``get_query_payload=True`` mode returns a
redacted parameter mapping.

``get_spectra`` verifies downloaded HDU lists with Astropy's ``"warn"`` mode
by default and accepts another FITS verification option through ``verify``.
The caller must close the HDU lists after use:

.. doctest::

  >>> spectra = lamost.get_spectra(176604010, resolution='low')  # doctest: +SKIP
  >>> try:  # doctest: +SKIP
  ...     obsid = spectra[0][0].header['OBSID']
  ... finally:
  ...     for spectrum in spectra:
  ...         spectrum.close()

``download_catalog`` downloads the named catalog product and returns its
local path. It uses ``verify='exception'`` by default and replaces the
destination only after verification with the selected option succeeds. With
``overwrite=False``, an existing file is returned without modification or
revalidation; the server is still contacted to resolve the filename. Failed
writes or verification preserve an existing destination and clean up the
unique temporary download. Some archive products fail strict verification
even though Astropy can open them; inspect those warnings before changing
``verify``. FITS verification checks file structure, not scientific pixel
quality. HTTP responses are closed after downloading or skipping a product.

Anonymous non-streaming requests use the inherited ``BaseQuery`` response
cache unless ``cache=False`` is available and supplied. ``get_spectra`` and
``get_images`` also use that response cache for anonymous requests; they do
not create a persistent spectrum library. Authenticated requests and streaming
catalog downloads bypass the response cache. ``download_catalog`` can skip an
existing destination but does not resume a partial download. No automatic
cross-release fallback is performed.

Local Spectrum Processing
-------------------------

``parse_lrs_spectrum`` reads the supported single-HDU image or two-HDU table
layout and returns wavelength, flux, and median-filtered flux arrays for
windows of 7 and 15 pixels, with zero padding at the edges. Unsupported
layouts raise ``ValueError``. LRS parsing does not reject nonfinite or
nonpositive wavelength/flux values; apply quality checks before analysis.
``parse_mrs_spectrum`` reads the spectral-band extensions and
returns a dictionary keyed by extension name, with ``wavelength`` and
``flux`` arrays in each entry. Wavelengths are in angstroms; fluxes retain
the archive file's units and normalization. Consult the FITS headers before
interpreting them as calibrated flux.

MRS tables support two explicit layouts: one row containing vector ``FLUX``
and ``WAVELENGTH`` columns, or successive rows containing scalar ``FLUX`` and
``LOGLAM`` pixels. Historical logarithmic wavelengths are converted using
``10**LOGLAM`` in double precision. Coadd names such as ``COADD_B`` and exposure
names such as ``B-83692603`` are retained. The parser does not sort pixels,
apply a radial-velocity correction, repair observation times, or change the
wavelength frame. Empty or unsupported extensions, mismatched arrays, ambiguous
wavelength columns, and duplicate extension names raise ``ValueError``.
``LOGLAM`` values must be finite, and wavelengths in either layout must be
finite and positive. Conversion overflow or underflow to zero is rejected.
Wavelength-value errors identify the file, extension, and first invalid pixel
(zero-based index).
Flux values are preserved; pixel quality selection belongs to the analysis.

For several local files, ``parse_mrs_spectra`` applies the same strict parser
and records each file's outcome:

.. code-block:: python

   from astroquery.nadc.lamost import parse_mrs_spectra

   spectra, manifest = parse_mrs_spectra(['first.fits', 'second.fits'])
   print(manifest['Local Path', 'Status', 'Message'])
   for spectrum, row in zip(spectra, manifest):
       if row['Status'] == 'COMPLETE':
           print(list(spectrum))

Both outputs preserve input order and repeated paths. A failed file has
``Status='ERROR'``, an exception type and reason in ``Message``, and ``None``
in the corresponding ``spectra`` entry. Successful files have
``Status='COMPLETE'`` and an empty message. File I/O, truncated-file, and
validation errors are recorded without stopping later files; unexpected
exceptions propagate. No input file or pixel is repaired or discarded.

Failures and Diagnostics
------------------------

Invalid query parameters or unknown catalog columns raise
`~astroquery.exceptions.InvalidQueryError`. Authentication failures raise
`~astroquery.exceptions.LoginError`; configure a valid token and create a
new client as described above. Other HTTP failures raise `requests.HTTPError`;
error payloads returned with HTTP success raise
`~astroquery.exceptions.RemoteServiceError`. Both preserve available, redacted
error details. A missing endpoint or unsupported release is a service/address
problem and does not by itself establish that a token is required.

Malformed responses, missing requested columns, failed datatype conversions,
and VOTables that would truncate data raise
`~astroquery.exceptions.TableParseError`. Query and JSON metadata endpoints
also reject empty bodies and HTML pages, including mislabeled HTML. A text
response with column headers and zero data rows is valid, as are supported
empty JSON and VOTable results.

When response parsing fails, ``client.response`` retains a redacted diagnostic
response. HTML and empty-body errors include the HTTP status, redacted URL,
content type, and body length. Payload inspection does not test authentication,
schema validity, or service availability.

Exceptions
----------

.. autoexception:: astroquery.exceptions.InvalidQueryError

.. autoexception:: astroquery.exceptions.LoginError

.. autoexception:: astroquery.exceptions.RemoteServiceError

.. autoexception:: astroquery.exceptions.TableParseError

Reference/API
=============

.. automodapi:: astroquery.nadc.lamost
    :no-inheritance-diagram:
