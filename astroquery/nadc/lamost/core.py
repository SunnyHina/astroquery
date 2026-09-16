# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
LAMOST Spectroscopic Survey Query Tool
=======================================

This module provides the core implementation for querying LAMOST data.
"""

# Standard library
from collections.abc import Mapping
from copy import copy
import csv
import os
import re
import tempfile
import warnings
from io import BytesIO, StringIO
from numbers import Integral
from urllib.parse import quote, quote_plus

# Third party
import astropy.units as u
import astropy.coordinates as coord
from astropy.io import ascii, fits, votable
from astropy.io.votable.exceptions import W46
from astropy.table import MaskedColumn, Table
import numpy as np
from requests import HTTPError, Response

# Local imports
from ._response_utils import response_looks_like_html, sanitize_votable_content
from ._utils import (
    _append_min_constraint,
    _append_range_constraint,
    _api_error_summary,
    _configured_token_from_env as _configured_token_from_env_base,
    _oauth_redirect_url,
    _strip_optional_quotes,
    _successful_response_error,
)
from ...query import BaseQuery
from ... import log
from ...utils import commons
from ...exceptions import InvalidQueryError, LoginError, RemoteServiceError, TableParseError
from . import conf


__all__ = ['Lamost', 'LamostClass']


_TOKEN_ENV_VARS = (
    "ASTROQUERY_LAMOST_TOKEN",
    "ASTROQUERY_NADC_LAMOST_TOKEN",
    "NADC_LAMOST_TOKEN",
    "CHINAVO_LAMOST_TOKEN",
    "ASTROQUERY_LAMOST_ACCESS_TOKEN",
    "ASTROQUERY_NADC_LAMOST_ACCESS_TOKEN",
    "NADC_LAMOST_ACCESS_TOKEN",
    "CHINAVO_LAMOST_ACCESS_TOKEN",
)


_SPECTRAL_FIELDS = {
    'low': {
        'snr': 'snrg',
        'teff': 'teff',
        'logg': 'logg',
        'feh': 'feh',
    },
    'medium': {
        'snr': 'snr',
        'teff': 'teff_lasp',
        'logg': 'logg_lasp',
        'feh': 'feh_lasp',
    },
}

_SCHEMA_COLUMN_KEYS = ('column_name', 'colname', 'column', 'name')
_SCHEMA_DATATYPE_KEYS = ('datatype', 'data_type', 'type', 'dbtype', 'dtype')
_SCHEMA_UNIT_KEYS = ('unit', 'units')


def _configured_token_from_env():
    return _configured_token_from_env_base(_TOKEN_ENV_VARS)


def _table_response_text(response):
    """Decode table text without ever interpreting an empty body as a path."""
    text = response.content.decode('utf-8-sig')
    if not text.strip():
        raise TableParseError(
            "LAMOST returned an empty response body, not a valid table. "
            "A zero-row text result must still contain a column header."
        )
    return text


def _csv_delimiter(text):
    """Use only unquoted delimiters in the header to recognize legacy pipes."""
    delimiters = set()
    quoted = False
    for character in text.lstrip():
        if character == '"':
            quoted = not quoted
        elif not quoted:
            if character in '\r\n':
                break
            if character in ',|':
                delimiters.add(character)
    if len(delimiters) > 1:
        raise ValueError("Ambiguous table header: both comma and pipe delimiters occur outside quotes.")
    return delimiters.pop() if delimiters else ','


class _CsvData(ascii.Csv.data_class):
    def process_lines(self, lines):
        # Inputs are complete CSV records. Whitespace-only records can be data;
        # physical blank lines inside quoted fields must never be filtered out.
        return lines


class _CsvInputter(ascii.BaseInputter):
    def get_lines(self, table, newline=None):
        # Even a single header record may contain quoted newlines. BaseInputter
        # would split that record again instead of preserving the logical record.
        return table


class _CsvReader(ascii.Csv):
    data_class = _CsvData
    inputter_class = _CsvInputter


class LamostClass(BaseQuery):
    """
    Class for querying the LAMOST spectroscopic survey database.

    The LAMOST (Large Sky Area Multi-Object Fiber Spectroscopic Telescope)
    survey provides spectroscopic data for millions of stars and galaxies.
    This class provides methods to query the catalog, download spectra,
    and access metadata.

    Notes
    -----
    Configuration is read when an instance is created, including the
    module-level ``Lamost`` instance at import time. Changing ``conf`` does
    not update existing instances. Create a new `LamostClass` to apply it.
    Authenticated requests and streaming downloads bypass the disk cache.
    ``get_query_payload=True`` returns parameters with credentials redacted.

    Authentication failures raise `~astroquery.exceptions.LoginError`.
    Other HTTP failures raise `requests.HTTPError`, service error payloads
    raise `~astroquery.exceptions.RemoteServiceError`, and malformed tables
    raise `~astroquery.exceptions.TableParseError`. Diagnostics retain
    available error details with credentials redacted.
    """

    URL = conf.server
    TIMEOUT = conf.timeout

    def __init__(
        self,
        *,
        token=None,
        data_release=None,
        sub_version=None,
        pylamost_config=None,
    ):
        """
        Initialize a LAMOST query instance.

        Parameters
        ----------
        token : str, optional
            Authentication token for LAMOST API access. Preferred ways to
            provide it are passing ``token`` directly, configuring
            ``astroquery.nadc.lamost.conf.token`` / astroquery.cfg, setting an
            environment variable such as ``ASTROQUERY_NADC_LAMOST_TOKEN``.
            An explicit empty string forces anonymous access.
        data_release : str, optional
            Data release version (e.g., 'dr10', 'dr11', 'dr12').
            Defaults to value from conf.data_release.
        sub_version : str, optional
            API sub-version (e.g., 'v2.0', 'v1.0').
            Defaults to value from conf.sub_version.
        pylamost_config : str or path-like, optional
            Explicit path to a pylamost-style config file. The file is only
            read when this argument is provided and no token was found from
            ``token``, ``conf.token``, or the environment.

        Notes
        -----
        Explicit pylamost-style config loading expects this format::

            token=your_token_here
            # Comments starting with # are supported

        Examples
        --------
        >>> from astroquery.nadc.lamost import LamostClass
        >>> lm = LamostClass(token='', data_release='dr10', sub_version='v2.0')
        >>> lm = LamostClass(pylamost_config='~/pylamost.ini')  # doctest: +SKIP
        """
        super().__init__()
        self.URL = conf.server.rstrip('/')
        self.TIMEOUT = conf.timeout
        explicit_token = token is not None
        if token is not None:
            self.token = _strip_optional_quotes(token) or None
        else:
            self.token = (
                _strip_optional_quotes(conf.token or "")
                or _configured_token_from_env()
            )
        self.data_release = data_release or conf.data_release
        self.sub_version = sub_version or conf.sub_version

        if self.token is None and not explicit_token and pylamost_config is not None:
            self._detect_token(pylamost_config)

    def _redact(self, value):
        """Copy diagnostic values without authentication credentials."""
        if isinstance(value, Mapping):
            return {
                key: '<redacted>' if str(key).lower() in {
                    'token', 'access_token', 'authorization', 'cookie', 'set-cookie',
                } else self._redact(item)
                for key, item in value.items()
            }
        if isinstance(value, (list, tuple)):
            return type(value)(self._redact(item) for item in value)
        if isinstance(value, bytes):
            return self._redact(value.decode('utf-8', 'replace')).encode('utf-8')
        if not isinstance(value, str):
            return value
        if self.token:
            for secret in {self.token, quote(self.token, safe=''), quote_plus(self.token)}:
                value = value.replace(secret, '<redacted>')
        return re.sub(
            r'''(?i)(\b(?:access_)?token["']?\s*[:=]\s*["']?)[^&\s"'<>]+''',
            r'\1<redacted>', value,
        )

    def _diagnostic_request(self, request):
        if request is None:
            return None
        sanitized = copy(request)
        sanitized.url = self._redact(request.url)
        sanitized.headers = self._redact(request.headers)
        sanitized.body = self._redact(request.body)
        sanitized._cookies = None
        return sanitized

    def _diagnostic_response(self, response):
        """Keep a safe response for parser diagnostics, without a live stream."""
        # copy(response) consumes streaming bodies through Response.__getstate__.
        # Redirect histories can also contain the response itself.
        sanitized = Response()
        sanitized.status_code = response.status_code
        sanitized.encoding = response.encoding
        sanitized.reason = self._redact(response.reason)
        sanitized.url = self._redact(getattr(response, 'url', None))
        sanitized.headers = self._redact(response.headers)
        sanitized.request = self._diagnostic_request(getattr(response, 'request', None))
        content = getattr(response, '_content', None)
        if isinstance(content, bytes):
            sanitized._content = self._redact(content)
            sanitized._content_consumed = True
        return sanitized

    def _sanitize_exception(self, error):
        """Redact attached requests and the accessible exception chain."""
        pending, seen = [error], set()
        while pending:
            current = pending.pop()
            if id(current) in seen:
                continue
            seen.add(id(current))
            current.args = tuple(self._redact(str(arg)) for arg in current.args)
            for name in ('url', 'filename', 'filename2', 'doc', 'msg', 'reason'):
                value = getattr(current, name, None)
                if isinstance(value, (str, Exception)):
                    setattr(current, name, self._redact(str(value)))
            if getattr(current, 'request', None) is not None:
                current.request = self._diagnostic_request(current.request)
            if getattr(current, 'response', None) is not None:
                current.response = self._diagnostic_response(current.response)
            pending.extend(item for item in (current.__cause__, current.__context__) if item is not None)

    def _response_hook(self, response, *args, **kwargs):
        # BaseQuery's hook logs raw URLs, headers and bodies, including tokens.
        log.debug('LAMOST HTTP %s %s: %s', response.request.method,
                  self._redact(response.request.url), response.status_code)

    def _safe_cache(self, cache, *, stream=False):
        """
        Return a cache flag safe to pass to ``BaseQuery._request``.

        Notes
        -----
        - When using authenticated requests (token provided), caching is disabled
          to avoid persisting credentials to disk via request/response caching.
        - When streaming downloads, caching is disabled to avoid pickling large
          responses and to prevent issues with partially-consumed streams.
        """
        if stream:
            return False
        if self.token:
            return False
        return cache

    def _normalize_resolution(self, resolution):
        normalized = str(resolution).strip().lower()
        if normalized not in {'low', 'medium'}:
            raise InvalidQueryError("resolution must be one of: low, medium.")
        return normalized

    def _normalize_output_format(self, output_format, *, allowed):
        normalized = str(output_format).strip().lower().lstrip('.')
        if normalized not in allowed:
            raise InvalidQueryError(
                "output_format must be one of: {0}.".format(", ".join(allowed))
            )
        return normalized

    def _request_raise(self, method, url, *, params=None, json=None,
                       timeout=None, cache=True, stream=False):
        response = None
        try:
            response = self._request(
                method,
                url,
                params=params,
                json=json,
                timeout=timeout or self.TIMEOUT,
                cache=self._safe_cache(cache, stream=stream),
                stream=stream,
            )
            context = f"LAMOST {method} {self._redact(url)} (HTTP {response.status_code})"
            oauth_redirect = _oauth_redirect_url(response)
            if oauth_redirect is not None or response.status_code in (401, 403):
                detail = ('the service redirected to an OAuth login page'
                          if oauth_redirect is not None else _api_error_summary(response))
                raise LoginError(
                    f"{context}: Authentication required: {self._redact(detail)}. "
                    "Pass token to LamostClass(token=...), or set conf.token / "
                    "ASTROQUERY_NADC_LAMOST_TOKEN before creating a new instance."
                )

            try:
                response.raise_for_status()
            except HTTPError as error:
                detail = _api_error_summary(response) or str(error)
                error.args = (f"{context}: {self._redact(detail)}",)
                raise

            error = _successful_response_error(response)
            if error:
                raise RemoteServiceError(f"{context}: {self._redact(error)}")
            return response
        except Exception as error:
            if response is None:
                response = getattr(error, 'response', None)
            if response is not None:
                self.response = self._diagnostic_response(response)
                response.close()
            self._sanitize_exception(error)
            raise

    @staticmethod
    def _radius_in(radius, unit):
        try:
            angle = coord.Angle(radius, unit=unit)
            value = angle.to_value(unit)
            if not angle.isscalar or not np.isfinite(value) or value <= 0:
                raise ValueError
        except (TypeError, ValueError, u.UnitsError) as error:
            raise InvalidQueryError(
                f"radius must be a finite positive scalar angle; bare numbers are in {unit}."
            ) from error
        return float(value)

    @staticmethod
    def _page_integer(value, name, *, minimum=1):
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral) or value < minimum:
            raise InvalidQueryError(f"{name} must be an integer >= {minimum}.")
        return int(value)

    def _parse_table_response(self, response, *, verbose=False, column_schema=None):
        table = self._parse_result(response, verbose=verbose, column_schema=column_schema)
        self.table = table
        return table

    def _response_json(self, response):
        self._validate_data_response(response)
        try:
            return response.json()
        except ValueError as error:
            self.response = self._diagnostic_response(response)
            self._sanitize_exception(error)
            raise TableParseError(f"Failed to parse LAMOST JSON response: {error}") from error

    def _validate_data_response(self, response):
        """Reject empty bodies and HTML before interpreting table or metadata formats."""
        if not response.content.removeprefix(b'\xef\xbb\xbf').strip():
            message = (
                "LAMOST returned an empty response body, not a valid table or metadata response. "
                "A zero-row text result must still contain a column header."
            )
        elif response_looks_like_html(response):
            message = (
                "Server returned HTML instead of a table response; "
                "check authentication or upstream errors."
            )
        else:
            return

        diagnostic = self.response = self._diagnostic_response(response)
        raise TableParseError(
            f"{message} HTTP {diagnostic.status_code}; URL: {diagnostic.url}; "
            f"Content-Type: {diagnostic.headers.get('Content-Type', 'unknown')}; "
            f"body: {len(response.content)} bytes."
        )

    def _normalize_tables_metadata(self, data):
        if isinstance(data, dict):
            if 'tables' in data and isinstance(data['tables'], dict):
                return data

            if 'tables' in data and isinstance(data['tables'], list):
                table_entries = data['tables']
            elif any(isinstance(value, dict) for value in data.values()):
                return {'tables': data}
            else:
                table_entries = [data]
        elif isinstance(data, list):
            table_entries = data
        else:
            raise TableParseError(
                "Expected a JSON object or list for table metadata."
            )

        tables = {}
        for entry in table_entries:
            if not isinstance(entry, dict):
                continue
            table_name = (
                entry.get('table_name')
                or entry.get('name')
                or entry.get('table')
                or entry.get('tablename')
            )
            if table_name is None:
                continue
            tables[str(table_name)] = dict(entry)

        return {'tables': tables}

    @staticmethod
    def _schema_value(metadata, keys):
        if not isinstance(metadata, Mapping):
            return None
        for key in keys:
            value = metadata.get(key)
            if value not in (None, ''):
                return value
        return None

    def _catalog_schema(self, catalog_name, *, cache=True):
        metadata = self.get_tables_metadata(cache=cache)
        tables = metadata.get('tables', {})
        if catalog_name not in tables:
            raise InvalidQueryError(f"Unknown LAMOST catalog: {catalog_name!r}.")

        table_metadata = tables[catalog_name]
        if not isinstance(table_metadata, Mapping):
            raise TableParseError(
                f"Metadata for LAMOST catalog {catalog_name!r} is not an object."
            )

        schema = self._columns_schema(table_metadata.get('columns'))
        if not schema:
            raise TableParseError(
                f"Metadata for LAMOST catalog {catalog_name!r} did not contain columns."
            )
        return schema

    def _columns_schema(self, columns):
        schema = {}
        if isinstance(columns, Mapping):
            for name, column_metadata in columns.items():
                if isinstance(column_metadata, Mapping):
                    schema[str(name)] = dict(column_metadata)
                else:
                    schema[str(name)] = {'datatype': column_metadata}
        elif isinstance(columns, (list, tuple)):
            for column_metadata in columns:
                if isinstance(column_metadata, str):
                    schema[column_metadata] = {}
                    continue
                name = self._schema_value(column_metadata, _SCHEMA_COLUMN_KEYS)
                if name is not None:
                    schema[str(name)] = dict(column_metadata)

        return schema

    def _validate_catalog_query(self, catalog_name, *, columns=None,
                                column_constraints=None,
                                position_constraints=None,
                                sort_by=None, cache=True):
        schema = self._catalog_schema(catalog_name, cache=cache)
        requested = [columns] if isinstance(columns, str) else list(columns or ())
        referenced = [str(name) for name in requested]

        constraints = column_constraints or ()
        if isinstance(constraints, Mapping):
            constraints = [constraints]
        for constraint in constraints:
            if not isinstance(constraint, Mapping):
                raise InvalidQueryError("Each column constraint must be a mapping.")
            column_name = constraint.get('column_name')
            if not column_name:
                raise InvalidQueryError(
                    "Each column constraint must specify `column_name`."
                )
            referenced.append(str(column_name))

        if position_constraints:
            referenced.extend(('ra', 'dec'))
        if sort_by:
            referenced.append(str(sort_by))

        unknown = list(dict.fromkeys(name for name in referenced if name not in schema))
        if unknown:
            raise InvalidQueryError(
                "Unknown LAMOST catalog column(s): {0}.".format(
                    ", ".join(unknown)
                )
            )
        return schema

    def _apply_catalog_schema(self, table, schema):
        integer_types = {
            'bigint', 'int', 'integer', 'long', 'short', 'smallint', 'tinyint', 'int32', 'int64',
        }
        floating_types = {'decimal', 'double', 'double precision', 'float', 'numeric', 'real', 'float32', 'float64'}
        text_types = {'char', 'varchar', 'text', 'string', 'unicodechar'}

        for name, column_metadata in schema.items():
            if name not in table.colnames:
                continue

            datatype = self._schema_value(column_metadata, _SCHEMA_DATATYPE_KEYS)
            normalized_type = str(datatype or '').lower().split('(', 1)[0].strip()
            caster = (
                int if normalized_type in integer_types
                else float if normalized_type in floating_types
                else str if normalized_type in text_types
                else None
            )
            if caster is not None:
                column = table[name]
                values = []
                mask = []
                try:
                    for value in column:
                        missing = value is None or np.ma.is_masked(value)
                        if isinstance(value, (str, bytes)):
                            missing = missing or not value or (caster is not str and not value.strip())
                        if isinstance(value, bytes):
                            value = value.decode('utf-8')
                        if (caster is int and not missing and isinstance(value, (float, np.floating))
                                and (not np.isfinite(value) or value != np.trunc(value))):
                            raise ValueError("Non-integral value in an integer column.")
                        values.append(caster(value) if not missing else ('' if caster is str else 0))
                        mask.append(missing)
                    converted = MaskedColumn(
                        values, mask=mask, name=name,
                        dtype=np.int64 if caster is int else np.float64 if caster is float else str,
                        unit=column.unit, description=column.description, meta=column.meta,
                    )
                except (TypeError, ValueError, OverflowError) as exc:
                    raise TableParseError(
                        f"Column {name!r} could not be converted to schema "
                        f"datatype {normalized_type!r}."
                    ) from exc
                table.replace_column(name, converted)

            unit = self._schema_value(column_metadata, _SCHEMA_UNIT_KEYS)
            if unit not in (None, '') and table[name].unit is None:
                try:
                    table[name].unit = u.Unit(str(unit))
                except ValueError:
                    pass

    def _prepare_catalog_result(self, table, *, catalog_name, columns):
        requested = [columns] if isinstance(columns, str) else list(columns or ())
        requested = [str(name) for name in requested]
        if len(table) == 0 and not table.colnames and requested:
            table = Table(names=requested)

        missing = [name for name in requested if name not in table.colnames]
        if missing:
            raise TableParseError(
                "LAMOST query result omitted requested column(s): {0}.".format(
                    ", ".join(missing)
                )
            )

        if requested:
            remaining = [name for name in table.colnames if name not in requested]
            table = table[requested + remaining]

        table.meta['catalog'] = catalog_name
        table.meta['data_release'] = self.data_release
        table.meta['sub_version'] = self.sub_version
        return table

    def _normalize_unique_id_result(self, data):
        if not isinstance(data, dict):
            raise TableParseError(
                "Expected a JSON object for unique-id lookup results."
            )

        normalized = dict(data)
        unique_id = normalized.get('unique_id', normalized.get('uid'))

        related_obsids = normalized.get('related_obsids')
        if related_obsids is None:
            related_obsids = []
        elif isinstance(related_obsids, (list, tuple, set)):
            related_obsids = list(related_obsids)
        else:
            related_obsids = [related_obsids]

        low_obsids = normalized.get(
            'related_obsids_low',
            normalized.get('obsid-low', normalized.get('obsid_low', [])),
        )
        medium_obsids = normalized.get(
            'related_obsids_medium',
            normalized.get('obsid-medium', normalized.get('obsid_medium', [])),
        )

        for key, value in (
            ('related_obsids_low', low_obsids),
            ('related_obsids_medium', medium_obsids),
        ):
            if isinstance(value, (list, tuple, set)):
                normalized[key] = list(value)
            elif value in (None, ''):
                normalized[key] = []
            else:
                normalized[key] = [value]

        if not related_obsids:
            merged_obsids = []
            for sequence in (
                normalized['related_obsids_low'],
                normalized['related_obsids_medium'],
            ):
                for obsid in sequence:
                    if obsid not in merged_obsids:
                        merged_obsids.append(obsid)
            related_obsids = merged_obsids

        if unique_id is not None:
            normalized['unique_id'] = unique_id
        normalized['related_obsids'] = related_obsids
        return normalized

    def _detect_token(self, config_file):
        """
        Load authentication token from an explicit pylamost-style config file.

        If the config file doesn't exist or doesn't contain a token, no error
        is raised to allow queries for public data.
        """
        config = self._get_config(config_file)
        if config and 'token' in config:
            token_value = config['token'].strip()
            if token_value:  # Only set if token is not empty
                self.token = token_value

    def _get_config(self, config_file):
        """
        Read configuration from an explicit pylamost-style config file.

        The config file format supports:
        - key=value pairs
        - Comments starting with #
        - Empty lines (ignored)

        Returns
        -------
        dict or None
            Configuration dictionary with key-value pairs, or None if file doesn't exist
            or cannot be read.

        Examples
        --------
        Example config file content:
            # LAMOST API Token
            token=your_token_here
        """
        config_file = os.path.expanduser(os.fspath(config_file))
        if not os.path.exists(config_file):
            return None

        config = {}
        try:
            with open(config_file, 'r') as fh:
                for line_num, line in enumerate(fh, 1):
                    line = line.strip()
                    # Skip empty lines and comments
                    if not line or line.startswith('#'):
                        continue
                    # Parse key=value pairs
                    if '=' in line:
                        key, value = line.split('=', 1)
                        config[key.strip()] = value.strip()
                    else:
                        # Skip malformed lines silently for compatibility
                        continue
        except Exception as e:
            # If config file cannot be read, return None to allow public data queries
            warnings.warn(
                f"Could not read config file {config_file}: {e}. "
                "Continuing without token (public data only).",
                UserWarning
            )
            return None

        return config if config else None

    def get_dr_versions(self):
        """
        Get available Data Release versions.

        Returns
        -------
        list of dict
            Service-provided release metadata, including ``dr_version``,
            ``sub_version``, and ``public_status``. Additional fields depend
            on the service response.

        Examples
        --------
        >>> from astroquery.nadc.lamost import Lamost
        >>> versions = Lamost.get_dr_versions()  # doctest: +SKIP
        >>> for v in versions:  # doctest: +SKIP
        ...     print(f"{v['dr_version']}/{v['sub_version']}: {v['public_status']}")  # doctest: +SKIP
        """
        url = f"{self.URL.rstrip('/')}/dr_versions"
        response = self._request_raise('GET', url)
        data = self._response_json(response)
        return data.get('versions', [])

    def get_unique_id_and_related_obsids(self, *, obsid=None, ra=None, dec=None, radius=None,
                                         get_query_payload=False, cache=True):
        """
        Get unique ID and related observation IDs for a target.

        This method finds all observations of the same astronomical object,
        which is useful for tracking multiple observations of the same target.

        Parameters
        ----------
        obsid : str or int, optional
            Observation ID to query. Either obsid OR (ra, dec, radius) must be provided.
        ra : float, optional
            Right ascension in degrees. Required if obsid is not provided.
        dec : float, optional
            Declination in degrees. Required if obsid is not provided.
        radius : str, float, or astropy.units.Quantity, optional
            Positive scalar search radius, required if ``obsid`` is absent.
            Bare numbers are degrees; angle strings and angular quantities are accepted.
        get_query_payload : bool, optional
            Return redacted parameters without executing the request.
        cache : bool, optional
            If True, cache the query result. Default is True.

        Returns
        -------
        dict
            Dictionary containing normalized keys:
            - unique_id: Unique identifier for the astronomical object
            - related_obsids: Combined list of all related observation IDs
            - related_obsids_low: Related low-resolution observation IDs, if available
            - related_obsids_medium: Related medium-resolution observation IDs, if available

        Raises
        ------
        ValueError
            If neither obsid nor (ra, dec, radius) are provided.

        Examples
        --------
        >>> from astroquery.nadc.lamost import Lamost
        >>> # Query by obsid
        >>> result = Lamost.get_unique_id_and_related_obsids(obsid='101001')  # doctest: +SKIP
        >>> print(result['unique_id'])  # doctest: +SKIP
        >>> print(result['related_obsids'])  # doctest: +SKIP

        >>> # Query by coordinates
        >>> result = Lamost.get_unique_id_and_related_obsids(  # doctest: +SKIP
        ...     ra=10.0, dec=40.0, radius=0.001  # doctest: +SKIP
        ... )  # doctest: +SKIP
        """
        if obsid is None and (ra is None or dec is None or radius is None):
            raise ValueError(
                "Either 'obsid' OR all of ('ra', 'dec', 'radius') must be provided"
            )

        request_payload = {}

        if obsid is not None:
            request_payload['obsid'] = str(obsid)

        if ra is not None:
            request_payload['ra'] = float(ra)

        if dec is not None:
            request_payload['dec'] = float(dec)

        if radius is not None:
            request_payload['radius'] = self._radius_in(radius, u.deg)

        if self.token:
            request_payload['token'] = self.token

        if get_query_payload:
            return self._redact(request_payload)

        # Build URL
        url = f"{self.URL}/{self.data_release}/{self.sub_version}/get_unique_id_and_related_obsids"
        response = self._request_raise('GET', url, params=request_payload, cache=cache)
        return self._normalize_unique_id_result(self._response_json(response))

    def get_query_result_count(self, sqlid, *, ismed=None, cache=True):
        """
        Get the total count of results for a query identified by sqlid.

        This provisional compatibility interface requires an externally
        created server-side ID. The current OpenAPI does not document an ID
        creation workflow. A complete live pagination workflow has not been
        validated; ``query_catalog`` has separate page-based retrieval.

        Parameters
        ----------
        sqlid : int
            Externally created SQL query ID. ``query_sql`` does not return one.
        ismed : bool, optional
            Deprecated pylamost-style parameter (unused). Included for
            compatibility only.
        cache : bool, optional
            If True, cache the query result. Default is True.

        Returns
        -------
        int
            Total number of rows in the query result.

        Examples
        --------
        >>> from astroquery.nadc.lamost import Lamost
        >>> count = Lamost.get_query_result_count(12345)  # doctest: +SKIP
        >>> print(f"Total results: {count}")  # doctest: +SKIP
        """
        request_payload = {'sqlid': int(sqlid)}

        if self.token:
            request_payload['token'] = self.token

        # Build URL
        url = f"{self.URL}/{self.data_release}/{self.sub_version}/get_query_result_count"
        response = self._request_raise('GET', url, params=request_payload, cache=cache)
        data = self._response_json(response)
        total = data.get('total') if isinstance(data, Mapping) else data
        try:
            return self._page_integer(total, 'count', minimum=0)
        except InvalidQueryError as error:
            raise RemoteServiceError(
                "LAMOST returned an invalid result count; expected a nonnegative integer."
            ) from error

    def get_query_result_by_page(self, sqlid, count, *, rows=10000, page=1,
                                 output_format='json', fmt=None, cache=True):
        """
        Get a specific page of query results.

        This provisional compatibility method requires an externally created
        ``sqlid``; see `get_query_result_count` for the service limitation.

        Parameters
        ----------
        sqlid : int
            SQL query ID.
        count : int
            Total number of results (from get_query_result_count).
        rows : int, optional
            Positive number of rows per page. Keep this value fixed across
            pages, including the last page. Default is 10000.
        page : int, optional
            Page number (1-indexed). Default is 1.
        output_format : str, optional
            Output format: 'json', 'csv', 'votable', or 'txt'. Default is 'json'.
        fmt : str, optional
            Deprecated alias for output_format.
        cache : bool, optional
            If True, cache the query result. Default is True.

        Returns
        -------
        `~astropy.table.Table` or list
            Query results for the requested page. Returns empty list/table if
            page number exceeds total pages.

        Raises
        ------
        astroquery.exceptions.InvalidQueryError
            ``rows`` or ``page`` is not a positive integer, or ``count`` is
            not a nonnegative integer.
        astroquery.exceptions.RemoteServiceError
            The page length disagrees with ``count`` and ``rows``. Partial
            results are not returned.

        Examples
        --------
        >>> from astroquery.nadc.lamost import Lamost
        >>> count = Lamost.get_query_result_count(12345)  # doctest: +SKIP
        >>> page1 = Lamost.get_query_result_by_page(12345, count, rows=1000, page=1)  # doctest: +SKIP
        """
        rows = self._page_integer(rows, 'rows')
        page = self._page_integer(page, 'page')
        count = self._page_integer(count, 'count', minimum=0)
        if fmt is not None:
            output_format = fmt
        output_format = self._normalize_output_format(
            output_format,
            allowed=('json', 'csv', 'votable', 'txt'),
        )

        if page > (count // rows + (1 if count % rows else 0)):
            # Page exceeds total pages
            if output_format == 'json':
                return []
            else:
                return Table()

        request_payload = {
            'sqlid': int(sqlid),
            'rows': int(rows),
            'page': int(page),
            'output.fmt': output_format
        }

        if self.token:
            request_payload['token'] = self.token

        # Build URL
        url = f"{self.URL}/{self.data_release}/{self.sub_version}/get_query_result"
        response = self._request_raise('GET', url, params=request_payload, cache=cache)

        # Parse based on format
        result = self._response_json(response) if output_format == 'json' else self._parse_result(response)
        expected = min(rows, count - (page - 1) * rows)
        if not isinstance(result, (list, Table)) or len(result) != expected:
            raise RemoteServiceError(
                f"LAMOST page {page} did not contain the expected {expected} rows; "
                "the result is incomplete or the reported count changed."
            )
        return result

    def get_query_result(self, sqlid, *, output_format='json', fmt=None,
                         page_size=10000, cache=True):
        """
        Get complete query results by automatically fetching all pages.

        This is a convenience method that handles pagination automatically
        for large result sets. It is provisional and requires an externally
        created ``sqlid``; see `get_query_result_count` for the service limitation.

        Parameters
        ----------
        sqlid : int
            SQL query ID.
        output_format : str, optional
            Output format: 'json', 'csv', 'votable', or 'txt'. Default is 'json'.
        fmt : str, optional
            Deprecated alias for output_format.
        page_size : int, optional
            Number of rows per page. Default is 10000.
        cache : bool, optional
            If True, cache the query results. Default is True.

        Returns
        -------
        `~astropy.table.Table` or list
            Results from every page of the supplied ``sqlid``. JSON returns
            a list of records; other formats return a table.

        Raises
        ------
        astroquery.exceptions.RemoteServiceError
            The service returns an invalid count or an unexpected page length.
        astroquery.exceptions.InvalidQueryError
            ``page_size`` is not a positive integer.

        Examples
        --------
        >>> from astroquery.nadc.lamost import Lamost
        >>> full_result = Lamost.get_query_result(12345)  # doctest: +SKIP
        >>> print(f"Retrieved {len(full_result)} rows")  # doctest: +SKIP
        """
        page_size = self._page_integer(page_size, 'page_size')
        if fmt is not None:
            output_format = fmt
        output_format = self._normalize_output_format(
            output_format,
            allowed=('json', 'csv', 'votable', 'txt'),
        )

        # Get total count
        count = self._page_integer(self.get_query_result_count(sqlid, cache=cache), 'count', minimum=0)

        # Calculate number of pages
        total_pages = count // page_size + (1 if count % page_size else 0)

        if output_format == 'json':
            # For JSON, accumulate list of dicts
            result = []
            for page in range(1, total_pages + 1):
                page_result = self.get_query_result_by_page(
                    sqlid, count, rows=page_size, page=page,
                    output_format=output_format, cache=cache
                )
                result.extend(page_result)
            return result
        else:
            # For table formats, stack tables
            from astropy.table import vstack
            tables = []
            for page in range(1, total_pages + 1):
                page_result = self.get_query_result_by_page(
                    sqlid, count, rows=page_size, page=page,
                    output_format=output_format, cache=cache
                )
                if len(page_result) > 0:
                    tables.append(page_result)

            if len(tables) == 0:
                return Table()
            elif len(tables) == 1:
                return tables[0]
            else:
                return vstack(tables)

    def download_query_result(self, sqlid, filename, *, output_format='csv',
                              fmt=None, page_size=10000, cache=True):
        """
        Download complete query results to a file.

        This provisional compatibility method requires an externally created
        ``sqlid``; see `get_query_result_count` for the service limitation.

        Parameters
        ----------
        sqlid : int
            SQL query ID.
        filename : str
            Path to save the results.
        output_format : str, optional
            Output format: 'csv' (default), 'json', 'votable', or 'txt'.
            Nonempty TXT exports contain tab-separated values and a column-name header.
        fmt : str, optional
            Deprecated alias for output_format.
        page_size : int, optional
            Number of rows per page for fetching. Default is 10000.
        cache : bool, optional
            If True, cache the query results. Default is True.

        Returns
        -------
        str
            Path to the saved file.

        Notes
        -----
        The complete result is fetched once in the required transport format.
        CSV, JSON, and TXT exports use JSON transport to preserve character
        values such as identifiers with leading zeros. VOTable exports use
        VOTable transport to retain its column metadata.
        For zero rows, JSON writes ``[]``; CSV/TXT have no column names to
        write because their source is an empty JSON record list.
        Output is written to a temporary file on the destination filesystem,
        then replaces the destination after writing succeeds. A query or write
        failure leaves an existing destination unchanged.

        Examples
        --------
        >>> from astroquery.nadc.lamost import Lamost
        >>> filepath = Lamost.download_query_result(12345, 'results.csv')  # doctest: +SKIP
        """
        import csv
        import json

        if fmt is not None:
            output_format = fmt
        output_format = self._normalize_output_format(
            output_format,
            allowed=('json', 'csv', 'votable', 'txt'),
        )

        query_format = 'votable' if output_format == 'votable' else 'json'
        results = self.get_query_result(
            sqlid, output_format=query_format, page_size=page_size, cache=cache
        )

        directory = os.path.dirname(os.path.abspath(filename))
        with tempfile.TemporaryDirectory(dir=directory) as temp_dir:
            temp_filename = os.path.join(temp_dir, 'result')
            if output_format == 'votable':
                results.write(temp_filename, format='votable')
            else:
                with open(temp_filename, 'w', newline='', encoding='utf-8') as f:
                    if output_format == 'csv':
                        if len(results) > 0:
                            writer = csv.DictWriter(f, fieldnames=results[0].keys())
                            writer.writeheader()
                            writer.writerows(results)
                    elif output_format == 'json':
                        json.dump(results, f, indent=2)
                    else:
                        Table(rows=results).write(f, format='ascii.tab')
            os.replace(temp_filename, filename)

        return filename

    def get_tables_metadata(self, *, cache=True):
        """
        Get metadata for all available tables in the data release.

        This method retrieves information about the database schema, including
        table names, column names, data types, descriptions, and units.

        Parameters
        ----------
        cache : bool, optional
            If True, cache the query result. Default is True.

        Returns
        -------
        dict
            Dictionary containing metadata for all tables under ``result['tables']``.

        Examples
        --------
        >>> from astroquery.nadc.lamost import Lamost
        >>> metadata = Lamost.get_tables_metadata()  # doctest: +SKIP
        >>> for table_name, table_info in metadata['tables'].items():  # doctest: +SKIP
        ...     print(f"Table: {table_name}")  # doctest: +SKIP
        ...     print(f"  Columns: {table_info.get('columns', [])}")  # doctest: +SKIP
        """
        request_params = {}

        if self.token:
            request_params['token'] = self.token

        # Build URL
        url = f"{self.URL}/{self.data_release}/{self.sub_version}/tables"
        response = self._request_raise('GET', url, params=request_params, cache=cache)
        return self._normalize_tables_metadata(self._response_json(response))

    def get_tap_url(self, *, cache=True):
        """
        Get the IVOA TAP (Table Access Protocol) service URL.

        This URL can be used with TAP clients for advanced queries following
        the IVOA TAP standard.

        Parameters
        ----------
        cache : bool, optional
            If True, cache the query result. Default is True.

        Returns
        -------
        dict
            Dictionary containing the TAP service URL and related information.

        Examples
        --------
        >>> from astroquery.nadc.lamost import Lamost
        >>> tap_info = Lamost.get_tap_url()  # doctest: +SKIP
        >>> print(tap_info['tap_url'])  # doctest: +SKIP
        """
        request_params = {}

        if self.token:
            request_params['token'] = self.token

        # Build URL - TAP doesn't use resolution path
        url = f"{self.URL}/{self.data_release}/{self.sub_version}/voservice/tap_url"
        response = self._request_raise('GET', url, params=request_params, cache=cache)
        return self._response_json(response)

    def get_footprint(self, *, resolution='low', cache=True):
        """
        Get the LAMOST observation footprint image.

        This provides a visual representation of the sky coverage for
        the specified resolution.

        Parameters
        ----------
        resolution : str, optional
            Spectral resolution: 'low' for LRS (default) or 'medium' for MRS.
        cache : bool, optional
            If True, cache the image. Default is True.

        Returns
        -------
        bytes
            Footprint image data (typically PNG or JPEG format).

        Examples
        --------
        >>> from astroquery.nadc.lamost import Lamost
        >>> footprint_img = Lamost.get_footprint(resolution='low')  # doctest: +SKIP
        >>> # Save to file
        >>> with open('footprint.png', 'wb') as f:  # doctest: +SKIP
        ...     f.write(footprint_img)  # doctest: +SKIP
        """
        resolution = self._normalize_resolution(resolution)
        request_params = {}

        if self.token:
            request_params['token'] = self.token

        # Build URL
        url = self._build_url('footprint', resolution=resolution)
        response = self._request_raise('GET', url, params=request_params, cache=cache)
        return response.content

    def _build_url(self, endpoint, resolution='low'):
        """
        Build complete API URL for a given endpoint.

        Parameters
        ----------
        endpoint : str
            API endpoint path (e.g., 'voservice/conesearch', 'spectrum/fits')
        resolution : str, optional
            Spectral resolution: 'low' for LRS or 'medium' for MRS.
            Default is 'low'.

        Returns
        -------
        str
            Complete URL for the API request
        """
        resolution = self._normalize_resolution(resolution)
        res_path = 'mrs' if resolution == 'medium' else 'lrs'
        return f"{self.URL}/{self.data_release}/{self.sub_version}/{res_path}/{endpoint}"

    def _request_query_region(self, coordinates, radius, *,
                              resolution='low', output_format='csv', fmt=None,
                              get_query_payload=False, cache=True):
        """
        Query LAMOST catalog using cone search around given coordinates.
        """
        if fmt is not None:
            import warnings
            warnings.warn(
                "The 'fmt' parameter is deprecated. Use 'output_format' instead.",
                DeprecationWarning,
                stacklevel=2
            )
            output_format = fmt
        resolution = self._normalize_resolution(resolution)
        output_format = self._normalize_output_format(
            output_format,
            allowed=('votable', 'json', 'csv'),
        )

        c = commons.parse_coordinates(coordinates)

        request_payload = {
            'ra': c.icrs.ra.deg,
            'dec': c.icrs.dec.deg,
            'sr': self._radius_in(radius, u.deg),
            'output.fmt': output_format
        }

        if self.token:
            request_payload['token'] = self.token

        if get_query_payload:
            return self._redact(request_payload)

        url = self._build_url('voservice/conesearch', resolution=resolution)
        return self._request_raise('GET', url, params=request_payload, cache=cache)

    def query_region(self, coordinates, radius, *,
                     resolution='low', output_format='csv', fmt=None,
                     get_query_payload=False, cache=True, verbose=False):
        """Query the LAMOST catalog around sky coordinates.

        Parameters
        ----------
        coordinates : str or astropy.coordinates object
            Cone-search center.
        radius : str, float, or astropy.units.Quantity
            Positive scalar cone-search radius. Bare numbers are degrees;
            angle strings such as ``'5 arcsec'`` and angular quantities are accepted.
        resolution : {"low", "medium"}, optional
            Spectral-resolution service to query.
        output_format : {"votable", "json", "csv"}, optional
            Output format requested from the service. The default is ``'csv'``
            because some service VOTables declare string
            fields too short to represent the returned identifiers.
        fmt : str, optional
            Deprecated alias for ``output_format``.
        get_query_payload : bool, optional
            Return the request payload instead of executing the request.
        cache : bool, optional
            Whether to use astroquery's request cache.
        verbose : bool, optional
            Emit parser diagnostics for table responses.

        Returns
        -------
        astropy.table.Table or dict
            Query result table with catalog types and units from metadata,
            or redacted parameters when ``get_query_payload=True``. Table
            metadata identifies the catalog, release, and sub-version.

        Raises
        ------
        astroquery.exceptions.InvalidQueryError
            The radius is invalid or the requested format is unsupported.
        astroquery.exceptions.TableParseError
            The response is malformed or would truncate string values.

        Notes
        -----
        This method retrieves one response. CSV avoids the known VOTable
        defect but does not guarantee completeness for every release or size.
        """
        response = self._request_query_region(
            coordinates,
            radius,
            resolution=resolution,
            output_format=output_format,
            fmt=fmt,
            get_query_payload=get_query_payload,
            cache=cache,
        )
        if get_query_payload:
            return response
        catalog_name = self._spectral_catalog_name(resolution, None)
        schema = self._catalog_schema(catalog_name, cache=cache)
        table = self._parse_table_response(response, verbose=verbose, column_schema=schema)
        return self._prepare_catalog_result(table, catalog_name=catalog_name, columns=None)

    def _request_query_ssap(self, coordinates, radius, *,
                            resolution='low', output_format='csv', fmt=None,
                            get_query_payload=False, cache=True):
        """Query LAMOST using IVOA Simple Spectral Access Protocol (SSAP)."""
        if fmt is not None:
            import warnings
            warnings.warn(
                "The 'fmt' parameter is deprecated. Use 'output_format' instead.",
                DeprecationWarning,
                stacklevel=2
            )
            output_format = fmt
        resolution = self._normalize_resolution(resolution)
        output_format = self._normalize_output_format(
            output_format,
            allowed=('votable', 'json', 'csv'),
        )

        c = commons.parse_coordinates(coordinates)

        request_payload = {
            'pos': f"{c.icrs.ra.deg},{c.icrs.dec.deg}",
            'size': self._radius_in(radius, u.deg),
            'output.fmt': output_format
        }

        if self.token:
            request_payload['token'] = self.token

        if get_query_payload:
            return self._redact(request_payload)

        url = self._build_url('voservice/ssap', resolution=resolution)
        return self._request_raise('GET', url, params=request_payload, cache=cache)

    def query_ssap(self, coordinates, radius, *,
                   resolution='low', output_format='csv', fmt=None,
                   get_query_payload=False, cache=True, verbose=False):
        """Query LAMOST with the IVOA Simple Spectral Access Protocol.

        Parameters
        ----------
        coordinates : str or astropy.coordinates object
            Search center.
        radius : str, float, or astropy.units.Quantity
            Positive scalar search radius. Bare numbers are degrees;
            angle strings and angular quantities are accepted.
        resolution : {"low", "medium"}, optional
            Spectral-resolution service to query.
        output_format : {"votable", "json", "csv"}, optional
            Output format requested from the service. Default is ``'csv'``.
        fmt : str, optional
            Deprecated alias for ``output_format``.
        get_query_payload : bool, optional
            Return the request payload instead of executing the request.
        cache : bool, optional
            Whether to use astroquery's request cache.
        verbose : bool, optional
            Emit parser diagnostics for table responses.

        Returns
        -------
        astropy.table.Table or dict
            SSAP result table, or request payload when
            ``get_query_payload=True``.
        """
        response = self._request_query_ssap(
            coordinates,
            radius,
            resolution=resolution,
            output_format=output_format,
            fmt=fmt,
            get_query_payload=get_query_payload,
            cache=cache,
        )
        if get_query_payload:
            return response
        return self._parse_table_response(response, verbose=verbose)

    def _request_query_sql(self, sql, *, output_format='json',
                           get_query_payload=False, cache=True):
        """Execute a raw SQL query on the LAMOST database."""
        output_format = self._normalize_output_format(
            output_format,
            allowed=('json', 'csv', 'votable', 'txt'),
        )
        request_payload = {
            'sql': sql,
            'output.fmt': output_format
        }

        if self.token:
            request_payload['token'] = self.token

        if get_query_payload:
            return self._redact(request_payload)

        url = f"{self.URL}/{self.data_release}/{self.sub_version}/sql"
        return self._request_raise('GET', url, params=request_payload, cache=cache)

    def query_sql(self, sql, *, output_format='json', column_schema=None,
                  get_query_payload=False, cache=True, verbose=False):
        """Execute a raw SQL query on the LAMOST database.

        Parameters
        ----------
        sql : str
            SQL statement accepted by the LAMOST service.
        output_format : {"json", "csv", "votable", "txt"}, optional
            Output format requested from the service.
        column_schema : dict, optional
            Result column names mapped to metadata dictionaries, for example
            ``{'temperature': {'datatype': 'double', 'unit': 'K'}}``.
            Use this for SQL aliases and expressions without response column
            metadata. Unannotated JSON values retain their service types.
            With a schema, CSV/TXT fields are read as strings before conversion;
            columns without a declared datatype remain strings. Character
            identifiers retain leading zeros. Only declared units are attached.
        get_query_payload : bool, optional
            Return the GET parameters with credentials redacted instead of
            executing the request. ``column_schema`` is local to the parser.
        cache : bool, optional
            Whether to use astroquery's request cache.
        verbose : bool, optional
            Emit parser diagnostics for table responses.

        Returns
        -------
        astropy.table.Table or dict
            SQL result table, or request payload when
            ``get_query_payload=True``.

        Raises
        ------
        astroquery.exceptions.InvalidQueryError
            ``column_schema`` is not a mapping of column metadata dictionaries.
        astroquery.exceptions.TableParseError
            A declared column is absent or cannot be converted to its datatype.
        """
        if column_schema is not None and (
            not isinstance(column_schema, Mapping)
            or any(not isinstance(value, Mapping) for value in column_schema.values())
        ):
            raise InvalidQueryError("column_schema must map result column names to metadata dictionaries.")
        response = self._request_query_sql(
            sql,
            output_format=output_format,
            get_query_payload=get_query_payload,
            cache=cache,
        )
        if get_query_payload:
            return response
        table = self._parse_table_response(response, verbose=verbose, column_schema=column_schema)
        missing = set(column_schema or ()) - set(table.colnames)
        if missing:
            raise TableParseError(f"LAMOST SQL result omitted schema column(s): {', '.join(sorted(missing))}.")
        return table

    def _request_query_catalog(self, catalog_name, *,
                               column_constraints=None,
                               position_constraints=None,
                               columns=None,
                               sort_by=None,
                               sort_order='asc',
                               max_rows=100,
                               page=1,
                               output_format='json',
                               get_query_payload=False,
                               cache=True):
        """Advanced parametric table query with constraints."""
        output_format = self._normalize_output_format(
            output_format,
            allowed=('json', 'csv', 'votable', 'txt'),
        )
        request_payload = {
            'rows': max_rows,
            'page': page,
            'output.fmt': output_format,
            'order': sort_order
        }

        if column_constraints:
            request_payload['column_constraints'] = column_constraints

        if position_constraints:
            request_payload['pos'] = position_constraints
            request_payload['pos_group'] = 'ra,dec'

        if columns:
            request_payload['showcol'] = columns

        if sort_by:
            request_payload['sort'] = sort_by

        if get_query_payload:
            if self.token:
                return self._redact({
                    'json': request_payload,
                    'params': {'token': self.token},
                })
            return self._redact(request_payload)

        url = f"{self.URL}/{self.data_release}/{self.sub_version}/query/{catalog_name}"
        return self._request_raise(
            'POST',
            url,
            json=request_payload,
            params={'token': self.token} if self.token else None,
            cache=cache,
        )

    def query_catalog(self, catalog_name, *,
                      column_constraints=None,
                      position_constraints=None,
                      columns=None,
                      sort_by=None,
                      sort_order='asc',
                      max_rows=100,
                      page=1,
                      output_format='json',
                      get_query_payload=False,
                      cache=True,
                      verbose=False):
        """Query a LAMOST catalog with structured constraints.

        Parameters
        ----------
        catalog_name : str
            LAMOST catalog endpoint name.
        column_constraints : list of dict, optional
            Service-format column constraints.
        position_constraints : dict, optional
            Service-format spatial constraint. Executed cone queries require
            ``cone_nearestonly=True`` because this endpoint cannot guarantee
            all matches. Use ``query_region(..., output_format='csv')`` for
            ordinary field searches.
        columns : iterable of str, optional
            Columns to include in the result.
        sort_by : str, optional
            Column used for sorting.
        sort_order : {"asc", "desc"}, optional
            Sort order.
        max_rows : int, optional
            Maximum rows in this page. Additional pages are not fetched automatically.
        page : int, optional
            One-based page number.
        output_format : {"json", "csv", "votable", "txt"}, optional
            Output format requested from the service.
        get_query_payload : bool, optional
            Return the request payload instead of executing the request.
        cache : bool, optional
            Whether to use astroquery's request cache.
        verbose : bool, optional
            Emit parser diagnostics for table responses.

        Returns
        -------
        astropy.table.Table or dict
            Query result table, or request payload when
            ``get_query_payload=True``.

        Raises
        ------
        astroquery.exceptions.RemoteServiceError
            A cone query requesting all matches is executed.
        astroquery.exceptions.InvalidQueryError
            Catalog names or columns disagree with the service metadata.
        astroquery.exceptions.TableParseError
            Requested columns are missing or values disagree with their datatypes.
        """
        if columns is not None:
            columns = [columns] if isinstance(columns, str) else list(columns)
        if column_constraints is not None and not isinstance(column_constraints, Mapping):
            column_constraints = list(column_constraints)

        schema = None
        if not get_query_payload:
            cone = (
                position_constraints.get('cone')
                if isinstance(position_constraints, Mapping)
                else None
            )
            if isinstance(cone, Mapping) and cone.get('cone_nearestonly') is not True:
                raise RemoteServiceError(
                    "The LAMOST structured cone-query service cannot guarantee all matches "
                    "when nearest_only=False. For a field query use "
                    "query_region(..., output_format='csv'); use nearest_only=True "
                    "only when a single nearest match is intended."
                )
            schema = self._validate_catalog_query(
                catalog_name,
                columns=columns,
                column_constraints=column_constraints,
                position_constraints=position_constraints,
                sort_by=sort_by,
                cache=cache,
            )

        response = self._request_query_catalog(
            catalog_name,
            column_constraints=column_constraints,
            position_constraints=position_constraints,
            columns=columns,
            sort_by=sort_by,
            sort_order=sort_order,
            max_rows=max_rows,
            page=page,
            output_format=output_format,
            get_query_payload=get_query_payload,
            cache=cache,
        )
        if get_query_payload:
            return response
        result_schema = schema if not columns else {name: schema[name] for name in columns}
        table = self._parse_result(response, verbose=verbose, column_schema=result_schema)
        table = self._prepare_catalog_result(
            table,
            catalog_name=catalog_name,
            columns=columns,
        )
        self.table = table
        return table

    def _build_structured_cone_constraint(self, coordinates, radius, *, nearest_only=False):
        c = commons.parse_coordinates(coordinates)

        return {
            'cone': {
                'racenter': float(c.icrs.ra.to_value(u.deg)),
                'deccenter': float(c.icrs.dec.to_value(u.deg)),
                'radius': self._radius_in(radius, u.arcsec),
                'cone_nearestonly': bool(nearest_only),
            },
        }

    def _spectral_catalog_name(self, resolution, catalog_name):
        resolution = self._normalize_resolution(resolution)
        if catalog_name is not None:
            return catalog_name
        return 'med_combined' if resolution == 'medium' else 'combined'

    def _spectral_fields(self, resolution):
        resolution = self._normalize_resolution(resolution)
        return _SPECTRAL_FIELDS[resolution]

    def _spectral_column_constraints(
        self,
        *,
        resolution='low',
        snr_min=None,
        snr_column=None,
        teff_range=None,
        logg_range=None,
        feh_range=None,
    ):
        fields = self._spectral_fields(resolution)
        if snr_column is None:
            snr_column = fields['snr']

        constraints = []
        if snr_min is not None:
            if not snr_column:
                raise InvalidQueryError("snr_column must be provided when snr_min is set.")
            _append_min_constraint(constraints, snr_column, snr_min)
        _append_range_constraint(constraints, fields['teff'], teff_range)
        _append_range_constraint(constraints, fields['logg'], logg_range)
        _append_range_constraint(constraints, fields['feh'], feh_range)
        return constraints or None

    def query_spectra(
        self,
        coordinates=None,
        radius=None,
        *,
        resolution='low',
        catalog_name=None,
        snr_min=None,
        snr_column=None,
        teff_range=None,
        logg_range=None,
        feh_range=None,
        columns=None,
        nearest_only=False,
        sort_by=None,
        sort_order='asc',
        max_rows=100,
        page=1,
        output_format='json',
        get_query_payload=False,
        cache=True,
        verbose=False,
    ):
        """Query LAMOST spectral catalog rows with common quality filters.

        Parameters
        ----------
        coordinates : str or astropy.coordinates object, optional
            Center position for a cone search.
        radius : str, float, or astropy.units.Quantity, optional
            Positive scalar cone-search radius, required with ``coordinates``.
            Bare numbers are arcseconds; angle strings and angular quantities are accepted.
        resolution : {"low", "medium"}, optional
            Spectral-resolution catalog family.
        catalog_name : str, optional
            Explicit catalog endpoint name. Field-name defaults still follow
            ``resolution``.
        snr_min : float, optional
            Minimum signal-to-noise ratio.
        snr_column : str, optional
            Column used for the SNR filter. Defaults to ``snrg`` for LRS and
            ``snr`` for MRS.
        teff_range, logg_range, feh_range : sequence of float, optional
            Inclusive stellar-parameter ranges: temperature in kelvin, log10
            surface gravity in cm/s2, and [Fe/H] in dex, respectively.
        columns : iterable of str, optional
            Columns to include in the result.
        nearest_only : bool, optional
            Return only the nearest spatial match. With coordinates and
            ``False``, execution raises `~astroquery.exceptions.RemoteServiceError`
            because the service cannot guarantee all matches. Payload inspection
            remains available; use ``query_region(..., output_format='csv')``
            for ordinary field searches.
        sort_by : str, optional
            Column used for sorting.
        sort_order : {"asc", "desc"}, optional
            Sort order.
        max_rows : int, optional
            Maximum rows per page.
        page : int, optional
            One-based page number.
        output_format : {"json", "csv", "votable", "txt"}, optional
            Output format requested from the service.
        get_query_payload : bool, optional
            Return the request payload instead of executing the request.
        cache : bool, optional
            Whether to use astroquery's request cache.
        verbose : bool, optional
            Emit parser diagnostics for table responses.

        Returns
        -------
        astropy.table.Table or dict
            Query result table, or request payload when
            ``get_query_payload=True``.
        """
        if (coordinates is None) ^ (radius is None):
            raise InvalidQueryError("coordinates and radius must be provided together.")

        position_constraints = None
        if coordinates is not None:
            position_constraints = self._build_structured_cone_constraint(
                coordinates,
                radius,
                nearest_only=nearest_only,
            )

        return self.query_catalog(
            self._spectral_catalog_name(resolution, catalog_name),
            column_constraints=self._spectral_column_constraints(
                resolution=resolution,
                snr_min=snr_min,
                snr_column=snr_column,
                teff_range=teff_range,
                logg_range=logg_range,
                feh_range=feh_range,
            ),
            position_constraints=position_constraints,
            columns=columns,
            sort_by=sort_by,
            sort_order=sort_order,
            max_rows=max_rows,
            page=page,
            output_format=output_format,
            get_query_payload=get_query_payload,
            cache=cache,
            verbose=verbose,
        )

    def query_stellar_parameters(
        self,
        coordinates=None,
        radius=None,
        *,
        resolution='low',
        catalog_name=None,
        snr_min=None,
        snr_column=None,
        teff_range=None,
        logg_range=None,
        feh_range=None,
        columns=None,
        nearest_only=False,
        sort_by=None,
        sort_order='asc',
        max_rows=100,
        page=1,
        output_format='json',
        get_query_payload=False,
        cache=True,
        verbose=False,
    ):
        """Query LAMOST rows focused on stellar atmospheric parameters.

        Parameters
        ----------
        coordinates : str or astropy.coordinates object, optional
            Center position for a cone search.
        radius : str, float, or astropy.units.Quantity, optional
            Positive scalar cone-search radius, required with ``coordinates``.
            Bare numbers are arcseconds; angle strings and angular quantities are accepted.
        resolution : {"low", "medium"}, optional
            Spectral-resolution catalog family.
        catalog_name : str, optional
            Explicit catalog endpoint name. Field-name defaults still follow
            ``resolution``.
        snr_min : float, optional
            Minimum signal-to-noise ratio.
        snr_column : str, optional
            Column used for the SNR filter and included by default. Defaults
            to ``snrg`` for LRS and ``snr`` for MRS.
        teff_range, logg_range, feh_range : sequence of float, optional
            Inclusive stellar-parameter ranges: temperature in kelvin, log10
            surface gravity in cm/s2, and [Fe/H] in dex, respectively.
        columns : iterable of str, optional
            Columns to include in the result. Defaults to core stellar
            atmospheric-parameter columns.
        nearest_only : bool, optional
            Return only the nearest spatial match. With coordinates and
            ``False``, execution raises `~astroquery.exceptions.RemoteServiceError`
            because the service cannot guarantee all matches. Payload inspection
            remains available; use ``query_region(..., output_format='csv')``
            for ordinary field searches.
        sort_by : str, optional
            Column used for sorting.
        sort_order : {"asc", "desc"}, optional
            Sort order.
        max_rows : int, optional
            Maximum rows per page.
        page : int, optional
            One-based page number.
        output_format : {"json", "csv", "votable", "txt"}, optional
            Output format requested from the service.
        get_query_payload : bool, optional
            Return the request payload instead of executing the request.
        cache : bool, optional
            Whether to use astroquery's request cache.
        verbose : bool, optional
            Emit parser diagnostics for table responses.

        Returns
        -------
        astropy.table.Table or dict
            Query result table, or request payload when
            ``get_query_payload=True``.
        """
        fields = self._spectral_fields(resolution)
        if snr_column is None:
            snr_column = fields['snr']

        if columns is None:
            columns = [
                'obsid',
                'ra',
                'dec',
                fields['teff'],
                fields['logg'],
                fields['feh'],
            ]
            if snr_column:
                columns = [*columns, snr_column]

        return self.query_spectra(
            coordinates=coordinates,
            radius=radius,
            resolution=resolution,
            catalog_name=catalog_name,
            snr_min=snr_min,
            snr_column=snr_column,
            teff_range=teff_range,
            logg_range=logg_range,
            feh_range=feh_range,
            columns=columns,
            nearest_only=nearest_only,
            sort_by=sort_by,
            sort_order=sort_order,
            max_rows=max_rows,
            page=page,
            output_format=output_format,
            get_query_payload=get_query_payload,
            cache=cache,
            verbose=verbose,
        )

    def query_repeat_observations(
        self,
        *,
        obsid=None,
        coordinates=None,
        radius=None,
        ra=None,
        dec=None,
        get_query_payload=False,
        cache=True,
    ):
        """Query related LAMOST observation IDs for one target.

        Parameters
        ----------
        obsid : str or int, optional
            Observation ID to resolve.
        coordinates : str or astropy.coordinates object, optional
            Target position. Mutually exclusive with ``ra`` and ``dec``.
        radius : str, float, or astropy.units.Quantity, optional
            Positive scalar search radius. Bare numbers are degrees;
            angle strings and angular quantities are accepted.
        ra, dec : float, optional
            Target coordinates in degrees.
        get_query_payload : bool, optional
            Return the request payload instead of executing the request.
        cache : bool, optional
            Whether to use astroquery's request cache.

        Returns
        -------
        dict
            Related-observation payload, or request payload when
            ``get_query_payload=True``.
        """
        has_coordinate_values = coordinates is not None or ra is not None or dec is not None or radius is not None
        if obsid is not None and has_coordinate_values:
            raise InvalidQueryError("Provide either obsid or coordinates/radius, not both.")

        if coordinates is not None:
            if ra is not None or dec is not None:
                raise InvalidQueryError("Use coordinates or ra/dec, not both.")
            c = commons.parse_coordinates(coordinates)
            ra = float(c.icrs.ra.to_value(u.deg))
            dec = float(c.icrs.dec.to_value(u.deg))

        return self.get_unique_id_and_related_obsids(
            obsid=obsid,
            ra=ra,
            dec=dec,
            radius=radius,
            get_query_payload=get_query_payload,
            cache=cache,
        )

    def _request_metadata(self, obsid, *,
                          resolution='low', get_query_payload=False, cache=True):
        """Get metadata/information for a specific spectrum."""
        resolution = self._normalize_resolution(resolution)
        request_payload = {'obsid': str(obsid)}

        if self.token:
            request_payload['token'] = self.token

        if get_query_payload:
            return self._redact(request_payload)

        url = self._build_url('spectrum/info', resolution=resolution)
        return self._request_raise('GET', url, params=request_payload, cache=cache)

    def get_metadata(self, obsid, *, resolution='low',
                     get_query_payload=False, cache=True, verbose=False):
        """Get metadata for a specific LAMOST spectrum.

        Parameters
        ----------
        obsid : str or int
            Observation ID of the spectrum.
        resolution : {"low", "medium"}, optional
            Spectral-resolution service to query.
        get_query_payload : bool, optional
            Return the request payload instead of executing the request.
        cache : bool, optional
            Whether to use astroquery's request cache.
        verbose : bool, optional
            Emit parser diagnostics for table responses.

        Returns
        -------
        astropy.table.Table or dict
            Metadata table, or request payload when
            ``get_query_payload=True``.
        """
        response = self._request_metadata(
            obsid,
            resolution=resolution,
            get_query_payload=get_query_payload,
            cache=cache,
        )
        if get_query_payload:
            return response
        return self._parse_table_response(response, verbose=verbose)

    def get_spectra(self, obsid, *, resolution='low', get_query_payload=False,
                    verify='warn'):
        """Download spectrum FITS data for an observation ID.

        Parameters
        ----------
        obsid : str or int
            Observation ID of the spectrum.
        resolution : {"low", "medium"}, optional
            Spectral-resolution service to query.
        get_query_payload : bool, optional
            Return redacted request parameters instead of fetching files.
        verify : str, optional
            FITS verification option passed to the ``verify`` method of the
            `astropy.io.fits.HDUList`. Default is ``"warn"``.

        Returns
        -------
        list of astropy.io.fits.HDUList or dict
            In-memory FITS objects, or redacted parameters when
            ``get_query_payload=True``. The caller must close the HDU lists.
            HTTP responses are closed before returning.
        """
        url_list = self.get_spectrum_list(obsid, resolution=resolution,
                                          get_query_payload=get_query_payload)
        if get_query_payload:
            return url_list

        with self._request_raise('GET', url_list[0]) as response:
            with commons.get_readable_fileobj(BytesIO(response.content), encoding='binary') as source:
                spectrum = fits.HDUList.fromstring(source.read())
        try:
            spectrum.verify(verify)
        except Exception:
            spectrum.close()
            raise
        return [spectrum]

    def get_spectrum_list(self, obsid, *, resolution='low', get_query_payload=False):
        """
        Get list of spectrum FITS file URLs without downloading.

        Parameters
        ----------
        obsid : str or int
            Observation ID of the spectrum.
        resolution : str, optional
            Spectral resolution: 'low' for LRS (default) or 'medium' for MRS.
        get_query_payload : bool, optional
            If True, return the request parameters dict without executing the query.

        Returns
        -------
        list of str or dict
            List of FITS file URLs, or redacted parameters when
            ``get_query_payload=True``. Download URLs contain the token for
            authenticated clients and must not be logged or shared.
        """
        resolution = self._normalize_resolution(resolution)
        request_payload = {'obsid': str(obsid)}

        if self.token:
            request_payload['token'] = self.token

        if get_query_payload:
            return self._redact(request_payload)

        # Build URL for FITS download
        from urllib.parse import urlencode
        base_url = self._build_url('spectrum/fits', resolution=resolution)
        url = f"{base_url}?{urlencode(request_payload)}"

        # Return as single-item list for consistency with get_images pattern
        return [url]

    def get_fits_csv(self, obsid, *, resolution='low', ismed=None, cache=True):
        """
        Get spectrum data in CSV format.

        This method retrieves the spectrum data as CSV text, which is useful
        for quick preview and data processing without downloading FITS files.

        Parameters
        ----------
        obsid : str or int
            Observation ID of the spectrum.
        resolution : str, optional
            Spectral resolution: 'low' for LRS (default) or 'medium' for MRS.
        ismed : bool, optional
            Deprecated pylamost-style resolution flag. If True, uses MRS.
        cache : bool, optional
            If True, cache the query result. Default is True.

        Returns
        -------
        str
            CSV-formatted spectrum data as text.

        Examples
        --------
        >>> from astroquery.nadc.lamost import Lamost
        >>> csv_data = Lamost.get_fits_csv('101001', resolution='low')  # doctest: +SKIP
        >>> print(csv_data[:100])  # doctest: +SKIP
        """
        if ismed is not None:
            resolution = 'medium' if ismed else 'low'
        resolution = self._normalize_resolution(resolution)

        request_payload = {'obsid': str(obsid)}

        if self.token:
            request_payload['token'] = self.token

        # Build URL for CSV export
        url = self._build_url('spectrum/fits2csv', resolution=resolution)
        response = self._request_raise('GET', url, params=request_payload, cache=cache)
        return response.text

    def get_images(self, obsid, *, resolution='low', get_query_payload=False):
        """Download spectrum PNG images for an observation ID.

        Parameters
        ----------
        obsid : str or int
            Observation ID of the spectrum.
        resolution : {"low", "medium"}, optional
            Spectral-resolution service to query.
        get_query_payload : bool, optional
            Return redacted request parameters instead of fetching files.

        Returns
        -------
        list of bytes or dict
            PNG image bytes, or redacted parameters when
            ``get_query_payload=True``. HTTP responses are closed before returning.
        """
        url_list = self.get_image_list(obsid, resolution=resolution,
                                       get_query_payload=get_query_payload)
        if get_query_payload:
            return url_list

        with self._request_raise('GET', url_list[0]) as response:
            return [response.content]

    def get_image_list(self, obsid, *, resolution='low', get_query_payload=False):
        """
        Get list of spectrum PNG image URLs without downloading.

        Parameters
        ----------
        obsid : str or int
            Observation ID of the spectrum.
        resolution : str, optional
            Spectral resolution: 'low' for LRS (default) or 'medium' for MRS.
        get_query_payload : bool, optional
            If True, return the request parameters dict without executing the query.

        Returns
        -------
        list of str or dict
            List of PNG image URLs, or redacted parameters when
            ``get_query_payload=True``. Download URLs contain the token for
            authenticated clients and must not be logged or shared.
        """
        resolution = self._normalize_resolution(resolution)
        request_payload = {'obsid': str(obsid)}

        if self.token:
            request_payload['token'] = self.token

        if get_query_payload:
            return self._redact(request_payload)

        # Build URL for PNG download
        from urllib.parse import urlencode
        base_url = self._build_url('spectrum/png', resolution=resolution)
        url = f"{base_url}?{urlencode(request_payload)}"

        # Return as single-item list
        return [url]

    def _resolve_download_path(self, response, save_dir, filename, default_name):
        import os
        import re

        if filename:
            if os.path.isabs(filename):
                return filename
            return os.path.join(save_dir, filename)

        content_disp = response.headers.get('Content-Disposition', '')
        filename_match = re.search(r'filename=([^;]+)', content_disp, re.IGNORECASE)
        if filename_match:
            resolved = filename_match.group(1).strip(' "\'')
            resolved = os.path.basename(resolved)
        else:
            resolved = default_name

        return os.path.join(save_dir, resolved)

    def download_catalog(self, catalog_name, *, resolution='low',
                         save_dir='./', savedir=None, ismed=None,
                         filename=None, overwrite=True, cache=True,
                         verify='exception'):
        """
        Download a catalog FITS product.

        Parameters
        ----------
        catalog_name : str
            Name of the catalog to download.
        resolution : str, optional
            Spectral resolution: 'low' for LRS (default) or 'medium' for MRS.
        save_dir : str, optional
            Directory to save the downloaded file. Default is current directory.
        savedir : str, optional
            Deprecated alias for save_dir.
        ismed : bool, optional
            Deprecated pylamost-style resolution flag. If True, uses MRS.
        filename : str, optional
            Override the output filename. Relative paths are resolved against
            save_dir. If None, use server-provided filename and fall back to
            "{catalog_name}.fits.gz".
        overwrite : bool, optional
            If True, overwrite existing file. Default is True.
        cache : bool, optional
            Retained for compatibility. Streaming downloads bypass the cache.
        verify : str, optional
            FITS verification option passed to
            the ``verify`` method of `astropy.io.fits.HDUList`. The default, ``"exception"``,
            rejects non-compliant catalog products before publication.

        Returns
        -------
        str
            Path to the downloaded catalog file, or the existing file when
            ``overwrite=False``. Existing files are not revalidated in this case.

        Notes
        -----
        The response is closed on success, skip, and failure. A temporary file
        with a unique name is validated before replacing the destination;
        a failed write or FITS
        verification preserves an existing destination.

        Examples
        --------
        >>> from astroquery.nadc.lamost import Lamost
        >>> filepath = Lamost.download_catalog('catalog_v1', resolution='low')  # doctest: +SKIP
        """
        if savedir is not None:
            save_dir = savedir

        if ismed is not None:
            resolution = 'medium' if ismed else 'low'
        resolution = self._normalize_resolution(resolution)

        request_params = {'name': catalog_name}

        if self.token:
            request_params['token'] = self.token

        # Build URL
        url = self._build_url('catalog', resolution=resolution)

        # Download file
        with self._request_raise(
            'GET',
            url,
            params=request_params,
            cache=cache,
            stream=True,
        ) as response:
            filepath = self._resolve_download_path(
                response, save_dir, filename, f"{catalog_name}.fits.gz"
            )
            os.makedirs(os.path.dirname(filepath) or save_dir, exist_ok=True)
            if os.path.exists(filepath) and not overwrite:
                return filepath

            try:
                directory = os.path.dirname(os.path.abspath(filepath))
                with tempfile.TemporaryDirectory(dir=directory) as temp_dir:
                    temp_filepath = os.path.join(temp_dir, 'catalog')
                    with open(temp_filepath, 'wb') as f:
                        for chunk in response.iter_content(chunk_size=8192):
                            if chunk:
                                f.write(chunk)
                    with fits.open(temp_filepath) as hdul:
                        hdul.verify(verify)
                    os.replace(temp_filepath, filepath)
            except Exception as error:
                self._sanitize_exception(error)
                raise

        return filepath

    def _parse_result(self, response, *, verbose=False, column_schema=None):
        """
        Parse response into an astropy Table based on content type.

        Supports votable, json, csv, and txt formats based on Content-Type header.

        Parameters
        ----------
        response : `requests.Response`
            HTTP response from query.
        verbose : bool, optional
            If False, suppress VOTable warnings. Default is False.

        Returns
        -------
        table : `~astropy.table.Table`
            Parsed table from the response.
        """
        try:
            self._validate_data_response(response)
            error = _successful_response_error(response)
            if error:
                raise RemoteServiceError(f"LAMOST API error: {self._redact(error)}")

            content_type = response.headers.get('Content-Type', '').lower()
            if 'json' in content_type:
                table = self._parse_json_result(response)
            elif 'csv' in content_type:
                table = self._parse_csv_result(response, column_schema=column_schema)
            elif 'text/plain' in content_type:
                table = self._parse_txt_result(response, column_schema=column_schema)
            else:
                table = self._parse_votable_result(response, verbose=verbose)

            schema = self._columns_schema(table.meta.get('columns'))
            schema.update(column_schema or {})
            if not len(table) and not table.colnames and schema:
                table = Table(names=list(schema), meta=table.meta)
            self._apply_catalog_schema(table, schema)
            return table
        except Exception as error:
            self.response = self._diagnostic_response(response)
            self._sanitize_exception(error)
            raise

    def _parse_votable_result(self, response, *, verbose=False):
        """
        Parse VOTable response into an astropy Table.

        Parameters
        ----------
        response : `requests.Response`
            HTTP response from query.
        verbose : bool, optional
            If False, suppress VOTable warnings. Default is False.

        Returns
        -------
        table : `~astropy.table.Table`
            Parsed table from the response.
        """
        try:
            self._validate_data_response(response)

            content = sanitize_votable_content(
                response.content,
                fix_invalid_date=True,
                fix_missing_field_datatype=True,
                fix_empty_arraysize=True,
            )
            tf = BytesIO(content)
            with warnings.catch_warnings(record=True) as caught_warnings:
                warnings.simplefilter('always')
                votable_file = votable.parse(tf, verify='warn')
                first_table = votable_file.get_first_table()
                table = first_table.to_table(use_names_over_ids=True)

            if any(issubclass(item.category, W46) for item in caught_warnings):
                raise TableParseError(
                    "LAMOST VOTable declares a string arraysize shorter than "
                    "the returned value; refusing to return truncated data."
                )
            if verbose:
                for item in caught_warnings:
                    warnings.warn(item.message)
            return table
        except TableParseError as ex:
            self.response = response
            self.table_parse_error = ex
            raise
        except Exception as ex:
            # Store response for debugging
            self.response = response
            self.table_parse_error = ex
            raise TableParseError(
                f"Failed to parse response as VOTable: {str(ex)}"
            )

    def _parse_json_result(self, response):
        """
        Parse JSON response into an astropy Table.

        Parameters
        ----------
        response : `requests.Response`
            HTTP response from query with JSON content.

        Returns
        -------
        table : `~astropy.table.Table`
            Parsed table from the JSON response.
        """
        try:
            data = self._response_json(response)

            # Handle different JSON structures
            if isinstance(data, list):
                # List of dictionaries
                if len(data) == 0:
                    return Table()
                return Table(data)
            elif isinstance(data, dict):
                # Single result or structured response
                if 'columns' in data and 'rows' in data:
                    columns = data['columns']
                    rows = data['rows']
                    if not isinstance(rows, list):
                        raise ValueError("JSON response 'rows' must be a list.")

                    column_names = []
                    if isinstance(columns, (list, tuple)):
                        for column in columns:
                            if isinstance(column, str):
                                column_names.append(column)
                            else:
                                name = self._schema_value(column, _SCHEMA_COLUMN_KEYS)
                                if name is not None:
                                    column_names.append(str(name))

                    if rows and all(isinstance(row, Mapping) for row in rows):
                        table = Table(rows)
                    elif rows and column_names:
                        table = Table(rows=rows, names=column_names)
                    elif rows:
                        table = Table(rows=rows)
                    else:
                        table = Table(names=column_names)

                    if column_names:
                        ordered = [name for name in column_names if name in table.colnames]
                        remaining = [name for name in table.colnames if name not in ordered]
                        table = table[ordered + remaining]
                    table.meta['columns'] = columns
                    if 'total' in data:
                        table.meta['total'] = data['total']
                    return table
                elif 'data' in data:
                    # Structured response with 'data' field
                    return Table(data['data'])
                else:
                    # Single result
                    return Table([data])
            else:
                raise ValueError(f"Unexpected JSON structure: {type(data)}")

        except Exception as ex:
            self.response = response
            self.json_parse_error = ex
            raise TableParseError(
                f"Failed to parse response as JSON: {str(ex)}"
            )

    def _parse_csv_result(self, response, *, column_schema=None):
        """
        Parse CSV or unambiguous historical pipe-delimited table text.

        Parameters
        ----------
        response : `requests.Response`
            HTTP response from query with CSV content.

        Returns
        -------
        table : `~astropy.table.Table`
            Parsed table from the CSV response.
        """
        try:
            text = _table_response_text(response)
            delimiter = _csv_delimiter(text)
            # Astropy's CSV reader pads short rows and renames duplicate headers.
            # Validate the structure first so damaged responses cannot look valid.
            # StringIO splits only on CSV line endings, unlike str.splitlines(),
            # which also splits characters such as vertical tabs inside fields.
            lines = list(StringIO(text, newline=''))
            records = csv.reader(lines, delimiter=delimiter, strict=True)
            header = None
            logical_records = []
            previous_line = 0
            for record in records:
                original = ''.join(lines[previous_line:records.line_num])
                previous_line = records.line_num
                if not record or (header is None and len(record) == 1 and not record[0].strip()):
                    continue
                if header is None:
                    header = [name.strip() for name in record]
                    if any(not name for name in header) or len(set(header)) != len(header):
                        raise ValueError("Table header contains empty or duplicate column names.")
                elif len(record) != len(header):
                    raise ValueError(
                        f"Table record ending at line {records.line_num} has {len(record)} fields; "
                        f"expected {len(header)}."
                    )
                logical_records.append(original)
            if header is None:
                raise ValueError("Table response has no column header.")
            reader = _CsvReader()
            reader.names = header
            reader.header.splitter.delimiter = delimiter
            reader.header.splitter.process_line = None
            reader.data.splitter.delimiter = delimiter
            reader.data.splitter.process_line = None
            reader.data.splitter.process_val = None
            reader.data.splitter.skipinitialspace = False
            # Match ascii.read's empty-field masks and numeric inference, while
            # retaining strings until an explicit catalog schema is applied.
            reader.data.fill_values = [('', '0')]
            reader.outputter.converters = {'*': [ascii.convert_numpy(str)]} if column_schema else {}
            return reader.read(logical_records)
        except Exception as ex:
            self.response = response
            raise TableParseError(
                f"Failed to parse response as CSV: {str(ex)}"
            )

    def _parse_txt_result(self, response, *, column_schema=None):
        """
        Parse plain text response into an astropy Table.

        Parameters
        ----------
        response : `requests.Response`
            HTTP response from query with plain text content.

        Returns
        -------
        table : `~astropy.table.Table`
            Parsed table from the text response.
        """
        try:
            text = _table_response_text(response)
            return ascii.read(
                text.splitlines(keepends=True),
                format='tab',
                comment=r'\s*#',
                guess=False,
                fast_reader=False,
                converters={'*': [ascii.convert_numpy(str)]} if column_schema else {},
            )
        except Exception as ex:
            self.response = response
            raise TableParseError(
                f"Failed to parse response as plain text: {str(ex)}"
            )


# Singleton instance for module-level access
Lamost: LamostClass = LamostClass()


# Utility functions for FITS spectrum processing
