"""Retry failed downloads via `urllib.request.urlopen`.

Remote files used in the example notebooks, e.g., PEtab problems from
cdn.jsdelivr.net, are intermittently unavailable (HTTP 403, 429, 5xx, dropped
connections). Failed jsDelivr GitHub downloads with an explicit branch or
commit (`/gh/<owner>/<repo>@<ref>/...`) are retried from
raw.githubusercontent.com, other failed HTTP(S) downloads are retried after a
delay.
This module is loaded at interpreter startup if its directory is on
`PYTHONPATH` (see `scripts/runNotebook.sh` and `doc/conf.py`). It covers
`urlopen`, `urlretrieve`, and pandas/PEtab reading from URLs.
"""

import http.client
import re
import sys
import time
import urllib.error
import urllib.request

_RETRY_HTTP_CODES = {403, 429, 500, 502, 503, 504}
_MAX_ATTEMPTS = 5
# Without a ref, or with a version (range), jsDelivr serves a release, whose
# GitHub ref may differ, so these are not rewritten.
_JSDELIVR_GH = re.compile(
    r"https://cdn\.jsdelivr\.net/gh/([^/@]+)/([^/@]+)@([^/]+)/(.*)"
)

_urlopen = urllib.request.urlopen


def _github_url(url):
    """The raw.githubusercontent.com URL for a jsDelivr GitHub URL, or None."""
    if isinstance(url, urllib.request.Request):
        if (new_url := _github_url(url.full_url)) is None:
            return None
        # A new request, as a failed attempt adds headers such as `Host`
        # to the original one
        return urllib.request.Request(
            new_url,
            data=url.data,
            headers=url.headers,
            method=url.get_method(),
        )

    if not (m := _JSDELIVR_GH.fullmatch(url)):
        return None
    owner, repo, ref, path = m.groups()
    if ref == "latest" or re.fullmatch(r"v?\d+(\.\d+)*", ref):
        return None
    return f"https://raw.githubusercontent.com/{owner}/{repo}/{ref}/{path}"


def _urlopen_with_retries(url, *args, **kwargs):
    full_url = url.full_url if isinstance(url, urllib.request.Request) else url
    if not full_url.startswith(("http://", "https://")):
        return _urlopen(url, *args, **kwargs)

    for attempt in range(1, _MAX_ATTEMPTS + 1):
        try:
            return _urlopen(url, *args, **kwargs)
        except urllib.error.HTTPError as e:
            if e.code not in _RETRY_HTTP_CODES or attempt == _MAX_ATTEMPTS:
                raise
            reason = f"HTTP {e.code}"
            e.close()
        except (OSError, http.client.HTTPException) as e:
            # `URLError`s, but also, e.g., a dropped connection or a timeout
            # while waiting for the response, which `urlopen` doesn't wrap
            if attempt == _MAX_ATTEMPTS:
                raise
            reason = getattr(e, "reason", None) or repr(e)
        if (github_url := _github_url(url)) is not None:
            print(
                f"Download failed ({reason}), retrying from GitHub.",
                file=sys.stderr,
            )
            url = github_url
            continue
        delay = 2**attempt
        print(
            f"Download failed ({reason}), retrying in {delay} s.",
            file=sys.stderr,
        )
        time.sleep(delay)


urllib.request.urlopen = _urlopen_with_retries
