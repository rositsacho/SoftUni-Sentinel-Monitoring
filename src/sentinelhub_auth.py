"""
sentinelhub_auth.py - OAuth2 authentication for the Copernicus Data Space
Ecosystem (CDSE) Sentinel Hub.

This module exchanges a CDSE OAuth2 *client* (client_id + client_secret) for
a short-lived Bearer access token using the OAuth2 "client_credentials"
grant, as described by the CDSE documentation:

    https://documentation.dataspace.copernicus.eu/APIs/SentinelHub/Overview/Authentication.html

Credentials are read locally from the project's .env file with python-dotenv
(keys: CDSE_CLIENT_ID, CDSE_CLIENT_SECRET). The secret and the resulting
access token are NEVER printed, logged or written to disk by this module.
Only non-sensitive diagnostics (HTTP status codes, expiry) are surfaced.
"""

from __future__ import annotations

import os
import time
from dataclasses import dataclass
from pathlib import Path

import requests
from dotenv import load_dotenv

# CDSE identity provider (Keycloak realm "CDSE").
TOKEN_URL = (
    "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/"
    "protocol/openid-connect/token"
)

# Project root = one level above this file (src/ -> project root).
PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ENV_PATH = PROJECT_ROOT / ".env"

# Refresh the token a little before it actually expires so an in-flight
# request never uses a token that expires mid-call.
EXPIRY_SAFETY_MARGIN_S = 30

# Network timeout for the token request (seconds).
DEFAULT_TIMEOUT_S = 30


class SentinelAuthError(RuntimeError):
    """Raised when CDSE OAuth2 authentication fails."""


class MissingCredentialsError(SentinelAuthError):
    """Raised when CDSE_CLIENT_ID / CDSE_CLIENT_SECRET are not available."""


def load_credentials(env_path: str | Path | None = None) -> tuple[str, str]:
    """Return (client_id, client_secret) loaded from .env / environment.

    Raises MissingCredentialsError if either variable is absent or blank.
    The returned values are never logged by this function.
    """
    path = Path(env_path) if env_path is not None else DEFAULT_ENV_PATH
    if path.exists():
        load_dotenv(path)
    else:
        # Fall back to a .env in the current working directory, then to the
        # plain process environment. Loading twice is harmless.
        load_dotenv()

    client_id = (os.getenv("CDSE_CLIENT_ID") or "").strip()
    client_secret = (os.getenv("CDSE_CLIENT_SECRET") or "").strip()

    missing = [
        name
        for name, value in (
            ("CDSE_CLIENT_ID", client_id),
            ("CDSE_CLIENT_SECRET", client_secret),
        )
        if not value
    ]
    if missing:
        raise MissingCredentialsError(
            "Missing CDSE credential(s): "
            + ", ".join(missing)
            + ". Add them to the project .env file (see .env.example / README)."
        )
    return client_id, client_secret


def _safe_error_text(response: requests.Response) -> str:
    """Extract a short, non-sensitive error description from a response.

    OAuth error bodies contain an error code/description - never the secret -
    but we still keep the text short to avoid flooding logs.
    """
    try:
        payload = response.json()
    except ValueError:
        return (response.text or "").strip()[:300]

    if isinstance(payload, dict):
        for key in ("error_description", "message", "error", "reason"):
            value = payload.get(key)
            if value:
                return str(value)[:300]
        return str(payload)[:300]
    return str(payload)[:300]


@dataclass
class _Token:
    """In-memory access token with its (monotonic) expiry instant."""

    access_token: str
    expires_at: float

    def is_valid(self) -> bool:
        return time.monotonic() < self.expires_at


class SentinelHubAuth:
    """Caches a CDSE access token and refreshes it on demand.

    Usage:
        auth = SentinelHubAuth()
        headers = auth.auth_headers()   # {'Authorization': 'Bearer ...'}
    """

    def __init__(
        self,
        env_path: str | Path | None = None,
        timeout: float = DEFAULT_TIMEOUT_S,
        session: requests.Session | None = None,
    ) -> None:
        self._env_path = env_path
        self._timeout = timeout
        self._session = session or requests.Session()
        self._token: _Token | None = None
        self._expires_in_s: float | None = None

    @property
    def token_expires_in_s(self) -> float | None:
        """Seconds of validity reported by the last successful token call."""
        return self._expires_in_s

    def get_token(self) -> str:
        """Return a valid Bearer token, refreshing it when necessary."""
        if self._token is None or not self._token.is_valid():
            self._token = self._request_token()
        return self._token.access_token

    def auth_headers(self) -> dict[str, str]:
        """Return request headers carrying the Bearer access token."""
        return {"Authorization": f"Bearer {self.get_token()}"}

    def _request_token(self) -> _Token:
        client_id, client_secret = load_credentials(self._env_path)
        payload = {
            "grant_type": "client_credentials",
            "client_id": client_id,
            "client_secret": client_secret,
        }
        try:
            response = self._session.post(
                TOKEN_URL, data=payload, timeout=self._timeout
            )
        except requests.RequestException as exc:
            # exc string never contains the secret; safe to surface the type.
            raise SentinelAuthError(
                f"Could not reach the CDSE token endpoint: {exc}"
            ) from exc

        if response.status_code != 200:
            raise SentinelAuthError(
                f"CDSE authentication failed (HTTP {response.status_code}): "
                f"{_safe_error_text(response)}"
            )

        try:
            data = response.json()
        except ValueError as exc:
            raise SentinelAuthError(
                "CDSE token endpoint returned a non-JSON response."
            ) from exc

        token = data.get("access_token")
        if not token:
            raise SentinelAuthError(
                "CDSE token response did not contain an 'access_token' field."
            )

        try:
            expires_in = float(data.get("expires_in", 600))
        except (TypeError, ValueError):
            expires_in = 600.0
        self._expires_in_s = expires_in

        return _Token(
            access_token=token,
            expires_at=time.monotonic()
            + max(expires_in - EXPIRY_SAFETY_MARGIN_S, 1.0),
        )


if __name__ == "__main__":
    # Local, offline sanity check: verifies the module imports and that the
    # credential lookup behaves as expected WITHOUT contacting CDSE.
    import sys

    try:
        cid, _secret = load_credentials()
        print(
            "Credentials found. CDSE_CLIENT_ID length: "
            f"{len(cid)} chars (value not shown)."
        )
    except MissingCredentialsError as exc:
        print(f"Credentials not configured: {exc}")
        sys.exit(1)

    auth = SentinelHubAuth()
    print(
        "Auth module ready. Call auth.auth_headers() to obtain a token "
        "(network request, not performed here)."
    )
