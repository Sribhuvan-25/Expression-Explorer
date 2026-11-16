"""
Test-wide guard: no test may open a network connection (CLAUDE.md:
"Don't let tests hit the network").

Why this exists: `GET /datasets` loads every registered dataset that is
not warming, and swallows load failures by design (METHODS.md 7.9). So a
test calling it on a machine with network access silently downloaded
every real cohort -- GDC TARGET, DepMap, GEO -- which was ~1,150 s of a
~1,230 s CI run, with nothing failing to say so. With connections refused
here, those loaders fail in milliseconds and the endpoint's behaviour
under test is unchanged; and any new test that genuinely tries to reach
the network fails loudly instead of quietly slowing CI.

Loopback stays allowed (the FastAPI TestClient is in-process anyway, and
a local Postgres in docker-compose listens on localhost). A test that
deliberately exercises an unreachable host (the database fall-back test)
gets the same refusal it expects.
"""
import ipaddress
import socket

import pytest

_LOOPBACK_NAMES = {"localhost", "localhost.localdomain"}


class NetworkBlockedInTests(OSError):
    pass


def _is_loopback(address) -> bool:
    host = address[0] if isinstance(address, tuple) else address
    if isinstance(host, bytes):
        host = host.decode()
    if not isinstance(host, str):
        return False
    if host in _LOOPBACK_NAMES:
        return True
    try:
        return ipaddress.ip_address(host).is_loopback
    except ValueError:
        return False


@pytest.fixture(autouse=True)
def _block_network(monkeypatch):
    real_connect = socket.socket.connect
    real_connect_ex = socket.socket.connect_ex
    real_create_connection = socket.create_connection
    real_getaddrinfo = socket.getaddrinfo

    def guarded_connect(self, address):
        if self.family in (socket.AF_INET, socket.AF_INET6) and not _is_loopback(address):
            raise NetworkBlockedInTests(f"Network access blocked in tests: {address!r}")
        return real_connect(self, address)

    def guarded_connect_ex(self, address):
        if self.family in (socket.AF_INET, socket.AF_INET6) and not _is_loopback(address):
            raise NetworkBlockedInTests(f"Network access blocked in tests: {address!r}")
        return real_connect_ex(self, address)

    def guarded_create_connection(address, *args, **kwargs):
        if not _is_loopback(address):
            raise NetworkBlockedInTests(f"Network access blocked in tests: {address!r}")
        return real_create_connection(address, *args, **kwargs)

    def guarded_getaddrinfo(host, *args, **kwargs):
        # Refusing name resolution too stops a client from spending its
        # DNS timeout before the connect is refused.
        if host is not None and not _is_loopback(host):
            raise NetworkBlockedInTests(f"Network access blocked in tests: {host!r}")
        return real_getaddrinfo(host, *args, **kwargs)

    monkeypatch.setattr(socket.socket, "connect", guarded_connect)
    monkeypatch.setattr(socket.socket, "connect_ex", guarded_connect_ex)
    monkeypatch.setattr(socket, "create_connection", guarded_create_connection)
    monkeypatch.setattr(socket, "getaddrinfo", guarded_getaddrinfo)
